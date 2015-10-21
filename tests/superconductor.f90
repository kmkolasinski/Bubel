! ------------------------------------------------------ !
! Quantulaba - superconductor.f90 - Krzysztof Kolasinski 2015
!
! Example based on Kwant tutorial:
! http://kwant-project.org/doc/1.0/tutorial/tutorial5
! ------------------------------------------------------ !
program supercond
use modunits     ! unit conversion tools
use modscatter   ! eigen values and transport
use modshape
implicit none
character(*),parameter :: output_folder = "superconductor_output/"
type(qscatter)             :: qt
type(qshape)               :: lead_shape
doubleprecision            :: zeros_vector(200)
doubleprecision            :: a_Emin,a_Emax
integer                    :: no_expected_states
integer ,parameter         :: nx = 30
integer ,parameter         :: ny = 10
doubleprecision,parameter  :: dx = 1.0 ! [nm]
integer , dimension(nx,ny) :: gindex ! converts local index (i,j) to global index
integer                    :: i,j,k,m
doubleprecision            :: lead_translation_vec(3),Ef
character(300) :: line
 ! Initalize system
call qt%init_system()
! ----------------------------------------------------------
! 1. Create mesh - loop over width and height of the lattice
! ----------------------------------------------------------
k      = 0
gindex = 0
do i = 1 , nx
do j = 1 , ny
    ! Initalize atom structure with position of the atom.
    call qt%qatom%init((/ (i-1) * dx , (j-1) * dx + (i-1) * dx * 0.0 , 0.0 * dx /),no_in_states=2)
    ! Add atom to the system.
    call qt%qsystem%add_atom(qt%qatom)
    k           = k + 1
    gindex(i,j) = k
enddo
enddo
! ----------------------------------------------------------
! 2. Construct logical connections between sites on the mesh.
! ----------------------------------------------------------
qt%qnnbparam%box = (/2*dx,2*dx,0.0D0/) ! do not search for the sites far than (-dx:+dx) direction
call qt%qsystem%make_lattice(qt%qnnbparam,c_matrix=coupling)


! ----------------------------------------------------------
! 3. Find eigenvalues of the system.
! You must provide range of energy (Emin,Emax) and
! number of expected states in that region. Eigen problem
! is solved with FEAST library.
! ----------------------------------------------------------
no_expected_states = 100
! we choose small enery window
a_Emin =0.0
a_Emax =0.2
call qt%qsystem%calc_eigenproblem(a_Emin,a_Emax,no_expected_states)
! Write founded eigenstates to file if there are states in the range (Emin,Emax)
if(qt%qsystem%no_eigenvalues > 0) then
zeros_vector = 0
open(unit=222,file=output_folder//"eigenvecs.dat")
do i = 1 , nx
do j = 1 , ny
    if(qt%qsystem%atoms(gindex(i,j))%bActive)then
    ! get unique ID of current site
    k = qt%qsystem%atoms(gindex(i,j))%globalIDs(2)  ! take spin = 2
    m = qt%qsystem%atoms(gindex(i,j))%globalIDs(1)  ! take spin = 1
    ! save total density
    write(222,"(500e20.6)"),qt%qsystem%atoms(gindex(i,j))%atom_pos,abs(qt%qsystem%eigenvecs(k,:))**2+abs(qt%qsystem%eigenvecs(m,:))**2
    else
    ! fill with zeros empty spaces
    write(222,"(500e20.6)"),qt%qsystem%atoms(gindex(i,j))%atom_pos,zeros_vector(1:qt%qsystem%no_eigenvalues)
    endif
enddo
    write(222,*),"" ! GNUPLOT necessary line break
enddo
close(222)
endif ! end of if are there eigenstates?
! Those can be viewed in Viewer
call qt%qsystem%save_data(output_folder//"eigenvectors.xml",array2d=abs(transpose(qt%qsystem%eigenvecs))**2)

! ----------------------------------------------------------
! 4. Use generated mesh to calculate the band structure
! in the region of homogenous lead.
! ----------------------------------------------------------
lead_translation_vec = (/  dx , 0.0D0 , 0.0D0 /)
call lead_shape%init_range_3d((/ -0.5*dx , 0.0D0 , 0.0D0 /),lead_translation_vec)
call qt%add_lead(lead_shape,lead_translation_vec)

lead_translation_vec = (/ -dx , 0.0D0 , 0.0D0 /)
call lead_shape%init_range_3d((/ (nx+0.5-1)*dx , 0.0D0 , 0.0D0 /),lead_translation_vec)
call qt%add_lead(lead_shape,lead_translation_vec)

a_Emin =-10.0
a_Emax = 10.0
call qt%leads(1)%bands(output_folder//"bands.dat",-M_PI,+M_PI,M_PI/160.0,a_Emin,a_Emax)


! Save lattice to file to see if constructed system is OK!
call qt%save_system(output_folder//"system.xml")



Ef = 1.0
QSYS_DEBUG_LEVEL = 1 ! tell me more
call qt%calculate_modes(Ef)
call qt%solve(1,Ef)

! Save calculated electron density to file
do i = 1 , size(qt%qauxvec)
    qt%qauxvec(i) = sum(qt%densities(:,i))
enddo
call qt%qsystem%save_data(output_folder//"densities.xml",array2d=qt%densities,array1d=qt%qauxvec)


print*,"Performing energy scan..."
open(unit=111,file=output_folder//"T.dat")
QSYS_DEBUG_LEVEL = 0 ! show me less
do Ef = -3.0 , 3.025 , 0.015
    ! Update hamiltonian elemenents value
    call qt%qsystem%update_lattice(c_matrix=coupling)
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)

    print*,"Energy:",Ef,sum(qt%Tn(:))
    write(111,"(100f20.6)"),Ef,sum(qt%Tn(:))
enddo
close(111)


! ----------------------------------------------------------
! X. Clean memory...
! ----------------------------------------------------------

call qt%destroy_system()
print*,"Generating plots..."
print*,"Plotting band structure..."
call system("cd "//output_folder//"; ./plot_bands.py")
print*,"Plotting eigenvectors..."
call system("cd "//output_folder//"; ./plot_eigenvecs.py")
print*,"Plotting Transmission..."
call system("cd "//output_folder//"; ./plot_T.py")
print*,"Use Viewer program to see the structure and created leads."
contains

! ---------------------------------------------------------------------------
! Implement coupling
! Taken from Kwant tutorial: http://kwant-project.org/doc/1.0/tutorial/tutorial5
! ---------------------------------------------------------------------------
logical function coupling(atomA,atomB,coupling_mat)
    use modcommons
    implicit none
    type(qatom) :: atomA,atomB
    complex*16  :: coupling_mat(:,:) ! you must overwrite this variable
    ! local variables
    integer         :: xdiff,ydiff
    doubleprecision :: dydiff,dxdiff,t0,mu,Delta

    ! Calculate distance between atoms in units of dx.
    dxdiff = (atomA%atom_pos(1)-atomB%atom_pos(1))/dx
    dydiff = (atomA%atom_pos(2)-atomB%atom_pos(2))/dx
    ! Convert it to integers
    xdiff = NINT(dxdiff)
    ydiff = NINT(dydiff)
    ! default return value
    coupling       = .false.
    coupling_mat   = 0.0

    t0    = 1.0 ! hoping energy
    mu    = 0.7 ! chemical potential
    Delta = 0.1 ! superconducting order parameter

    if( xdiff == 0 .and. ydiff == 0 ) then
        coupling     = .true.
        coupling_mat = (4*t0 - mu)* MAT_SZ + Delta*MAT_SX
    else if( abs(xdiff) ==  1 .and. ydiff == 0 ) then
        coupling     = .true.
        coupling_mat = -t0 *(MAT_SZ)
    else if( xdiff ==  0 .and. abs(ydiff) == 1 ) then
        coupling     = .true.
        coupling_mat = -t0 *(MAT_SZ)
    endif
end function coupling
end program supercond
