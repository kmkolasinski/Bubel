! ------------------------------------------------------ !
! Quantulaba - strange_graphene.f90 - Krzysztof Kolasinski 2015
!
! Graphene unit cell consists of two atoms A and B, those
! atoms can be viewed as a two spin states of some quantum
! object. In that case Graphene lattice can be presented
! as a simple square lattice with propper coupling between
! nearest neightbours.
! ------------------------------------------------------ !

program transporter

use modunits ! unit conversion tools
use modsys   ! eigen values
use modlead  ! bandgap structure
use modshape
implicit none
character(*),parameter :: output_folder = "strange_graphene_output/"
type(qsys)                 :: qsystem
type(qshape)               :: lead_shape
type(qlead)                :: lead
doubleprecision            :: zeros_vector(200)
doubleprecision            :: a_dx,a_Emin,a_Emax
integer                    :: no_expected_states
integer ,parameter         :: nx = 50
integer ,parameter         :: ny = 50
doubleprecision,parameter  :: dx = 5.0 ! [nm]
doubleprecision,parameter  :: alpha30 =  30.0/180.0*M_PI

integer , dimension(nx,ny) :: gindex ! converts local index (i,j) to global index
integer :: i,j,k,m
doubleprecision            :: lead_translation_vec(3),x,y,z
character(300) :: line


! Initalize system
call qsystem%init()
! ----------------------------------------------------------
! 1. Create mesh - loop over width and height of the lattice
! ----------------------------------------------------------
k      = 0
gindex = 0
do i = 1 , nx
do j = 1 , ny
    ! Initalize atom structure with position of the atom. Now instead of atoms A and B we have
    ! one atom with two spins.
    call qsystem%qatom%init((/ (i-1) * dx , (j-1) * dx + (i-1) * dx * 0.0 , 0.0 * dx /),no_in_states=2)
    ! Add atom to the system.
    call qsystem%add_atom(qsystem%qatom)
    k           = k + 1
    gindex(i,j) = k
enddo
enddo
! ----------------------------------------------------------
! 2. Construct logical connections between sites on the mesh.
! ----------------------------------------------------------
qsystem%qnnbparam%box = (/2*dx,2*dx,0.0D0/)
call qsystem%make_lattice(qsystem%qnnbparam,c_matrix=connect_matrix)
! Save lattice to file to see if constructed system is OK!
call qsystem%save_lattice(output_folder//"lattice.xml")

! ----------------------------------------------------------
! 3. Find eigenvalues of the system.
! You must provide range of energy (Emin,Emax) and
! number of expected states in that region. Eigen problem
! is solved with FEAST library.
! ----------------------------------------------------------
no_expected_states = 40
! we choose small window
a_Emin =-3.0+0.1
a_Emax =-2.97+0.1
call qsystem%calc_eigenproblem(a_Emin,a_Emax,no_expected_states)

! Write founded eigenstates to file if there are states in the range (Emin,Emax)
if(qsystem%no_eigenvalues > 0) then
zeros_vector = 0
open(unit=222,file=output_folder//"eigenvecs.dat")
do i = 1 , nx
do j = 1 , ny
    ! Remapping to graphen lattice
    x = (qsystem%atoms(gindex(i,j))%atom_pos(1) + qsystem%atoms(gindex(i,j))%atom_pos(2)*sin(alpha30))/dx
    y = (qsystem%atoms(gindex(i,j))%atom_pos(2) - sin(alpha30)*qsystem%atoms(gindex(i,j))%atom_pos(1))/dx
    z = 0
    if(qsystem%atoms(gindex(i,j))%bActive)then
        ! get unique ID of current site
        k = qsystem%atoms(gindex(i,j))%globalIDs(2)  ! spin = 2
        m = qsystem%atoms(gindex(i,j))%globalIDs(1)  ! spin = 1
        write(222,"(500e20.6)"),x,y,z,abs(qsystem%eigenvecs(k,:))**2+abs(qsystem%eigenvecs(m,:))**2
    else
        ! fill with zeros empty spaces
        write(222,"(500e20.6)"),x,y,z,zeros_vector(1:qsystem%no_eigenvalues)
    endif
enddo
write(222,*),"" ! GNUPLOT necessary line break
enddo
close(222)
endif ! end of if are there eigenstates?



! ----------------------------------------------------------
! 4. Use generated mesh to calculate the band structure
! in the region of homogenous lead.
! ----------------------------------------------------------
! Setup shape and initialize lead
lead_translation_vec = (/  dx , 0.0D0 , 0.0D0 /)
call lead_shape%init_range_3d((/ -0.5*dx , 0.0D0 , 0.0D0 /),lead_translation_vec)
call lead%init_lead(lead_shape,lead_translation_vec,qsystem%atoms)
call lead%save_lead(output_folder//"lead.xml")
a_Emin =-3.0
a_Emax = 3.0
call lead%bands(output_folder//"bands.dat",-M_PI,+M_PI,M_PI/160.0,a_Emin,a_Emax)

! ----------------------------------------------------------
! X. Clean memory...
! ----------------------------------------------------------
call lead%destroy()
call qsystem%destroy()

print*,"Generating plots..."
print*,"Plotting band structure..."
call system("cd "//output_folder//"; ./plot_bands.py")
print*,"Plotting eigenvectors..."
call system("cd "//output_folder//"; ./plot_eigenvecs.py")
print*,"Use Viewer program to see the structure and crated lead."
contains

! ----------------------------------------------------------
! Define coupling between atoms
! ----------------------------------------------------------
logical function connect_matrix(atomA,atomB,coupling_mat)
    use modatom
    implicit none
    type(qatom) :: atomA,atomB
    complex*16  :: coupling_mat(:,:) ! you must overwrite this variable
    ! local variables
    integer         :: xdiff,ydiff
    doubleprecision :: dydiff,dxdiff,t0
    ! Calculate distance between atoms in units of dx.
    dxdiff = (atomA%atom_pos(1)-atomB%atom_pos(1))/dx
    dydiff = (atomA%atom_pos(2)-atomB%atom_pos(2))/dx
    ! Convert it to integers
    xdiff = NINT(dxdiff)
    ydiff = NINT(dydiff)

    ! default return value
    connect_matrix = .false.
    coupling_mat   = 0.0
    t0 = -1
    if( xdiff == 0 .and. ydiff == 0 ) then
        connect_matrix      = .true.
        coupling_mat = t0 * MAT_SX
    else if( abs(xdiff) ==  1 .and. ydiff == 0 ) then
        connect_matrix = .true.
        coupling_mat = 0.5*t0 *(MAT_SX + xdiff*II*MAT_SY)
    else if( xdiff ==  0 .and. abs(ydiff) == 1 ) then
        connect_matrix = .true.
        coupling_mat = 0.5*t0 *(MAT_SX + ydiff*II*MAT_SY)
    endif
end function connect_matrix

end program transporter
