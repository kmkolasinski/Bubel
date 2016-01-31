! ------------------------------------------------------ !
! Quantulaba - superconductor_transport.f90 - Krzysztof Kolasinski 2015
!
! Example based on Kwant tutorial i.e. how to deal
! with multiple lattices in Quantulaba:
! http://kwant-project.org/doc/1/tutorial/tutorial5
! Here we create two lattices displaced in z direction
! for which the one with z=0 is responsible for -1/2 spin value
! and z=dx +1/2.
! ------------------------------------------------------ !
program supercond
use modunits     ! unit conversion tools
use modscatter   ! eigen values and transport
use modshape
use modsys
implicit none
character(*),parameter :: output_folder = "supercnd_trans_output/"
type(qscatter)             :: qt
type(qshape)               :: lead_shape
type(qsys)                 :: lead_e,lead_h ! storage for atoms in the electron and hole leads

doubleprecision            :: a_Emin,a_Emax

integer ,parameter         :: nx = 10
integer ,parameter         :: ny = 10
doubleprecision,parameter  :: dx = 1.0 ! [nm]

integer                    :: i,j,k,m
doubleprecision            :: lead_translation_vec(3),Ef,T
character(300) :: line
 ! Initalize system
call qt%init_system()
call lead_e%init()
call lead_h%init()

! ----------------------------------------------------------
! 1. Create mesh - loop over width and height of the lattice
! ----------------------------------------------------------

do i = 1 , nx
do j = 1 , ny
    ! Initalize atom structure with position of the atom.
    call qt%qatom%init((/ (i-1) * dx , (j-1) * dx + (i-1) * dx * 0.0 , 0.0D0 /),flag=1)
    call qt%qsystem%add_atom(qt%qatom)
    ! store electrons lead "atoms"
    if( i == 1 ) call lead_e%add_atom(qt%qatom)
    ! Add second lattice (above, which will responsible for +1/2 spin)
    call qt%qatom%init((/ (i-1) * dx , (j-1) * dx + (i-1) * dx * 0.0 , 1.0D0 /),flag=2)
    call qt%qsystem%add_atom(qt%qatom)
    ! store holes lead "atoms"
    if( i == 1 ) call lead_h%add_atom(qt%qatom)
enddo
enddo
! ----------------------------------------------------------
! 2. Construct logical connections between sites on the mesh.
! ----------------------------------------------------------
qt%qnnbparam%box = (/2*dx,2*dx,2*dx/) ! do not search for the sites far than (-dx:+dx) direction
call qt%qsystem%make_lattice(qt%qnnbparam,c_matrix=coupling)


! ----------------------------------------------------------
! 3. Add leads to the system
! ----------------------------------------------------------
lead_translation_vec = (/  dx , 0.0D0 , 0.0D0 /)
call lead_shape%init_atoms_list(lead_e%atoms)
call qt%add_lead(lead_shape,lead_translation_vec)

lead_translation_vec = (/  dx , 0.0D0 , 0.0D0 /)
call lead_shape%init_atoms_list(lead_h%atoms)
call qt%add_lead(lead_shape,lead_translation_vec)

lead_translation_vec = (/ -dx , 0.0D0 , 0.0D0 /)
call lead_shape%init_range_3d((/ (nx+0.5-1)*dx , 0.0D0 , 0.0D0 /),lead_translation_vec)
call qt%add_lead(lead_shape,lead_translation_vec)


! ----------------------------------------------------------
! 4. Use generated mesh to calculate the band structure
! in the region of homogenous leads: electrons,holes,superconductor
! ----------------------------------------------------------
a_Emin =-10.0
a_Emax = 10.0
call qt%leads(1)%bands(output_folder//"bands_e.dat",-M_PI,+M_PI,M_PI/160.0,a_Emin,a_Emax)
call qt%leads(2)%bands(output_folder//"bands_h.dat",-M_PI,+M_PI,M_PI/160.0,a_Emin,a_Emax)
call qt%leads(3)%bands(output_folder//"bands_s.dat",-M_PI,+M_PI,M_PI/160.0,a_Emin,a_Emax)


! Save lattice to file to see if constructed system is OK!
call qt%save_system(output_folder//"system.xml")

! ----------------------------------------------------------
! 5. Solve scattering problem for a given incident energy
! and save obtained electron density to file.
! ----------------------------------------------------------
Ef = 0.1
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
do Ef = 0.0 , 0.2 , 0.002
    ! Update hamiltonian elemenents value
    call qt%qsystem%update_lattice(c_matrix=coupling)
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)

    ! Current =  N - Ree - Reh
    T =  qt%leads(1)%no_in_modes - qt%leads(1)%totalT + qt%leads(2)%totalT
    print*,"Energy:",Ef,T
    write(111,"(100f20.6)"),Ef,T
enddo
close(111)


! ----------------------------------------------------------
! X. Clean memory...
! ----------------------------------------------------------

call qt%destroy_system()
call lead_e%destroy()
call lead_h%destroy()
print*,"Generating plots..."
print*,"Plotting band structure..."
call system("cd "//output_folder//"; python plot_bands.py")
print*,"Plotting Transmission..."
call system("cd "//output_folder//"; python plot_T.py")
print*,"Use Viewer program to see the structure and created leads."
contains

! ---------------------------------------------------------------------------
! Implement coupling between atoms
! This was taken from Kwant tutorial:
! http://kwant-project.org/doc/1/tutorial/tutorial5#lattice-description-using-different-lattices
! ---------------------------------------------------------------------------
logical function coupling(atomA,atomB,coupling_mat)
    use modcommons
    implicit none
    type(qatom) :: atomA,atomB
    complex*16  :: coupling_mat(:,:) ! you must overwrite this variable
    ! local variables
    integer         :: xdiff,ydiff
    integer         :: flagA,flagB

    doubleprecision :: dydiff,dxdiff,t0,mu,Delta,barrier

    ! Calculate distance between atoms in units of dx.
    dxdiff = (atomA%atom_pos(1)-atomB%atom_pos(1))/dx
    dydiff = (atomA%atom_pos(2)-atomB%atom_pos(2))/dx
    ! Convert it to integers
    xdiff = NINT(dxdiff)
    ydiff = NINT(dydiff)
    ! default return value
    coupling       = .false.
    coupling_mat   = 0.0
    flagA =  atomA%flag
    flagB =  atomB%flag

    t0    = 1.0 ! hoping energy
    mu    = 0.4 ! chemical potential
    Delta = 0.1 ! superconducting order parameter
    ! Potential barrier definition
    barrier = 0
    if( atomA%atom_pos(1) >= nx/2.0-1 .and. atomA%atom_pos(1) < nx/2.0 ) barrier = 1.5

    if( xdiff == 0 .and. ydiff == 0 ) then
        coupling     = .true.
        coupling_mat = (4*t0 - mu + barrier ) * MAT_SZ(flagA,flagB) + Delta*MAT_SX(flagA,flagB)
        ! do not create hopping for normal conductor part
        if(flagA /= flagB .and. atomA%atom_pos(1) < nx/2 ) coupling     = .false.
    else if( abs(xdiff) ==  1 .and. ydiff == 0 .and. flagA == flagB ) then
        coupling     = .true.
        coupling_mat = -t0 *(MAT_SZ(flagA,flagB))
    else if( xdiff ==  0 .and. abs(ydiff) == 1 .and. flagA == flagB ) then
        coupling     = .true.
        coupling_mat = -t0 *(MAT_SZ(flagA,flagB))
    endif
end function coupling
end program supercond
