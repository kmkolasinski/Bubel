! ------------------------------------------------------ !
! Bubel - rashba_2deg.f90 - Krzysztof Kolasinski 2015
!
! We want to construct hamiltonian for effective mass
! Schroedinger equation with Rashba interaction added.
! and magnetic field
! ------------------------------------------------------ !

program rashba

use modunits ! unit conversion tools
use modshape
use modscatter  ! transport and eigen problem
implicit none
character(*),parameter :: output_folder = "rashba_2deg_output/"
type(qscatter)             :: qt
type(qshape)               :: lead_shape

doubleprecision            :: a_dx,a_Emin,a_Emax
integer                    :: no_expected_states,flagA,flagB
integer ,parameter         :: nx = 30
integer ,parameter         :: ny = 20
doubleprecision,parameter  :: dx = 1.0 ! [nm]
integer                    :: i,j,k
doubleprecision            :: lead_translation_vec(3),Ef
character(300) :: line

 ! Use atomic units in effective band model -> see modunit.f90 for more details
 call modunits_set_GaAs_params()
 a_dx = dx * L2LA ! convert dx from nm to atomic units
 ! Initalize system
 call qt%init_system()
 ! ----------------------------------------------------------
 ! 1. Create mesh - loop over width and height of the lattice
 ! ----------------------------------------------------------

 do i = 1 , nx
 do j = 1 , ny
    ! Initalize atom structure with position of the atom. Note that here atom has two spin states
    call qt%qatom%init((/ (i-1) * dx , (j-1) * dx + (i-1) * dx * 0.0 , 0.0 * dx /),no_in_states=2)
    ! Add atom to the system.
    call qt%qsystem%add_atom(qt%qatom)
 enddo
 enddo

 ! ----------------------------------------------------------
 ! 2. Construct logical connections between sites on the mesh.
! ----------------------------------------------------------
! Set criterium for the nearest neightbours "radius" search algorithm.
! Same as above qt%qnnbparam is a auxiliary variable of type(nnb_params) - more details in modsys.f90
! This structure is responsible for different criteria of nearest neighbour searching
qt%qnnbparam%box = (/2*dx,2*dx,0.0D0/) ! do not search for the sites far than (-dx:+dx) direction

! Setup connections between sites with provided by you function "connect", see below for example.
! here we use matrix representation for more than one orbital
call qt%qsystem%make_lattice(qt%qnnbparam,c_matrix=coupling)

! ----------------------------------------------------------
! 3. Use generated mesh to calculate the band structure
! in the region of homogenous lead.
! ----------------------------------------------------------
! Setup shape and initialize lead
call lead_shape%init_rect(SHAPE_RECTANGLE_XY,-0.5*dx,0.5*dx,-0.5*dx,(ny+1)*dx)
! Lead needs to know where it is (lead_shape) and using this information it will
! create propper matrices and connections using list of atoms
lead_translation_vec = (/  dx , 0.0D0 , 0.0D0 /)
call qt%add_lead(lead_shape,lead_translation_vec)

call qt%leads(1)%save_lead(output_folder//"lead.xml")
a_Emin =-0.01 / A0 / 1000.0 ! converting from [meV] to atomic units
a_Emax =10.00 / A0 / 1000.0
call qt%leads(1)%bands(output_folder//"bands.dat",-M_PI,+M_PI,M_PI/100.0,a_Emin,a_Emax)


call lead_shape%init_rect(SHAPE_RECTANGLE_XY,(nx-0.5-1)*dx,(0.5+nx-1)*dx,-0.5*dx,(ny+1)*dx)
! Lead needs to know where it is (lead_shape) and using this information it will
! create propper matrices and connections using list of atoms
lead_translation_vec = (/  -dx , 0.0D0 , 0.0D0 /)
call qt%add_lead(lead_shape,lead_translation_vec)

! Save lattice to file to see if constructed system is OK!
call qt%save_system(output_folder//"system.xml")

!QSYS_DEBUG_LEVEL = 1 ! See more
open(unit=111,file=output_folder//"T.dat")
print*,"Performing energy scan..."
do Ef = 0.0 , 0.006 , 0.00005
    call qt%qsystem%update_lattice(c_matrix=coupling)
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)
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
print*,"Plotting Transmission..."
call system("cd "//output_folder//"; ./plot_T.py")
print*,"Use Viewer program to see the structure and created leads."
contains

! ---------------------------------------------------------------------------
! This coupling function
! ---------------------------------------------------------------------------
logical function coupling(atomA,atomB,coupling_mat)
    use modcommons
    implicit none
    type(qatom) :: atomA,atomB

    complex*16  :: coupling_mat(:,:) ! you must overwrite this variable
    ! local variables
    integer         :: xdiff,ydiff
    doubleprecision :: dydiff,dxdiff,t0,y,rs,bz

    ! Calculate distance between atoms in units of dx.
    dxdiff = (atomA%atom_pos(1)-atomB%atom_pos(1))/dx
    dydiff = (atomA%atom_pos(2)-atomB%atom_pos(2))/dx
    ! Convert it to integers
    xdiff = NINT(dxdiff)
    ydiff = NINT(dydiff)

    ! default return value
    coupling = .false.
    coupling_mat   = 0.0

    t0 = 1/(2*m_eff*a_dx**2)
    rs = 0.06/(2*a_dx) ! rashba coupling
    bz = 0.0001
    if( xdiff == 0 .and. ydiff == 0 ) then
        coupling      = .true.
        coupling_mat = 4*t0 * MAT_DIAG + bz * MAT_SZ
    else if( abs(xdiff) ==  1 .and. ydiff == 0 ) then
        coupling = .true.
        coupling_mat = -t0 * MAT_DIAG + xdiff*rs*II*MAT_SY
    else if( xdiff ==  0 .and. abs(ydiff) == 1 ) then
        coupling = .true.
        coupling_mat = -t0 * MAT_DIAG - ydiff*rs*II*MAT_SX
    endif

end function coupling

end program rashba
