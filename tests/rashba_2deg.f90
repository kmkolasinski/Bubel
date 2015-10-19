! ------------------------------------------------------ !
! Quantulaba - rashba_2deg.f90 - Krzysztof Kolasinski 2015
!
! We want to construct hamiltonian for effective mass
! Schroedinger equation with Rashba interaction added.
! We neglect the efect of the magnetic field.
! ------------------------------------------------------ !

program transporter

 use modunits ! unit conversion tools
 use modsys   ! eigen values
 use modlead  ! bandgap structure
 use modshape
 implicit none
 character(*),parameter :: output_folder = "rashba_2deg_output/"
 type(qsys)                 :: qsystem
 type(qshape)               :: lead_shape
 type(qlead)                :: lead
 doubleprecision            :: a_dx,a_Emin,a_Emax
 integer                    :: no_expected_states,flagA,flagB
 integer ,parameter         :: nx = 30
 integer ,parameter         :: ny = 20
 doubleprecision,parameter  :: dx = 1.0 ! [nm]
 integer                    :: i,j,k
 doubleprecision            :: lead_translation_vec(3)
 character(300) :: line

 ! Use atomic units in effective band model -> see modunit.f90 for more details
 call modunits_set_GaAs_params()
 a_dx = dx * L2LA ! convert dx from nm to atomic units
 ! Initalize system
 call qsystem%init()
 ! ----------------------------------------------------------
 ! 1. Create mesh - loop over width and height of the lattice
 ! ----------------------------------------------------------

 do i = 1 , nx
 do j = 1 , ny
    ! Initalize atom structure with position of the atom. Note that here atom has two spin states
    call qsystem%qatom%init((/ (i-1) * dx , (j-1) * dx + (i-1) * dx * 0.0 , 0.0 * dx /),no_in_states=2)
    ! Add atom to the system.
    call qsystem%add_atom(qsystem%qatom)
 enddo
 enddo

 ! ----------------------------------------------------------
 ! 2. Construct logical connections between sites on the mesh.
 ! ----------------------------------------------------------
 ! Set criterium for the nearest neightbours "radius" search algorithm.
 ! Same as above qsystem%qnnbparam is a auxiliary variable of type(nnb_params) - more details in modsys.f90
 ! This structure is responsible for different criteria of nearest neighbour searching
  qsystem%qnnbparam%box = (/2*dx,2*dx,0.0D0/) ! do not search for the sites far than (-dx:+dx) direction

 ! Setup connections between sites with provided by you function "connect", see below for example.
 !call qsystem%make_lattice(qsystem%qnnbparam,c_default=connect)
 call qsystem%make_lattice(qsystem%qnnbparam,c_matrix=connect_matrix)
 ! Save lattice to file to see if constructed system is OK!
 ! Use plot_lattice.py to see the results.
 call qsystem%save_lattice(output_folder//"lattice.xml")

 ! ----------------------------------------------------------
 ! 3. Use generated mesh to calculate the band structure
 ! in the region of homogenous lead.
 ! ----------------------------------------------------------
 ! Setup shape and initialize lead
 call lead_shape%init_rect(SHAPE_RECTANGLE_XY,-0.5*dx,0.5*dx,-0.5*dx,(ny+1)*dx)
 ! Lead needs to know where it is (lead_shape) and using this information it will
 ! create propper matrices and connections using list of atoms
 lead_translation_vec = (/  dx , 0.0D0 , 0.0D0 /)
 call lead%init_lead(lead_shape,lead_translation_vec,qsystem%atoms)

 call lead%save_lead(output_folder//"lead.xml")
 a_Emin =-0.01 / A0 / 1000.0 ! converting from [meV] to atomic units
 a_Emax =10.00 / A0 / 1000.0
 call lead%bands(output_folder//"bands.dat",-M_PI,+M_PI,M_PI/100.0,a_Emin,a_Emax)
 ! ----------------------------------------------------------
 ! X. Clean memory...
 ! ----------------------------------------------------------
 call lead%destroy()
 call qsystem%destroy()

 print*,"Generating plots..."
 print*,"Plotting band structure..."
 call system("cd "//output_folder//"; ./plot_bands.py")
 print*,"Use Viewer program to see the structure and crated lead."
 contains

! ---------------------------------------------------------------------------
! This function decides if site A (called here atomA) with spin s1 has hoping
! to atom B with spin s2, and what is the value of the coupling.
! If there is no interaction between them returns false, otherwise true.
! ---------------------------------------------------------------------------
logical function connect_matrix(atomA,atomB,coupling_mat)
    use modcommons
    implicit none
    type(qatom) :: atomA,atomB

    complex*16  :: coupling_mat(:,:) ! you must overwrite this variable
    ! local variables
    integer         :: xdiff,ydiff
    doubleprecision :: dydiff,dxdiff,t0,y,rs

    ! Calculate distance between atoms in units of dx.
    dxdiff = (atomA%atom_pos(1)-atomB%atom_pos(1))/dx
    dydiff = (atomA%atom_pos(2)-atomB%atom_pos(2))/dx
    ! Convert it to integers
    xdiff = NINT(dxdiff)
    ydiff = NINT(dydiff)

    ! default return value
    connect_matrix = .false.
    coupling_mat   = 0.0

    t0 = 1/(2*m_eff*a_dx**2)
    rs = 0.05/(2*a_dx) ! rashba coupling

    if( xdiff == 0 .and. ydiff == 0 ) then
        connect_matrix      = .true.
        coupling_mat = 4*t0 * MAT_DIAG
    else if( abs(xdiff) ==  1 .and. ydiff == 0 ) then
        connect_matrix = .true.
        coupling_mat = -t0 * MAT_DIAG + xdiff*rs*II*MAT_SY
    else if( xdiff ==  0 .and. abs(ydiff) == 1 ) then
        connect_matrix = .true.
        coupling_mat = -t0 * MAT_DIAG - ydiff*rs*II*MAT_SX
    endif

end function connect_matrix

end program transporter
