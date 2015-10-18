!! ------------------------------------------------------ !
!! Quantulaba - rashba_2deg.f90 - Krzysztof Kolasinski 2015
!!
!! We want to construct hamiltonian for effective mass
!! Schroedinger equation with Rashba interaction added.
!! We neglect the efect of the magnetic field.
!! ------------------------------------------------------ !
!
!program transporter
!
! use modunits ! unit conversion tools
! use modsys   ! eigen values
! use modlead  ! bandgap structure
! use modshape
! use modscatter
! implicit none
! character(*),parameter :: output_folder = ""
!
! type(qshape)               :: lead_shape
! type(qlead)                :: lead
! type(qscatter)             :: qt
! doubleprecision            :: a_dx,a_Emin,a_Emax,Ef
! integer                    :: no_expected_states,flagA,flagB
! integer ,parameter         :: nx = 200
! integer ,parameter         :: ny = 50
! doubleprecision,parameter  :: dx = 0.5 ! [nm]
! integer                    :: i,j,k,m
! doubleprecision            :: lead_translation_vec(3)
! character(300) :: line
!
! ! Use atomic units in effective band model -> see modunit.f90 for more details
! call modunits_set_GaAs_params()
! a_dx = dx * L2LA ! convert dx from nm to atomic units
! ! Initalize system
! call qt%init_system()
! ! ----------------------------------------------------------
! ! 1. Create mesh - loop over width and height of the lattice
! ! ----------------------------------------------------------
!
! do i = 1 , nx
! do j = 1 , ny
!    ! Initalize atom structure with position of the atom. Note that here atom has two spin states
!    call qt%qatom%init((/ (i-1) * dx , (j-1) * dx + (i-1) * dx * 0.0 , 0.0 * dx /),no_in_states=2)
!    ! Add atom to the system.
!    if( i > 75 .and. i < 125 .and. j < 40 ) qt%qatom%bActive = .false.
!    call qt%qsystem%add_atom(qt%qatom)
! enddo
! enddo
!
! ! ----------------------------------------------------------
! ! 2. Construct logical connections between sites on the mesh.
! ! ----------------------------------------------------------
! ! Set criterium for the nearest neightbours "radius" search algorithm.
! ! Same as above qsystem%qnnbparam is a auxiliary variable of type(nnb_params) - more details in modsys.f90
! ! This structure is responsible for different criteria of nearest neighbour searching
!  qt%qnnbparam%box = (/2*dx,2*dx,0.0D0/) ! do not search for the sites far than (-dx:+dx) direction
!
! ! Setup connections between sites with provided by you function "connect", see below for example.
! !call qsystem%make_lattice(qsystem%qnnbparam,c_default=connect)
! call qt%qsystem%make_lattice(qt%qnnbparam,c_matrix=connect_matrix)
! ! Save lattice to file to see if constructed system is OK!
! ! Use plot_lattice.py to see the results.
!! call qt%qsystem%save_lattice(output_folder//"lattice.xml")
!
! ! ----------------------------------------------------------
! ! 3. Use generated mesh to calculate the band structure
! ! in the region of homogenous lead.
! ! ----------------------------------------------------------
! ! Setup shape and initialize lead
! call lead_shape%init_rect(SHAPE_RECTANGLE_XY,-0.5*dx,0.5*dx,-0.5*dx,(ny+1)*dx)
! ! Lead needs to know where it is (lead_shape) and using this information it will
! ! create propper matrices and connections using list of atoms
! lead_translation_vec = (/  dx , 0.0D0 , 0.0D0 /)
! call qt%add_lead(lead_shape,lead_translation_vec)
!
!
! call lead_shape%init_rect(SHAPE_RECTANGLE_XY,(-0.5+nx-1)*dx,(0.5+nx-1)*dx,-0.5*dx,(ny+1)*dx)
! ! Lead needs to know where it is (lead_shape) and using this information it will
! ! create propper matrices and connections using list of atoms
! lead_translation_vec = (/ -dx , 0.0D0 , 0.0D0 /)
! call qt%add_lead(lead_shape,lead_translation_vec)
!
! call qt%save_system("system.xml")
!
!! call qt%leads(1)%save_lead(output_folder//"lead.xml")
! a_Emin =-0.01 / A0 / 1000.0 ! converting from [meV] to atomic units
! a_Emax =10.00 / A0 / 1000.0
!! call qt%leads(1)%bands(output_folder//"bands.dat",-M_PI,+M_PI,M_PI/100.0,a_Emin,a_Emax)
! ! ----------------------------------------------------------
! ! X. Clean memory...
! ! ----------------------------------------------------------
! Ef = 0.01
! call qt%calculate_modes(Ef)
! call qt%solve(1,Ef)
!
! k = qt%leads(1)%no_in_modes
!do i = 1 , qt%qsystem%no_atoms
!    if(qt%qsystem%atoms(i)%bActive) then
!        m = qt%qsystem%atoms(i)%globalIDs(1)
!        write(111,"(100e20.6)"),qt%qsystem%atoms(i)%atom_pos,sum(qt%densities(:,m)),qt%densities(:,m)
!    endif
!enddo
!
!do i = 1 , size(qt%qauxvec)
!    qt%qauxvec(i) = sum(qt%densities(:,i))
!enddo
!call qt%qsystem%save_data("densities.xml",array2d=qt%densities,array1d=qt%qauxvec)
!call qt%destroy_system()
!
! print*,"Generating plots..."
! print*,"Plotting band structure..."
!! call system("cd "//output_folder//"; ./plot_bands.py")
! print*,"Use Viewer program to see the structure and crated lead."
! contains
!
!! ---------------------------------------------------------------------------
!! This function decides if site A (called here atomA) with spin s1 has hoping
!! to atom B with spin s2, and what is the value of the coupling.
!! If there is no interaction between them returns false, otherwise true.
!! ---------------------------------------------------------------------------
!logical function connect_matrix(atomA,atomB,coupling_mat)
!    use modatom
!    implicit none
!    type(qatom) :: atomA,atomB
!
!    complex*16  :: coupling_mat(:,:) ! you must overwrite this variable
!    ! local variables
!    integer         :: xdiff,ydiff
!    doubleprecision :: dydiff,dxdiff,t0,y,rs
!
!    ! Calculate distance between atoms in units of dx.
!    dxdiff = (atomA%atom_pos(1)-atomB%atom_pos(1))/dx
!    dydiff = (atomA%atom_pos(2)-atomB%atom_pos(2))/dx
!    ! Convert it to integers
!    xdiff = NINT(dxdiff)
!    ydiff = NINT(dydiff)
!
!    ! default return value
!    connect_matrix = .false.
!    coupling_mat   = 0.0
!
!    t0 = 1/(2*m_eff*a_dx**2)
!    rs = 0.05/(2*a_dx) ! rashba coupling
!
!    if( xdiff == 0 .and. ydiff == 0 ) then
!        connect_matrix      = .true.
!        coupling_mat = 4*t0 * MAT_DIAG
!    else if( abs(xdiff) ==  1 .and. ydiff == 0 ) then
!        connect_matrix = .true.
!        coupling_mat = -t0 * MAT_DIAG + xdiff*rs*II*MAT_SY
!    else if( xdiff ==  0 .and. abs(ydiff) == 1 ) then
!        connect_matrix = .true.
!        coupling_mat = -t0 * MAT_DIAG - ydiff*rs*II*MAT_SX
!    endif
!
!end function connect_matrix
!
!end program transporter



! ------------------------------------------------------ !
! Quantulaba - carbon_nanotube.f90 - Krzysztof Kolasinski 2015
!
! Example of hamiltonian creation for simple carbon
! nanotube.
! ------------------------------------------------------ !

program transporter

 use modunits ! unit conversion tools
 use modsys   ! eigen values
 use modlead  ! bandgap structure
 use modshape
 use modscatter
 implicit none
 character(*),parameter :: output_folder = ""
 type(qsys)                 :: tmpsystem
 type(qshape)               :: lead_shape,ashape

 type(qscatter)             :: qt
 doubleprecision            :: posA(3),posB(3),deltaR(3),Emin,Emax,Ef
 integer                    :: flagA,flagB
 logical                    :: activeA
 integer :: i,j,k , M
 doubleprecision            :: lead_translation_vec(3)
 character(300) :: line

 ! Initalize system
 call qt%init_system()
 call tmpsystem%init()

!  ----------------------------------------------------------
!  1. Create mesh - read unit cell from file
!  ----------------------------------------------------------
open(unit=3,file=output_folder//"nt-17-0-1.xyz")
read(3,*) line
read(3,*) line
! ----------------------------------
! Reading unit cell
! ----------------------------------
do i=1,17*2
    read(3,*) line, (posA(j), j=1, 3)
    read(3,*) line, (posB(j), j=1, 3)

    call tmpsystem%qatom%init(posA,flag=1)
    if(i==1) tmpsystem%qatom%flag_aux0 = 1
    call tmpsystem%add_atom(tmpsystem%qatom)
    call tmpsystem%qatom%init(posB,flag=2)
!    if(i==2) tmpsystem%qatom%bActive = .false.
    call tmpsystem%add_atom(tmpsystem%qatom)
enddo
close(3)

call ashape%init_atoms_list(tmpsystem%atoms(1:tmpsystem%no_atoms))


! ----------------------------------
! Translating unit cell X-times
! ----------------------------------
M = 20
do k = 1 , M
    deltaR = (/0.0,0.0,(k-1)*4.26/)
    do i=1,tmpsystem%no_atoms
        posA  = tmpsystem%atoms(i)%atom_pos
        flagA = tmpsystem%atoms(i)%flag
        activeA = tmpsystem%atoms(i)%bActive
        call qt%qatom%init(posA+deltaR,flag=flagA,bActive=activeA)
        qt%qatom%flag_aux0  = tmpsystem%atoms(i)%flag_aux0
!        if(k>10 .and. k < 12 .and. i<3) qt%qatom%bActive = .false.
        call qt%qsystem%add_atom(qt%qatom)
    enddo
enddo


!  ----------------------------------------------------------
!  2. Set nearest neightbour search radius and make Hmiltonian
!  ----------------------------------------------------------
qt%qnnbparam%distance   = 2.0
qt%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_DISTANCE
call qt%qsystem%make_lattice(qt%qnnbparam,c_simple=connect_simple)
!call qt%qsystem%save_lattice(output_folder//"lattice.xml")

!  ----------------------------------------------------------
!  3. Set up unit cell
!  ----------------------------------------------------------

lead_translation_vec = (/  0.0D0 , 0.0D0 , 4.26D0 /)
posA = (/0.0,0.0,-0.5/)
posB = lead_translation_vec ! direction and range
call lead_shape%init_range_3d(posA,posB)
call qt%add_lead(lead_shape,lead_translation_vec)
!call qt%leads(1)%save_lead(output_folder//"lead.xml")


lead_translation_vec = (/  0.0D0 , 0.0D0 , -4.26D0 /)
posA = (/0.0,0.0,-0.5/) - lead_translation_vec*M
posB = lead_translation_vec ! direction and range
call lead_shape%init_range_3d(posA,posB)
call qt%add_lead(lead_shape,lead_translation_vec)

call qt%save_system("system.xml")
Ef = 0.14;
call qt%calculate_modes(Ef)
call qt%solve(2,Ef)


k = qt%leads(1)%no_in_modes
do i = 1 , qt%qsystem%no_atoms
    if(qt%qsystem%atoms(i)%bActive) then
        m = qt%qsystem%atoms(i)%globalIDs(1)
        write(111,"(100e20.6)"),qt%qsystem%atoms(i)%atom_pos,sum(qt%densities(:,m)),qt%densities(:,m)
    endif
enddo

!  ----------------------------------------------------------
!  4. Calculate band structure
!  ----------------------------------------------------------
Emin = -3.0
Emax =  3.0
call qt%leads(1)%bands(output_folder//"bands.dat",-M_PI,+M_PI,M_PI/60.0,Emin,Emax)


do i = 1 , size(qt%qauxvec)
    qt%qauxvec(i) = sum(qt%densities(:,i))
enddo
call qt%qsystem%save_data("densities.xml",array2d=qt%densities,array1d=qt%qauxvec)
call qt%save_system("system.xml")
call qt%destroy_system()
!call qsystem%destroy()
call tmpsystem%destroy()

print*,"Generating plots..."
!print*,"Plotting band structure..."
!call system("cd "//output_folder//"; ./plot_bands.py")
!print*,"Use Viewer program to see the structure and crated lead."
 contains

! ---------------------------------------------------------------------------
! This function decides if site A (called here atomA) with spin s1 has hoping
! to atom B with spin s2, and what is the value of the coupling.
! If there is no interaction between them returns false, otherwise true.
! ---------------------------------------------------------------------------
logical function connect_simple(atomA,atomB,coupling_val)
    use modatom
    implicit none
    type(qatom) :: atomA,atomB
    complex*16  :: coupling_val ! you must overwrite this variable

    ! default return value
    connect_simple = .false.
    coupling_val   = 0.0
    if( atomA%flag /= atomB%flag ) then
        connect_simple = .true.
        coupling_val = 1.0
    endif

    if( atomA%flag_aux0 == 1 .and. (atomA%flag == atomB%flag) ) then
        coupling_val = coupling_val + 0.00000001
        connect_simple = .true.
    endif

end function connect_simple
end program transporter



!! ------------------------------------------------------ !
!! Quantulaba - simple_graphene.f90 - Krzysztof Kolasinski 2015
!!
!! ------------------------------------------------------ !
!program transporter
! use modscatter
! use modsys
! use modshape
! use modunits
! implicit none
! type(qscatter) :: qt
! type(qshape) :: rect_shape
! character(*),parameter :: output_folder = ""
! integer :: i , j, k , N , m
! doubleprecision,parameter :: alpha30 =  30.0/180.0*M_PI
! doubleprecision,parameter :: vecs_armchair(2,2) =  (/  (/ 1.0D0,0.0D0 /) , (/ sin(alpha30) , cos(alpha30) /) /)
! doubleprecision,parameter :: atoms_armchair(2,2) =  (/  (/ 0.0D0,0.0D0 /) , (/ 0.0D0 , 1.0D0/sqrt(3.0) /) /)
!
! doubleprecision,parameter :: vecs_zigzag(2,2) =   (/  (/ (3.0/2.0)/sqrt(3.0D0)  ,0.5D0 /) , (/ -(3.0/2.0)/sqrt(3.0D0) , 0.5D0 /) /)
! doubleprecision,parameter :: atoms_zigzag(2,2) =  (/  (/ 0.0D0,0.0D0 /) , (/ 1.0D0/sqrt(3.0)  ,0.0D0  /) /)
!
! doubleprecision,parameter :: pos_offset(2) =  (/ -5.0D0,0.0D0 /)
! doubleprecision,parameter :: pos_offset_zigzag(2) =  (/ -5.0D0,-10.0D0 /)
! doubleprecision           :: atom_pos(3),lead_translation_vec(3),  Ef , zero_Vector(100)
! integer ,parameter        :: atom_A = 1 , atom_B = 2 , atom_C = 3
!doubleprecision :: xmax
! integer :: atom
! integer :: gindex(30,31)
!
!
!
! call qt%init_system()
!
! ! --------------------------------------------------------------------------
! ! ARMCHAIR test
! ! --------------------------------------------------------------------------
!! k = 0
!! atom_pos = 0
!! gindex   = 0
!! do i = 1 , size(gindex,1)
!! do j = 1 , size(gindex,2)
!!    do atom = atom_A , atom_B
!!    ! set atom position in space
!!    atom_pos(1:2) = atoms_armchair(:,atom) + (i-1) * vecs_armchair(:,1) +  (j-1) * vecs_armchair(:,2) + pos_offset
!!    ! cut some atoms to have rectangular flake
!!    if(atom_pos(1) > 0.2 .and. atom_pos(1) < 11.7 .and.  atom_pos(2) < 10.5 ) then
!!        call qt%qatom%init( (/ atom_pos(1) , atom_pos(2) , 0.0D0 /),flag=atom)
!!        call qt%qsystem%add_atom(qt%qatom)
!!        k = k + 1
!!        gindex(i,j) = k
!!    endif
!!    enddo ! end of atom loop
!!enddo
!!enddo
!!
!!qt%qnnbparam%distance   = 0.6
!!qt%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_DISTANCE
!!call qt%qsystem%make_lattice(connect,qt%qnnbparam)
!!! remove single bonds
!!do atom = 1 ,qt%qsystem%no_atoms
!!    if(qt%qsystem%atoms(atom)%no_bonds == 1)   qt%qsystem%atoms(atom)%bActive  = .false.
!!enddo
!
!
!!call qt%qsystem%make_lattice(connect,qt%qnnbparam)
!!call qt%qsystem%save_lattice(output_folder//"lattice.xml")
!!
!!! adding lead
!!call rect_shape%init_rect(SHAPE_RECTANGLE_XY,0.4D0,1.1D0,0.0D0,11.0D0)
!!call qt%add_lead(rect_shape,(/1.0D0,0.0D0,0.0D0/))
!!call qt%leads(1)%print_lead(output_folder//"lead.xml",qt%qsystem%atoms)
!!call qt%leads(1)%bands(output_folder//"bands.dat",-3.14D0,3.14D0,0.1D0,-15.0D0,15.0D0)
!
! ! --------------------------------------------------------------------------
! ! ZIGZAG test
! ! --------------------------------------------------------------------------
!xmax = 50.5
!k = 0
!gindex   = 0
!do i = -200 , 200
!do j = -200 , 200
!    do atom = atom_A , atom_B
!        ! set atom position in space
!        atom_pos(1:2) = atoms_zigzag(:,atom) + (i-1) * vecs_zigzag(:,1) +  (j-1) * vecs_zigzag(:,2) + pos_offset_zigzag
!        ! cut some atoms to have rectangular flake
!        if(atom_pos(1) > 0.2 .and. atom_pos(1) < xmax .and. &
!           atom_pos(2) > 0.0 .and. atom_pos(2) < 50.9 ) then
!            call qt%qatom%init( (/ atom_pos(1) , atom_pos(2) , 0.0D0 /),flag=atom)
!            if( atom_pos(1) > 20 .and. atom_pos(1) < 30.0 .and. atom_pos(2)  < 20.0 ) qt%qatom%bActive = .false.
!            if( atom_pos(1) > 20 .and. atom_pos(1) < 30.0 .and. atom_pos(2)  > 40.0 ) qt%qatom%bActive = .false.
!
!            call qt%qsystem%add_atom(qt%qatom)
!         !   k = k + 1
!        !    gindex(i,j) = k
!        endif
!    enddo ! end of atom loop
!enddo
!enddo
!
!qt%qnnbparam%distance   = 0.6
!qt%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_DISTANCE
!call qt%qsystem%make_lattice(qt%qnnbparam,c_default=connect)
!call qt%qsystem%save_lattice(output_folder//"lattice.xml")
!
!
!
!lead_translation_vec = (/ 3.0/sqrt(3.0D0) ,  0.0D0 , 0.0D0 /)
!call rect_shape%init_range_3d((/0.4D0,0.0D0,0D0/),lead_translation_vec)
!!call rect_shape%init_rect(SHAPE_RECTANGLE_XY,0.4D0,13.1D0,0.4D0,1.3D0)
!!lead_translation_vec = (/ 0.0D0 , 1.0D0 , 0.0D0 /)
!
!call qt%add_lead(rect_shape,lead_translation_vec)
!call qt%leads(1)%save_lead(output_folder//"lead.xml")
!call qt%leads(1)%bands(output_folder//"bands.dat",-3.14D0,3.14D0,0.1D0,-15.0D0,15.0D0)
!
!lead_translation_vec = (/ -3.0/sqrt(3.0D0) ,  0.0D0 , 0.0D0 /)
!call rect_shape%init_range_3d((/xmax+0.1D0,0.0D0,0D0/),lead_translation_vec)
!call qt%add_lead(rect_shape,lead_translation_vec)
!
!call qt%leads(2)%save_lead(output_folder//"lead2.xml")
!!call qt%leads(2)%bands(output_folder//"bands.dat",-3.14D0,3.14D0,0.1D0,-15.0D0,15.0D0)
!call qt%save_system("system.xml")
!!call qt%leads(1)%calculate_modes(0.3D0)
!Ef = 0.4;
!call qt%calculate_modes(Ef)
!call qt%solve(1,Ef)
!!print*,"asdasd"
!!do i = 1 , qt%qsystem
!
!k = qt%leads(1)%no_in_modes
!do i = 1 , qt%qsystem%no_atoms
!
!
!    if(qt%qsystem%atoms(i)%bActive) then
!        m = qt%qsystem%atoms(i)%globalIDs(1)
!!        print*,"m=",m,k
!        write(111,"(100e20.6)"),qt%qsystem%atoms(i)%atom_pos,sum(qt%densities(:,m)),qt%densities(:,m)
!    endif
!!do j = 1 , size(gindex,2)
!!    if(gindex(i,j)  == 0 ) cycle
!!    if(qt%qsystem%atoms(gindex(i,j))%bActive) then
!!    m = qt%qsystem%atoms(gindex(i,j))%globalIDs(1)
!!
!!    write(111,"(100e20.6)"),qt%qsystem%atoms(gindex(i,j))%atom_pos,sum(qt%densities(:,m)),qt%densities(:,m)
!!    else
!!    write(111,"(100e20.6)"),qt%qsystem%atoms(gindex(i,j))%atom_pos,zero_vector(1:k)
!!    endif
!!enddo
!!    write(111,*),""
!enddo
!call qt%destroy_system()
!
!
!print*,"Generating plots..."
!print*,"Plotting band structure..."
!!call system("cd "//output_folder//"; ./plot_bands.py")
!!print*,"Use Viewer program to see the structure and crated lead."
!
!contains
!
!logical function connect(atomA,atomB,s1,s2,coupling_val)
!    use modatom
!    implicit none
!    type(qatom) :: atomA,atomB
!    integer    :: s1,s2
!    complex*16 :: coupling_val
!    logical :: test
!    connect = .false.
!    ! In clean graphene atom from sublattice A always couples to atom from lattice B
!    test = .not.(atomA%flag == atomB%flag)
!    if(test) then
!        connect = .true.
!        coupling_val = 1.0D0
!    endif
!!    if(.not. test) then
!!        connect = .true.
!!        coupling_val = 0.0D0
!!    endif
!
!end function connect
!
!end program transporter



!
!
!! ------------------------------------------------------ !
!! Quantulaba - simple_sqlat.f90 - Krzysztof Kolasinski 2015
!!
!! We want to find few first eigenvalues of rectangular
!! system defined within effective mass Shroedinger
!! equation. Magnetic field will be included into account.
!! We assume that after the finite discretization nodes
!! of the mesh are separated by dx=5nm distance in
!! each direction.
!! ------------------------------------------------------ !
!
!program transporter
!
! use modunits ! unit conversion tools
! use modsys   ! eigen values
! use modlead  ! bandgap structure
! use modshape
! use modscatter
! implicit none
! character(*),parameter :: output_folder = "./"
!
! type(qshape)               :: lead_shape
!
! type(qscatter)             :: solver
! doubleprecision            :: zeros_vector(200),posA(3),posB(3),deltaR(3)
! doubleprecision            :: a_dx,a_Emin,a_Emax,a_Bz
! integer                    :: no_expected_states,flagA,flagB
! integer ,parameter         :: nx = 200
! integer ,parameter         :: ny = 100
! doubleprecision,parameter  :: dx = 1.0 ! [nm]
! integer , dimension(nx,ny) :: gindex ! converts local index (i,j) to global index
! integer                    :: i,j,k,m,l
! doubleprecision            :: lead_translation_vec(3),x,y,z,Ef,dens
! character(300) :: line
!
! ! Use atomic units in effective band model -> see modunit.f90 for more details
! call modunits_set_GaAs_params()
! a_dx = dx * L2LA ! convert it to atomic units
! a_Bz = BtoAtomicB(0.0D0) ! 1 Tesla to atomic units
!
! ! Initalize system
!
! call solver%init_system()
!
! ! ----------------------------------------------------------
! ! 1. Create mesh - loop over width and height of the lattice
! ! ----------------------------------------------------------
! k      = 0
! gindex = 0
! j = 1
! do i = 1 , nx
!! do j = 1 , ny
!    ! Initalize atom structure with position of the atom.
!    call solver%qatom%init((/ (i-1) * dx , (j-1) * dx + (i-1) * dx * 0.0 , 0.0 * dx /))
!    ! Add atom to the system.
!!    if( i > nx/2 .and. j < ny/2+20 ) solver%qatom%bActive = .false.
!!    if( i > nx/2+30 .and. j < ny/2+20 ) solver%qatom%bActive = .false.
!!    if( abs(i-nx/2) < 5 .and. j < 10 ) solver%qatom%bActive = .false.
!
!    call solver%qsystem%add_atom(solver%qatom)
!    k           = k + 1
!    gindex(i,j) = k
!! enddo
! enddo
!
!! ----------------------------------------------------------
!! 2. Construct logical connections between sites on the mesh.
!! ----------------------------------------------------------
!solver%qnnbparam%box = (/2*dx,2*dx,0.0D0/) ! do not search for the sites far than (-dx:+dx) direction
!call solver%qsystem%make_lattice(solver%qnnbparam,c_matrix=connect_matrix)
!
!! Save lattice to file to see if constructed system is OK!
!! Use plot_lattice.py to see the results.
!call solver%qsystem%save_lattice(output_folder//"lattice.xml")
!
!
!!! ! ----------------------------------------------------------
!!! ! 3. Find eigenvalues of the system.
!!! ! You must provide range of energy (Emin,Emax) and
!!! ! number of expected states in that region. Eigen problem
!!! ! is solved with FEAST library.
!!! ! ----------------------------------------------------------
!! no_expected_states = 100
!! ! we choose small window
!! a_Emin =0.0
!! a_Emax =0.2
!! call qsystem%calc_eigenproblem(a_Emin,a_Emax,no_expected_states)
!!
!! ! Write founded eigenstates to file if there are states in the range (Emin,Emax)
!! if(qsystem%no_eigenvalues > 0) then
!! zeros_vector = 0
!! open(unit=222,file=output_folder//"eigenvecs.dat")
!! do i = 1 , nx
!! do j = 1 , ny
!!    ! Remapping to graphen lattice
!!    x = qsystem%atoms(gindex(i,j))%atom_pos(1) !+ qsystem%atoms(gindex(i,j))%atom_pos(2)*sin(alpha30))/dx
!!    y = qsystem%atoms(gindex(i,j))%atom_pos(2) !- sin(alpha30)*qsystem%atoms(gindex(i,j))%atom_pos(1))/dx
!!    z = 0
!!    if(qsystem%atoms(gindex(i,j))%bActive)then
!!    ! get unique ID of current site
!!
!!    m = qsystem%atoms(gindex(i,j))%globalIDs(1)
!!    write(222,"(500e20.6)"),x,y,z,abs(qsystem%eigenvecs(m,:))**2
!!    else
!!    ! fill with zeros empty spaces
!!    write(222,"(500e20.6)"),x,y,z,zeros_vector(1:qsystem%no_eigenvalues)
!!    endif
!! enddo
!!    write(222,*),"" ! GNUPLOT necessary line break
!! enddo
!! close(222)
!! print*,qsystem%eigenvals
!! endif ! end of if are there eigenstates?
!
!! ----------------------------------------------------------
!! 4. Use generated mesh to calculate the band structure
!! in the region of homogenous lead.
!! ----------------------------------------------------------
! ! Lead needs to know where it is (lead_shape) and using this information it will
! ! create propper matrices and connections using list of atoms
!lead_translation_vec = (/  dx , 0.0D0 , 0.0D0 /)
!call lead_shape%init_range_3d((/ -0.5*dx , 0.0D0 , 0.0D0 /),lead_translation_vec)
!call solver%add_lead(lead_shape,lead_translation_vec)
!call solver%leads(1)%save_lead(output_folder//"lead1.xml")
!
!!a_Emin =-10.0
!!a_Emax = 10.0
!!call solver%leads(1)%bands(output_folder//"bands.dat",-M_PI,+M_PI,M_PI/160.0,a_Emin,a_Emax)
!
!lead_translation_vec = (/ -dx , 0.0D0 , 0.0D0 /)
!call lead_shape%init_range_3d((/ (+0.5+nx-1)*dx , 0.0D0 , 0.0D0 /),lead_translation_vec)
!call solver%add_lead(lead_shape,lead_translation_vec)
!call solver%leads(2)%save_lead(output_folder//"lead2.xml")
!
!call solver%save_system("system.xml")
!Ef = 0.1
!call solver%calculate_modes(Ef);
!
!
!call solver%solve(1,Ef);
!
!
!
!l = solver%leads(1)%no_in_modes
!
!j = 1
!do i = 1 , nx
!!do j = 1 , ny
!    if(solver%qsystem%atoms(gindex(i,j))%bActive) then
!    m = solver%qsystem%atoms(gindex(i,j))%globalIDs(1)
!    k = gindex(i,j)
!    dens = sum(solver%densities(:,m))
!    write(111,"(100e20.6)"),solver%qsystem%atoms(k)%atom_pos(1:2),dens,solver%densities(:,m)
!    else
!    write(111,"(100e20.6)"),solver%qsystem%atoms(gindex(i,j))%atom_pos(1:2),0.0D0,zeros_vector(1:l)
!    endif
!!enddo
!!    write(111,*),""
!enddo
!
!!lead_translation_vec = (/  -dx , 0.0D0 , 0.0D0 /)
!!
!!call lead_shape%init_range_3d((/ (nx+0.5-1)*dx , 0.0D0 , 0.0D0 /),lead_translation_vec)
!!call lead%init_lead(lead_shape,lead_translation_vec,qsystem%atoms)
!!call lead%print_lead(output_folder//"lead.xml",qsystem%atoms)
!!call lead%calculate_modes(0.08D0)
!
!
!
!!call lead%print_lead(output_folder//"lead.xml",qsystem%atoms)
!!a_Emin =-10.0
!!a_Emax = 10.0
!!call lead%bands(output_folder//"bands.dat",-M_PI,+M_PI,M_PI/160.0,a_Emin,a_Emax)
!
!! ! ----------------------------------------------------------
!! ! X. Clean memory...
!! ! ----------------------------------------------------------
! call solver%destroy_system()
!!!!
!!!
!!!
!
!!! print*,"Generating plots..."
!! print*,"Plotting band structure..."
!! call system("cd "//output_folder//"; ./plot_bands.py")
!! print*,"Plotting eigenvectors..."
!! call system("cd "//output_folder//"; ./plot_eigenvecs.py")
!!! print*,"Plotting lattice..."
!!! call system("cd "//output_folder//"; ./plot_lattice.py")
!!! print*,"Plotting lead..."
!!! call system("cd "//output_folder//"; ./plot_lead.py")
! contains
!
!! ---------------------------------------------------------------------------
!! Implement coupling
!! Taken from Kwant tutorial: http://kwant-project.org/doc/1.0/tutorial/tutorial5
!! ---------------------------------------------------------------------------
!logical function connect_matrix(atomA,atomB,coupling_mat)
!    use modatom
!    implicit none
!    type(qatom) :: atomA,atomB
!
!    complex*16  :: coupling_mat(:,:) ! you must overwrite this variable
!    ! local variables
!    integer         :: xdiff,ydiff
!    doubleprecision :: dydiff,dxdiff,t0,y,x
!
!    ! Calculate distance between atoms in units of dx.
!    dxdiff = -(atomA%atom_pos(1)-atomB%atom_pos(1))/dx
!    dydiff = -(atomA%atom_pos(2)-atomB%atom_pos(2))/dx
!    ! Convert it to integers
!    xdiff = NINT(dxdiff)
!    ydiff = NINT(dydiff)
!
!    x = atomA%atom_pos(1)
!
!    ! default return value
!    connect_matrix = .false.
!    coupling_mat   = 0.0
!    ! We assume that there is no spin taken into account so in our
!    ! case is always s1 = 1 and s2 = 1, thus s1 and s2 are not used
!    ! here.
!    t0    = 1.0
!
!
!    if( xdiff == 0 .and. ydiff == 0 ) then
!        connect_matrix      = .true.
!        coupling_mat = 2*t0 + 0.1*t0*exp(-0.04*(x-150.0)**2) + 0.1*t0*exp(-0.04*(x-50.0)**2)
!    else if( abs(xdiff) ==  1 .and. ydiff == 0 ) then
!        connect_matrix = .true.
!        y = (atomA%atom_pos(2) - ny/2*dx ) * L2LA
!        coupling_mat = -t0 * EXP(II*xdiff*y*0.00)
!    else if( xdiff ==  0 .and. abs(ydiff) == 1 ) then
!        connect_matrix = .true.
!        coupling_mat = -t0
!    endif
!
!end function connect_matrix
!
!end program transporter
