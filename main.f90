! ------------------------------------------------------ !
<<<<<<< HEAD
! Quantulaba - benchmark.f90 - Krzysztof Kolasinski 2016
!
! Test of speed and accuracy
=======
! Quantulaba - simple_graphene.f90 - Krzysztof Kolasinski 2015
!
>>>>>>> ab1ee0d20986ddd68c0c0bfa4bb7a9a073c82237
! ------------------------------------------------------ !

program graphene2
 use modscatter
 use modsys
 use modshape
 use modunits
 implicit none
<<<<<<< HEAD
 character(*),parameter :: output_folder = "plots/"
 type(qscatter)             :: qt
 type(qshape)               :: lead_shape
 doubleprecision            :: zeros_vector(200)
 doubleprecision            :: a_dx,a_Emin,a_Emax,a_Bz
 integer                    :: no_expected_states
 integer          :: nx = 50
 integer          :: ny = 50
 doubleprecision,parameter  :: dx = 5.0 ! [nm]

 integer :: i,j,k
 doubleprecision            :: lead_translation_vec(3),Ef,Bz
 doubleprecision :: T_AUTO,T_GGEV,T_SCHUR,T_EXACT ! transmissions
 doubleprecision :: Tm_AUTO,Tm_GGEV,Tm_SCHUR ! timers

! Use atomic units in effective band model -> see modunit.f90 for more details
call modunits_set_GaAs_params()
a_dx = dx * L2LA ! convert it to atomic units
a_Bz = BtoAtomicB(0.45D0) ! 1 Tesla to atomic units


ny = 20
call create_system(ny,3)

Ef = 5.0/E0/1000.0
call qt%calculate_modes(Ef)
call qt%solve(1,Ef)
! Save calculated electron density to file
do i = 1 , size(qt%qauxvec)
    qt%qauxvec(i) = sum(qt%densities(:,i))
enddo
call qt%qsystem%save_data(output_folder//"densities.xml",array2d=qt%densities,array1d=qt%qauxvec)

=======
 type(qscatter) :: qt
 type(qshape) :: lead_area


 character(*),parameter    :: output_folder = "plots/"
 doubleprecision,parameter :: alpha30 =  30.0/180.0*M_PI
 doubleprecision,parameter :: vecs_armchair(2,2) =  (/  (/ 1.0D0,0.0D0 /) , (/ sin(alpha30) , cos(alpha30) /) /)
 doubleprecision,parameter :: atoms_armchair(2,2) =  (/  (/ 0.0D0,0.0D0 /) , (/ 0.0D0 , 1.0D0/sqrt(3.0) /) /)
 doubleprecision,parameter :: pos_offset(2) =  (/ -5.0D0,0.0D0 /)
 doubleprecision           :: atom_pos(3),lead_translation_vec(3),range_base(3),range_dir(3),Ef
 integer,parameter         :: atom_A = 1 , atom_B = 2
 integer           :: atom,i,j
 integer,parameter :: nx = 50  , ny = 5

 call qt%init_system()
!QSYS_FORCE_SCHUR_DECOMPOSITION = .true.
! --------------------------------------------------------------------------
! Graphene test
! 1. Generate atoms positions:
! --------------------------------------------------------------------------
do i = 1 , nx
do j = 1 , ny
    do atom = atom_A , atom_B
    ! set atom position in space
        atom_pos(1:2)   = atoms_armchair(:,atom) + &
                   (i-1) * vecs_armchair(:,1) +  &
                   (j-1) * vecs_armchair(:,2)
        call qt%qatom%init( (/ atom_pos(1) , atom_pos(2) , 0.0D0 /),flag=atom)
        call qt%qsystem%add_atom(qt%qatom)
    enddo ! end of atom loop
enddo
enddo
! --------------------------------------------------------------------------
! 2. Initiate couplings between atoms: We use c_simple connection type
!    i.e. without including spin degree of freedom per atom
! --------------------------------------------------------------------------
qt%qnnbparam%distance   = 0.6
qt%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_DISTANCE
call qt%qsystem%make_lattice(qt%qnnbparam,c_simple=connect)
>>>>>>> ab1ee0d20986ddd68c0c0bfa4bb7a9a073c82237

! remove single bonds if necessary
do atom = 1 ,qt%qsystem%no_atoms
    if(qt%qsystem%atoms(atom)%no_bonds == 2)   qt%qsystem%atoms(atom)%bActive  = .false.
enddo
! 3. Calculate coupling again.
call qt%qsystem%make_lattice(qt%qnnbparam,c_simple=connect)

<<<<<<< HEAD
! Enable/Disable part of tests
if(.false.) then
! -------------------------------------------------------------------
! Test 1
! -------------------------------------------------------------------
ny = 50
call create_system(ny,1)

Ef = 5.0/E0/1000.0
call qt%calculate_modes(Ef)
call qt%solve(1,Ef)
! Save calculated electron density to file
do i = 1 , size(qt%qauxvec)
    qt%qauxvec(i) = sum(qt%densities(:,i))
enddo
call qt%qsystem%save_data(output_folder//"densities.xml",array2d=qt%densities,array1d=qt%qauxvec)



call qt%save_system(output_folder//"system_Test1.xml")
open(unit=111,file=output_folder//"Test1.dat")
print*,"Performing energy Test1..."
QSYS_DEBUG_LEVEL = 0 ! show more info
do Ef = 0.0 , 8.0/E0/1000.0 , 8.0/E0/1000.0/80
    ! Update hamiltonian elemenents value
    call qt%qsystem%update_lattice(c_simple=connect)

    QSYS_FORCE_ZGGEV_TO_FIND_MODES = .false.
    QSYS_FORCE_SCHUR_DECOMPOSITION = .false.
    Tm_AUTO = get_clock()
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)
    T_AUTO  =sum(qt%Tn(:))+sum(qt%Rn(:))
    T_EXACT = qt%leads(1)%no_in_modes
    Tm_AUTO = get_clock() - Tm_AUTO

    QSYS_FORCE_ZGGEV_TO_FIND_MODES = .true.
    QSYS_FORCE_SCHUR_DECOMPOSITION = .false.
    Tm_GGEV = get_clock()
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)
    T_GGEV  = sum(qt%Tn(:))+sum(qt%Rn(:))
    Tm_GGEV = get_clock() - Tm_GGEV

    QSYS_FORCE_ZGGEV_TO_FIND_MODES = .false.
    QSYS_FORCE_SCHUR_DECOMPOSITION = .true.
    Tm_SCHUR = get_clock()
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)
    T_SCHUR  =sum(qt%Tn(:))+sum(qt%Rn(:))
    Tm_SCHUR = get_clock() - Tm_SCHUR
    print*,Ef*1000*E0,sum(qt%Tn(:))
    write(111,"(100e20.10)"),Ef*1000*E0,T_EXACT,T_AUTO,T_GGEV,T_SCHUR,Tm_AUTO,Tm_GGEV,Tm_SCHUR
enddo
close(111)

call system("cd "//output_folder//"; python plot_Tests.py 1")


! -------------------------------------------------------------------
! Test 2
! -------------------------------------------------------------------
ny = 80
call create_system(ny,1)


call qt%save_system(output_folder//"system_Test2.xml")
open(unit=111,file=output_folder//"Test2.dat")
print*,"Performing energy Test2..."
QSYS_DEBUG_LEVEL = 0 ! show more info
do Ef = 0.0 , 8.0/E0/1000.0 , 8.0/E0/1000.0/80
    ! Update hamiltonian elemenents value
    call qt%qsystem%update_lattice(c_simple=connect)

    QSYS_FORCE_ZGGEV_TO_FIND_MODES = .false.
    QSYS_FORCE_SCHUR_DECOMPOSITION = .false.
    Tm_AUTO = get_clock()
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)
    T_AUTO  =sum(qt%Tn(:))+sum(qt%Rn(:))
    T_EXACT = qt%leads(1)%no_in_modes
    Tm_AUTO = get_clock() - Tm_AUTO

    QSYS_FORCE_ZGGEV_TO_FIND_MODES = .true.
    QSYS_FORCE_SCHUR_DECOMPOSITION = .false.
    Tm_GGEV = get_clock()
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)
    T_GGEV  = sum(qt%Tn(:))+sum(qt%Rn(:))
    Tm_GGEV = get_clock() - Tm_GGEV

    QSYS_FORCE_ZGGEV_TO_FIND_MODES = .false.
    QSYS_FORCE_SCHUR_DECOMPOSITION = .true.
    Tm_SCHUR = get_clock()
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)
    T_SCHUR  =sum(qt%Tn(:))+sum(qt%Rn(:))
    Tm_SCHUR = get_clock() - Tm_SCHUR
    print*,Ef*1000*E0,sum(qt%Tn(:))
    write(111,"(100e20.10)"),Ef*1000*E0,T_EXACT,T_AUTO,T_GGEV,T_SCHUR,Tm_AUTO,Tm_GGEV,Tm_SCHUR
enddo
close(111)

call system("cd "//output_folder//"; python plot_Tests.py 2")



! -------------------------------------------------------------------
! Test 3
! -------------------------------------------------------------------
ny = 100
call create_system(ny,1)


call qt%save_system(output_folder//"system_Test3.xml")
open(unit=111,file=output_folder//"Test3.dat")
print*,"Performing energy Test3..."
QSYS_DEBUG_LEVEL = 0 ! show more info
do Ef = 0.0 , 8.0/E0/1000.0 , 8.0/E0/1000.0/80
    ! Update hamiltonian elemenents value
    call qt%qsystem%update_lattice(c_simple=connect)

    QSYS_FORCE_ZGGEV_TO_FIND_MODES = .false.
    QSYS_FORCE_SCHUR_DECOMPOSITION = .false.
    Tm_AUTO = get_clock()
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)
    T_AUTO  =sum(qt%Tn(:))+sum(qt%Rn(:))
    T_EXACT = qt%leads(1)%no_in_modes
    Tm_AUTO = get_clock() - Tm_AUTO

    QSYS_FORCE_ZGGEV_TO_FIND_MODES = .true.
    QSYS_FORCE_SCHUR_DECOMPOSITION = .false.
    Tm_GGEV = get_clock()
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)
    T_GGEV  = sum(qt%Tn(:))+sum(qt%Rn(:))
    Tm_GGEV = get_clock() - Tm_GGEV

    QSYS_FORCE_ZGGEV_TO_FIND_MODES = .false.
    QSYS_FORCE_SCHUR_DECOMPOSITION = .true.
    Tm_SCHUR = get_clock()
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)
    T_SCHUR  =sum(qt%Tn(:))+sum(qt%Rn(:))
    Tm_SCHUR = get_clock() - Tm_SCHUR
    print*,Ef*1000*E0,sum(qt%Tn(:))
    write(111,"(100e20.10)"),Ef*1000*E0,T_EXACT,T_AUTO,T_GGEV,T_SCHUR,Tm_AUTO,Tm_GGEV,Tm_SCHUR
enddo
close(111)

call system("cd "//output_folder//"; python plot_Tests.py 3")
endif


if(.false.) then
! -------------------------------------------------------------------
! Test 4
! -------------------------------------------------------------------
ny = 100
nx = 30
call create_system(ny,2)
a_Bz = BtoAtomicB(0.45D0) ! 1 Tesla to atomic units
Ef = 5.0/E0/1000.0
call qt%calculate_modes(Ef)
=======
! --------------------------------------------------------------------------
! 3. Add translational lead to the system. We use range_3d object to define
!    area where the lead is located in the system.
! --------------------------------------------------------------------------
range_base = (/-0.5,0.0,0.0 /) ! initial position of range
range_dir  =  0.8*(/cos(alpha30),-sin(alpha30),0.0D0/) ! direction of the range (lenght contains distance)
call lead_area%init_range_3d(range_base,range_dir)
call qt%add_lead(lead_area,(/1.0D0,0.0D0,0.0D0/))
!call qt%leads(1)%bands(output_folder//"bands.dat",-3.14D0,3.14D0,0.1D0,-15.0D0,15.0D0)

range_base = (/49.4,0.0,0.0 /) ! initial position of range
call lead_area%init_range_3d(range_base,-range_dir)
call qt%add_lead(lead_area,(/-1.0D0,0.0D0,0.0D0/))


QSYS_DEBUG_LEVEL = 0

Ef = 2.0

call qt%calculate_modes(Ef)
>>>>>>> ab1ee0d20986ddd68c0c0bfa4bb7a9a073c82237
call qt%solve(1,Ef)

! Save calculated electron density to file
do i = 1 , size(qt%qauxvec)
    qt%qauxvec(i) = sum(qt%densities(:,i))
enddo
call qt%qsystem%save_data(output_folder//"densities.xml",array2d=qt%densities,array1d=qt%qauxvec)
<<<<<<< HEAD



call qt%save_system(output_folder//"system_Test4.xml")
open(unit=111,file=output_folder//"Test4.dat")
print*,"Performing energy Test4..."
QSYS_DEBUG_LEVEL = 0 ! show more info
do Ef = 0.0 , 8.0/E0/1000.0 , 8.0/E0/1000.0/40
    ! Update hamiltonian elemenents value
    call qt%qsystem%update_lattice(c_simple=connect)

    QSYS_FORCE_ZGGEV_TO_FIND_MODES = .false.
    QSYS_FORCE_SCHUR_DECOMPOSITION = .false.
    Tm_AUTO = get_clock()
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)
    T_AUTO  =sum(qt%Tn(:))+sum(qt%Rn(:))
    T_EXACT = qt%leads(1)%no_in_modes
    Tm_AUTO = get_clock() - Tm_AUTO

    QSYS_FORCE_ZGGEV_TO_FIND_MODES = .true.
    QSYS_FORCE_SCHUR_DECOMPOSITION = .false.
    Tm_GGEV = get_clock()
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)
    T_GGEV  = sum(qt%Tn(:))+sum(qt%Rn(:))
    Tm_GGEV = get_clock() - Tm_GGEV

    QSYS_FORCE_ZGGEV_TO_FIND_MODES = .false.
    QSYS_FORCE_SCHUR_DECOMPOSITION = .true.
    Tm_SCHUR = get_clock()
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)
    T_SCHUR  =sum(qt%Tn(:))+sum(qt%Rn(:))
    Tm_SCHUR = get_clock() - Tm_SCHUR
    print*,Ef*1000*E0,sum(qt%Tn(:))
    write(111,"(100e20.10)"),Ef*1000*E0,T_EXACT,T_AUTO,T_GGEV,T_SCHUR,Tm_AUTO,Tm_GGEV,Tm_SCHUR
enddo
close(111)
call system("cd "//output_folder//"; python plot_Tests.py 4")
endif ! end of if ennable/disable second test group


! ----------------------------------------------------------
! X. Clean memory...
! ----------------------------------------------------------
call lead_shape%destroy_shape()
call qt%destroy_system()
print*,"Use Viewer program to see the structure and created leads."
contains

! ---------------------------------------------------------------------------
! This function decides if site A (called here atomA)  has hoping
! to atom B, and what is the value of the coupling.
! If there is no interaction between them returns false, otherwise true.
! ---------------------------------------------------------------------------
=======
call qt%save_system(output_folder//"system.xml")
!
!print*,"Performing energy scan..."
!open(unit=111,file=output_folder//"T.dat")
!!QSYS_DEBUG_LEVEL = 1 ! show more info
!do Ef = -3.0 , 3.025 , 0.025
!    ! Update hamiltonian elemenents value
!    call qt%qsystem%update_lattice(c_simple=connect)
!    call qt%calculate_modes(Ef)
!    call qt%solve(1,Ef)
!
!    print*,"Energy:",Ef
!    write(111,"(100f20.6)"),Ef,sum(qt%Tn(:))
!enddo
!close(111)
!
!print*,"Generating plots..."
!print*,"Plotting band structure..."
!call system("cd "//output_folder//"; ./plot_bands.py")
!print*,"Plotting Transmission..."
!call system("cd "//output_folder//"; ./plot_T.py")
!print*,"Use Viewer program to see the structure and created leads."

call qt%destroy_system()
contains

>>>>>>> ab1ee0d20986ddd68c0c0bfa4bb7a9a073c82237
logical function connect(atomA,atomB,coupling_val)
    use modcommons
    implicit none
    type(qatom) :: atomA,atomB
<<<<<<< HEAD

    complex*16  :: coupling_val ! you must overwrite this variable
    ! local variables
    integer         :: xdiff,ydiff
    doubleprecision :: dydiff,dxdiff,t0,y

    ! Calculate distance between atoms in units of dx.
    dxdiff = (atomA%atom_pos(1)-atomB%atom_pos(1))/dx
    dydiff = (atomA%atom_pos(2)-atomB%atom_pos(2))/dx


    ! Convert it to integers
    xdiff = NINT(dxdiff)
    ydiff = NINT(dydiff)

    ! default return value
    connect = .false.

    ! hoping parameter
    t0 = 1/(2*m_eff*a_dx**2)

    if( xdiff == 0 .and. ydiff == 0 ) then
        connect      = .true.
        coupling_val = 4*t0
    else if( abs(xdiff) ==  1 .and. ydiff == 0 ) then
        connect = .true.
        ! magnetic Field enters the hamiltonian by phase: EXP(i*DX*Bx*y)
        y = (atomA%atom_pos(2) - ny/2*dx ) * L2LA ! convert to atomic units
        coupling_val = -t0 * EXP(II*xdiff*a_dx*a_Bz*y)
    else if( xdiff ==  0 .and. abs(ydiff) == 1 ) then
        connect = .true.
        coupling_val = -t0
    endif
end function connect

! Create system of with n_y and unit cell lenght l_uc
! The system will contain little hole in the middle
subroutine create_system(n_y,l_uc)
    integer :: n_y , l_uc
    call lead_shape%destroy_shape()
    call qt%destroy_system()
    call qt%init_system()
    ! ----------------------------------------------------------
    ! 1. Create mesh - loop over width and height of the lattice
    ! ----------------------------------------------------------
    do i = 1 , nx
    do j = 1 , n_y
        call qt%qatom%init((/ (i-1) * dx , (j-1) * dx  , 0.0 * dx /))
        if( sqrt(abs(i-nx/2.0)**2) + sqrt(abs(j-n_y/2.0)**2)  < 4. ) then
        else
        call qt%qsystem%add_atom(qt%qatom)
        endif
    enddo
    enddo
    ! ----------------------------------------------------------
    ! 2. Construct logical connections between sites on the mesh.
    ! ----------------------------------------------------------
    qt%qnnbparam%box = (/2*dx,2*dx,0.0D0/) ! do not search for the sites far than (-dx:+dx) direction
    ! Setup connections between sites with provided by you function "connect", see below for example.
    call qt%qsystem%make_lattice(qt%qnnbparam,c_simple=connect)
    ! ----------------------------------------------------------
    ! 3. Use generated mesh to calculate the band structure
    ! in the region of homogenous lead.
    ! ----------------------------------------------------------
    lead_translation_vec = (/  l_uc*dx , 0.0D0 , 0.0D0 /)
    call lead_shape%init_range_3d((/-0.5D0*dx,0.0D0,0.0D0/),lead_translation_vec)
    call qt%add_lead(lead_shape,lead_translation_vec)
end subroutine create_system

end program transporter
=======
    complex*16 :: coupling_val
    logical :: test
    connect = .false.

    ! In clean graphene atom from sublattice A always couples to atom from lattice B
    test = .not.(atomA%flag == atomB%flag)
    if(test) then
        connect      = .true.
        coupling_val = 1.1D0
    endif
end function connect
end program graphene2
>>>>>>> ab1ee0d20986ddd68c0c0bfa4bb7a9a073c82237
