! ------------------------------------------------------ !
! Bubel - simple_graphene.f90 - Krzysztof Kolasinski 2015
!
! ------------------------------------------------------ !

program graphene2
 use modscatter
 use modsys
 use modshape
 use modunits
 implicit none
 type(qscatter) :: qt
 type(qshape) :: lead_area


 character(*),parameter    :: output_folder = "simple_graphene2_output/"
 doubleprecision,parameter :: alpha30 =  30.0/180.0*M_PI
 doubleprecision,parameter :: vecs_armchair(2,2) =  (/  (/ 1.0D0,0.0D0 /) , (/ sin(alpha30) , cos(alpha30) /) /)
 doubleprecision,parameter :: atoms_armchair(2,2) =  (/  (/ 0.0D0,0.0D0 /) , (/ 0.0D0 , 1.0D0/sqrt(3.0) /) /)
 doubleprecision,parameter :: pos_offset(2) =  (/ -5.0D0,0.0D0 /)
 doubleprecision           :: atom_pos(3),lead_translation_vec(3),range_base(3),range_dir(3),Ef
 integer,parameter         :: atom_A = 1 , atom_B = 2
 integer           :: atom,i,j
 integer,parameter :: nx = 50  , ny = 10

 call qt%init_system()
QSYS_FORCE_SCHUR_DECOMPOSITION = .true.
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

! remove single bonds if necessary
do atom = 1 ,qt%qsystem%no_atoms
    if(qt%qsystem%atoms(atom)%no_bonds == 2)   qt%qsystem%atoms(atom)%bActive  = .false.
enddo
! 3. Calculate coupling again.
call qt%qsystem%make_lattice(qt%qnnbparam,c_simple=connect)

! --------------------------------------------------------------------------
! 3. Add translational lead to the system. We use range_3d object to define
!    area where the lead is located in the system.
! --------------------------------------------------------------------------
range_base = (/-0.5,0.0,0.0 /) ! initial position of range
range_dir  =  0.8*(/cos(alpha30),-sin(alpha30),0.0D0/) ! direction of the range (lenght contains distance)
call lead_area%init_range_3d(range_base,range_dir)
call qt%add_lead(lead_area,(/1.0D0,0.0D0,0.0D0/))
call qt%leads(1)%bands(output_folder//"bands.dat",-3.14D0,3.14D0,0.1D0,-15.0D0,15.0D0)

range_base = (/49.4,0.0,0.0 /) ! initial position of range
call lead_area%init_range_3d(range_base,-range_dir)
call qt%add_lead(lead_area,(/-1.0D0,0.0D0,0.0D0/))



call qt%save_system(output_folder//"system.xml")
Ef = 0.1

call qt%calculate_modes(Ef)
call qt%solve(1,Ef)

! Save calculated electron density to file
do i = 1 , size(qt%qauxvec)
    qt%qauxvec(i) = sum(qt%densities(:,i))
enddo
call qt%qsystem%save_data(output_folder//"densities.xml",array2d=qt%densities,array1d=qt%qauxvec)


print*,"Performing energy scan..."
open(unit=111,file=output_folder//"T.dat")
!QSYS_DEBUG_LEVEL = 1 ! show more info
do Ef = -3.0 , 3.025 , 0.025
    ! Update hamiltonian elemenents value
    call qt%qsystem%update_lattice(c_simple=connect)
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)

    print*,"Energy:",Ef
    write(111,"(100f20.6)"),Ef,sum(qt%Tn(:))
enddo
close(111)

print*,"Generating plots..."
print*,"Plotting band structure..."
call system("cd "//output_folder//"; ./plot_bands.py")
print*,"Plotting Transmission..."
call system("cd "//output_folder//"; ./plot_T.py")
print*,"Use Viewer program to see the structure and created leads."

call qt%destroy_system()
contains

logical function connect(atomA,atomB,coupling_val)
    use modcommons
    implicit none
    type(qatom) :: atomA,atomB
    complex*16 :: coupling_val
    logical :: test
    connect = .false.

    ! In clean graphene atom from sublattice A always couples to atom from lattice B
    test = .not.(atomA%flag == atomB%flag)
    if(test) then
        connect      = .true.
        coupling_val = 1.0D0
    endif
end function connect
end program graphene2
