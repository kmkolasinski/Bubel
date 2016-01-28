! ------------------------------------------------------ !
! Quantulaba - simple_graphene.f90 - Krzysztof Kolasinski 2015
!
! ------------------------------------------------------ !
program graphene
 use modscatter
 use modsys
 use modshape
 use modunits
 implicit none
 type(qscatter) :: qt
 type(qshape) :: rect_shape
 character(*),parameter :: output_folder = "simple_graphene_output/"
 integer :: i , j, k , N

 doubleprecision,parameter :: alpha30 =  30.0/180.0*M_PI
 doubleprecision,parameter :: vecs_armchair(2,2) =  (/  (/ 1.0D0,0.0D0 /) , (/ sin(alpha30) , cos(alpha30) /) /)
 doubleprecision,parameter :: atoms_armchair(2,2) =  (/  (/ 0.0D0,0.0D0 /) , (/ 0.0D0 , 1.0D0/sqrt(3.0) /) /)

 doubleprecision,parameter :: vecs_zigzag(2,2) =   (/  (/ (3.0/2.0)/sqrt(3.0D0)  ,0.5D0 /) , (/ -(3.0/2.0)/sqrt(3.0D0) , 0.5D0 /) /)
 doubleprecision,parameter :: atoms_zigzag(2,2) =  (/  (/ 0.0D0,0.0D0 /) , (/ 1.0D0/sqrt(3.0)  ,0.0D0  /) /)

 doubleprecision,parameter :: pos_offset(2) =  (/ -5.0D0,0.0D0 /)
 doubleprecision,parameter :: pos_offset_zigzag(2) =  (/ -5.0D0,-10.0D0 /)
 doubleprecision           :: atom_pos(3),lead_translation_vec(3),lead_start_pos(3),Ef
 integer ,parameter        :: atom_A = 1 , atom_B = 2 , atom_C = 3

 integer :: atom

call qt%init_system()
!--------------------------------------------------------------------
! ZIGZAG test
! --------------------------------------------------------------------------
QSYS_FORCE_SCHUR_DECOMPOSITION  = .true. ! use schur method to calculate modes which is more stable


do i = 1 , 100
do j = 1 , 500
    do atom = atom_A , atom_B
        ! set atom position in space
        atom_pos(1:2) = atoms_zigzag(:,atom) + &
                    (i-1) * vecs_zigzag(:,1) + &
                    (j-1) * vecs_zigzag(:,2) + pos_offset_zigzag
        ! cut some atoms to have rectangular flake
        if(atom_pos(1) > 0.2 .and. atom_pos(1) < 12.5 .and. &
           atom_pos(2) > 0.0 .and. atom_pos(2) < 70.9 ) then
            call qt%qatom%init( (/ atom_pos(1) , atom_pos(2) , 0.0D0 /),flag=atom)
            call qt%qsystem%add_atom(qt%qatom)
        endif
    enddo ! end of atom loop
enddo
enddo

qt%qnnbparam%distance   = 0.6
qt%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_DISTANCE
call qt%qsystem%make_lattice(qt%qnnbparam,c_simple=connect)

! Add first lead using rectangle area
call rect_shape%init_rect(SHAPE_RECTANGLE_XY,0.4D0,13.1D0,0.4D0,1.3D0)
lead_translation_vec = (/ 0.0D0 , 1.0D0 , 0.0D0 /)

call qt%add_lead(rect_shape,lead_translation_vec)
call qt%leads(1)%bands(output_folder//"bands.dat",-3.14D0,3.14D0,0.1D0,-15.0D0,15.0D0)

! Another lead using different approach.
lead_start_pos = (/0.0,70.9,0.0/) ! Start vector
lead_translation_vec = (/ 0.0D0 ,-1.0D0 , 0.0D0 /) ! Lead - Unit cell direction and width
call rect_shape%init_range_3d(lead_start_pos,lead_translation_vec)
call qt%add_lead(rect_shape,lead_translation_vec)

call qt%save_system(output_folder//"system.xml")

Ef = 2.0
QSYS_DEBUG_LEVEL = 1 ! show more info
! Force ZGGEV for graphen when looking for eigen modes
QSYS_USE_ZGGEV_TO_FIND_MODES = .true.

call qt%calculate_modes(Ef)
call qt%solve(1,Ef)
! Save calculated electron density to file
do i = 1 , size(qt%qauxvec)
    qt%qauxvec(i) = sum(qt%densities(:,i))
enddo
call qt%qsystem%save_data(output_folder//"densities.xml",array2d=qt%densities,array1d=qt%qauxvec)


print*,"Performing energy scan..."
open(unit=111,file=output_folder//"T.dat")
QSYS_DEBUG_LEVEL = 0 ! show more info
do Ef = -3.0+0.0001 , 3.025 , 0.025
    ! Update hamiltonian elemenents value
    call qt%qsystem%update_lattice(c_simple=connect)
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)

    print*,"Energy:",Ef
    write(111,"(100f20.6)"),Ef,sum(qt%Tn(:))
enddo
close(111)

call qt%destroy_system()


print*,"Generating plots..."
print*,"Plotting band structure..."
call system("cd "//output_folder//"; ./plot_bands.py")
print*,"Plotting Transmission..."
call system("cd "//output_folder//"; ./plot_T.py")
print*,"Use Viewer program to see the structure and created leads."

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
        connect = .true.
        coupling_val = 1.0D0
    endif
end function connect

end program graphene
