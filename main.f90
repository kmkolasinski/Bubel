program transporter
 use modscatter
 use modsys
 use modshape
 use modunits
 implicit none
 type(qscatter) :: qt
 type(qshape) :: rect_shape
 doubleprecision :: quad_positions(4,2)
 integer :: i , j, k , N
 doubleprecision,parameter :: alpha30 =  30.0/180.0*M_PI
 doubleprecision,parameter :: vecs_armchair(2,2) =  (/  (/ 1.0D0,0.0D0 /) , (/ sin(alpha30) , cos(alpha30) /) /)
 doubleprecision,parameter :: atoms_armchair(2,2) =  (/  (/ 0.0D0,0.0D0 /) , (/ 0.0D0 , 1.0D0/sqrt(3.0) /) /)

 doubleprecision,parameter :: vecs_zigzag(2,2) =   (/  (/ (3.0/2.0)/sqrt(3.0D0)  ,0.5D0 /) , (/ -(3.0/2.0)/sqrt(3.0D0) , 0.5D0 /) /)
 doubleprecision,parameter :: atoms_zigzag(2,2) =  (/  (/ 0.0D0,0.0D0 /) , (/ 1.0D0/sqrt(3.0)  ,0.0D0  /) /)

 doubleprecision,parameter :: pos_offset(2) =  (/ -5.0D0,0.0D0 /)
 doubleprecision,parameter :: pos_offset_zigzag(2) =  (/ -5.0D0,-10.0D0 /)
 doubleprecision :: atom_pos(3),lead_translation_vec(3)
 integer ,parameter :: atom_A = 1 , atom_B = 2 , atom_C = 3
 integer :: atom
 integer :: gindex(30,31)

 call qt%init_system()

 ! --------------------------------------------------------------------------
 ! ARMCHAIR test
 ! --------------------------------------------------------------------------
 k = 0
 atom_pos = 0
 gindex   = 0
 do i = 1 , size(gindex,1)
 do j = 1 , size(gindex,2)
    do atom = atom_A , atom_B
    ! set atom position in space
    atom_pos(1:2) = atoms_armchair(:,atom) + (i-1) * vecs_armchair(:,1) +  (j-1) * vecs_armchair(:,2) + pos_offset
    ! cut some atoms to have rectangular flake
!    if(atom_pos(1) > 0.2 .and. atom_pos(1) < 11.7 .and.  atom_pos(2) < 10.5 ) then
        call qt%qatom%init( (/ atom_pos(1) , atom_pos(2) , 0.0D0 /),flag=atom)
        call qt%qsystem%add_atom(qt%qatom)
        k = k + 1
        gindex(i,j) = k
!    endif
    enddo ! end of atom loop
enddo
enddo

qt%qnnbparam%distance   = 0.6
qt%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_DISTANCE
call qt%qsystem%make_lattice(connect,qt%qnnbparam)

! remove single bonds
do atom = 1 ,qt%qsystem%no_atoms
    if(qt%qsystem%atoms(atom)%no_bonds == 1)   qt%qsystem%atoms(atom)%bActive  = .false.
enddo
call qt%qsystem%make_lattice(connect,qt%qnnbparam)


call qt%qsystem%save_lattice("lattice.dat")

! adding lead
!
quad_positions(1,:) = (/-5.4 ,  0.6  /)
quad_positions(2,:) = (/-4.8 ,  0.2  /)
quad_positions(3,:) = (/ 10.4 , 26.0  /)
quad_positions(4,:) = (/ 9.63 , 26.44  /)

call rect_shape%init_convex_quad(quad_positions)

!!call rect_shape%init_rect(SHAPE_RECTANGLE_XY,0.4D0,1.1D0,0.0D0,11.0D0)
call qt%add_lead(rect_shape,(/1.0D0,0.0D0,0.0D0/))
call qt%leads(1)%print_lead("lead.dat",qt%qsystem%atoms)
call qt%leads(1)%bands("bands.dat",-3.14D0,3.14D0,0.1D0,-15.0D0,15.0D0)



! call qt%qsystem%calc_eigenproblem(0.0D0,0.2D0,150)
! zero_array = 0
! if(qt%qsystem%no_eigenvalues > 0) then
! do i = 1 , size(gindex,1)
! do j = 1 , size(gindex,2)
!    if(qt%qsystem%atoms(gindex(i,j))%bActive) then
!    k = qt%qsystem%atoms(gindex(i,j))%globalIDs(1)
!    write(222,"(500e20.6)"),qt%qsystem%atoms(gindex(i,j))%atom_pos(:),abs(qt%qsystem%eigenvecs(k,:))**2
!    else
!    write(222,"(500e20.6)"),qt%qsystem%atoms(gindex(i,j))%atom_pos(:),zero_array(1:qt%qsystem%no_eigenvalues)
!    endif
! enddo
!    write(222,*),""
! enddo
! endif

call qt%destroy_system()
contains

logical function connect(atomA,atomB,s1,s2,coupling_val)
    use modatom
    implicit none
    type(qatom) :: atomA,atomB
    integer    :: s1,s2
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

end program transporter
