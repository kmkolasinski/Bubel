program transporter
 use modscatter
 use modsys
 use modshape
 use modunits
 implicit none
 type(qscatter) :: qt
 type(qshape) :: rect_shape



 integer :: i , j, k , N
 doubleprecision :: dx,zero_array(100)
 doubleprecision,parameter :: lattice_vectors(2,2) =  (/  (/ 1.0D0,0.0D0 /) , (/ sin(30.0/180.0*M_PI) , cos(30.0/180.0*M_PI) /) /)
 doubleprecision,parameter :: atoms_positions(2,2) =  (/  (/ 0.0D0,0.0D0 /) , (/ 0.0D0 , 1.0D0/sqrt(3.0) /) /)
 doubleprecision,parameter :: pos_offset(2) =  (/ -5.0D0,0.0D0 /)
 doubleprecision :: atom_pos(3)
 integer ,parameter :: atom_A = 1 , atom_B = 2 , atom_C = 3
 integer :: atom
 integer :: gindex(19,11)
 print*,"Program start"


 print*,"latice vector 1=",lattice_vectors(:,1)
 print*,"latice vector 2=",lattice_vectors(:,2)
 call qt%init_system()
 dx = 1.0
 k = 0
 atom_pos = 0
 gindex   = 0
 do i = 1 , size(gindex,1)
 do j = 1 , size(gindex,2)
    do atom = atom_A , atom_B
    ! set atom position in space
    atom_pos(1:2) = atoms_positions(:,atom) + (i-1) * lattice_vectors(:,1) +  (j-1) * lattice_vectors(:,2) + pos_offset
    ! cut some atoms to have rectangular flake
    if(atom_pos(1) > 0.2 .and. atom_pos(1) < 11.7 ) then
        call qt%qatom%init( (/ atom_pos(1) , atom_pos(2) , 0.0D0 /),flag=atom)
        call qt%qsystem%add_atom(qt%qatom)
        k = k + 1
        gindex(i,j) = k
    endif
    enddo ! end of atom loop

    !if( abs(i - size(gindex,1)/2) > 5 .and. j > 20 ) qsystem%qatom%bActive = .false.
    !if( abs(i - size(gindex,1)/2.0)**2 + abs(j - size(gindex,2)/2.0)**2 > 200 ) qsystem%qatom%bActive = .false.
!    print*,qt%qatom%atom_pos


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

call rect_shape%init_rect(SHAPE_RECTANGLE_XY,0.4D0,1.1D0,0.0D0,9.0D0)
call qt%add_lead(rect_shape)

call qt%leads(1)%bands("rel.txt",-3.14D0,3.14D0,0.1D0,-20.0D0,20.0D0)


! do i = 1 , 50
!
!    !print*,i,qsystem%atoms(i)%no_bonds
!    do k = 1 , qsystem%atoms(i)%no_bonds
!        !print"(A,4i,2f6.2)"," ",k,qsystem%atoms(i)%bonds(k)%fromInnerID,qsystem%atoms(i)%bonds(k)%toAtomID,qsystem%atoms(i)%bonds(k)%toInnerID,qsystem%atoms(i)%bonds(k)%bondValue
!    enddo
! enddo

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


! do i = 1 , qsystem%no_atoms
!    k = qsystem%atoms(i)%globalIDs(1)
!    write(222,"(500e20.6)"),qsystem%atoms(i)%atom_pos(:),abs(qsystem%eigenvecs(k,:))**2
! enddo

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
    test = atomA%flag == atom_A .and. atomB%flag == atom_B
    test = test .or. atomB%flag == atom_A .and. atomA%flag == atom_B

    if(test) then
        connect = .true.
        coupling_val = 4.0D0
    endif



end function connect

end program transporter
