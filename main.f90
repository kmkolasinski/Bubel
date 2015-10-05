! TODO:
! 1. dodac mapowanie z lokalny na globalny
! 2. zrobic aby dzialo liczenie stanow wlasnych
! 3. zrobic rysowanie siatki

program transporter
 use modscatter
 use modshape
 implicit none
 type(qscatter) :: qt
 type(qshape) :: rect_shape



 integer :: i , j, k , N
 doubleprecision :: dx,zero_array(100)
 integer :: gindex(10,5)
 print*,"Program start"


 call qt%init_system()
 dx = 0.1
  k = 0
 do i = 1 , size(gindex,1)
 do j = 1 , size(gindex,2)
    call qt%qatom%init((/ (i-1) * dx , (j-1) * dx + (i-1) * dx * 0.0 , 0.0 * dx /))
    !if( abs(i - size(gindex,1)/2) > 5 .and. j > 20 ) qsystem%qatom%bActive = .false.
    !if( abs(i - size(gindex,1)/2.0)**2 + abs(j - size(gindex,2)/2.0)**2 > 200 ) qsystem%qatom%bActive = .false.
    call qt%qsystem%add_atom(qt%qatom)
    k = k + 1
    gindex(i,j) = k
 enddo
 enddo

 do i = 1 , size(gindex,1)
 do j = 1 , size(gindex,2)

 enddo
 enddo

! qsystem%atoms(gindex(15,5))%bActive = .false.

! print*,"System size:",qsystem%no_atoms

 qt%qnnbparam%max_distance = (/2*dx,2*dx,0.0D0/)
 call qt%qsystem%make_lattice(connect,qt%qnnbparam)


 call rect_shape%init_rect(SHAPE_RECTANGLE_XY,-0.5*dx,0.5*dx,-0.5*dx,size(gindex,2)*dx)


 call qt%qsystem%save_lattice("lattice.dat")
 call qt%add_lead(rect_shape)
! do i = 1 , 50
!
!    !print*,i,qsystem%atoms(i)%no_bonds
!    do k = 1 , qsystem%atoms(i)%no_bonds
!        !print"(A,4i,2f6.2)"," ",k,qsystem%atoms(i)%bonds(k)%fromInnerID,qsystem%atoms(i)%bonds(k)%toAtomID,qsystem%atoms(i)%bonds(k)%toInnerID,qsystem%atoms(i)%bonds(k)%bondValue
!    enddo
! enddo

 call qt%qsystem%calc_eigenproblem(0.0D0,0.2D0,50)
 zero_array = 0
 if(qt%qsystem%no_eigenvalues > 0) then
 do i = 1 , size(gindex,1)
 do j = 1 , size(gindex,2)
    if(qt%qsystem%atoms(gindex(i,j))%bActive) then
    k = qt%qsystem%atoms(gindex(i,j))%globalIDs(1)
    write(222,"(500e20.6)"),qt%qsystem%atoms(gindex(i,j))%atom_pos(:),abs(qt%qsystem%eigenvecs(k,:))**2
    else
    write(222,"(500e20.6)"),qt%qsystem%atoms(gindex(i,j))%atom_pos(:),zero_array(1:qt%qsystem%no_eigenvalues)
    endif
 enddo
    write(222,*),""
 enddo
 endif


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
    integer :: xdiff,ydiff
    doubleprecision :: dydiff
    xdiff = NINT((atomA%atom_pos(1)-atomB%atom_pos(1))/dx)
    ydiff = NINT((atomA%atom_pos(2)-atomB%atom_pos(2))/dx)
    dydiff = (atomA%atom_pos(2)-atomB%atom_pos(2))/dx
    connect = .false.

    if(s1 == s2) then
        if( xdiff == 0 .and. ydiff == 0 ) then
            connect = .true.
            coupling_val = 4.0D0
            if( abs(atomA%atom_pos(1) - 4.0) < 0.05 ) coupling_val = 4.5D0
            if( abs(atomA%atom_pos(1) - 6.0) < 0.05 ) coupling_val = 4.5D0
        else if( xdiff ==  1 .and. ydiff == 0 ) then
            connect = .true.
            coupling_val = -1.0D0
        else if( xdiff == -1 .and. ydiff == 0 ) then
            connect = .true.
            coupling_val = -1.0D0
        else if( xdiff ==  0 .and. ydiff == 1 ) then
            connect = .true.
            coupling_val = -1.0D0 !* EXP(+cmplx(0,1)*dydiff)
        else if( xdiff ==  0 .and. ydiff ==-1 ) then
            connect = .true.
            coupling_val = -1.0D0 !* EXP(+cmplx(0,1)*dydiff)
!        else if( xdiff ==  1 .and. ydiff == 1 ) then
!            connect = .true.
!            coupling_val = -0.01D0
!        else if( xdiff == -1 .and. ydiff ==-1 ) then
!            connect = .true.
!            coupling_val = -0.01D0
        endif
    endif

end function connect

end program transporter
