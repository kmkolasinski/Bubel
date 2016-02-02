MODULE modutils
use modalgs
use modcommons
implicit none
private
integer,parameter :: FD_SIZE = QTOOLS_FD_EXPANSION_MAX_ORDER
doubleprecision,dimension(0:FD_SIZE,-FD_SIZE:FD_SIZE) :: v_dnFdxn,v_dnFdyn,v_dnFdzn

public :: qtools_fd_template
public :: v_dnFdxn,v_dnFdyn,v_dnFdzn
public :: dnFdxn,dnFdyn,dnFdzn

contains


subroutine read_filed (UnitNum, FileName, NumRows, NumCols, Array )
  integer, intent (in) :: UnitNum
  character (len=*), intent (in) :: FileName
  integer, intent (in) :: NumRows, NumCols
  double precision, dimension (1:NumRows, 1:NumCols), intent (out) :: Array
  integer :: i, j
  open (unit=UnitNum, file=FileName, status='old', action='read' )
  do i=1, NumRows
     read (UnitNum, *) (Array (i, j), j=1,NumCols)
  end do
  close (UnitNum)
  return
end subroutine read_filed

subroutine read_filei (UnitNum, FileName, NumRows, NumCols, Array )
  integer, intent (in) :: UnitNum
  character (len=*), intent (in) :: FileName
  integer, intent (in) :: NumRows, NumCols
  integer, dimension (1:NumRows, 1:NumCols), intent (out) :: Array
  integer :: i, j
  open (unit=UnitNum, file=FileName, status='old', action='read' )
  do i=1, NumRows
     read (UnitNum, *) (Array (i, j), j=1,NumCols)
  end do
  close (UnitNum)
  return
end subroutine read_filei

integer function qtools_factorial(n) result(rval)
    integer :: n,i
    rval = 1
    do i = 1 , n
        rval = rval * i
    enddo
end function qtools_factorial

integer function qtools_Dij(i,j) result(rval)
    integer :: i,j
    rval = 0
    if(i==j) rval = 1
end function qtools_Dij
! ------------------------------------------------------------------- !
!               FINITE DIFFERENCE EXPANSION GENERATOR
! ------------------------------------------------------------------- !
! Calculate convolution vector for approximated finite difference
! derivative of order (deriv) and error order (order). You can use
! provided internal arrays: dnFdxn,dnFdyn,dnFdzn
! d_array(deriv,convolution) is filled for given deriv and values of
!       convolution are stored in range -m:m where
!       m = floor( (deriv + order - 1)/2 )
! ------------------------------------------------------------------- !
subroutine qtools_fd_template(d_array,deriv,order)
    doubleprecision,dimension(0:FD_SIZE,-FD_SIZE:FD_SIZE) :: d_array
    integer :: order,deriv
    integer :: f,i,n,m,p,d
    doubleprecision, allocatable,dimension(:,:) :: Cni


    ! clean values for given derivative

    d = deriv
    p = order


    m = floor( (p + d - 1.0)/2 )
    if(2*m + 1 /= p + d) then
        p = p + 1 ! fix order of derivative
    endif
    m = floor( (p + d - 1.0)/2 )

    if(deriv > FD_SIZE) then
        print*,"QTOOLS::Required derivative order (d) is to big"
        print*,"        Provided value:",deriv," Max. value:",FD_SIZE
        stop
    endif
    if(    m > FD_SIZE) then
        print*,"QTOOLS::Required derivative error order (p) is to big"
        print*,"        Provided value:",order," Max. value:",FD_SIZE
        stop
    endif



    ! Create coeficient matrix
    allocate(Cni(0:p + d - 1,-m:m))
    Cni = 0
    do n =  0 , p + d - 1
    do i = -m , m
        Cni(n,i) = i**n
    enddo
    enddo
    ! invert and get the values for p-th derivative
    call dalg_invmat(Cni)

    Cni = transpose(Cni)
    f = qtools_factorial(d)
    d_array(deriv,:) = 0
    do i = -m , m
         d_array(deriv,i) = f*Cni(deriv,i)
    enddo
    if(QSYS_DEBUG_LEVEL > 0) then
        print"(A,i3,A,i3)"," QTOOLS::Created FD template for d=",d," derivative of order p=",p
        print*,"        Calculated template (convolution vector):"
        print*,"        | ======================== |"
        print*,"        |  F(x +-i e ) |   VALUE   |"
        print*,"        | ======================== |"
        do i = -m , -1
             print"(A,i2,A,f10.6,A)","         | F( x -",-i,"*e ) |",d_array(deriv,i)," |"
        enddo
        print"(A,f10.6,A)","         |    F( x )    |",d_array(deriv,i)," |"
        do i =  1 , m
             print"(A,i2,A,f10.6,A)","         | F( x +",i,"*e ) |",d_array(deriv,i)," |"
        enddo
        print*,"        | ======================== |"
    endif
    deallocate(Cni)
end subroutine qtools_fd_template

doubleprecision function dnFdxn(n,dx,dy,dz) result(rval)
    integer :: n,dx,dy,dz
    rval = v_dnFdxn(n,dx)*qtools_Dij(0,dy)*qtools_Dij(0,dz)
end function dnFdxn

doubleprecision function dnFdyn(n,dx,dy,dz) result(rval)
    integer :: n,dx,dy,dz
    rval = v_dnFdyn(n,dy)*qtools_Dij(0,dx)*qtools_Dij(0,dz)
end function dnFdyn

doubleprecision function dnFdzn(n,dx,dy,dz) result(rval)
    integer :: n,dx,dy,dz
    rval = v_dnFdzn(n,dz)*qtools_Dij(0,dx)*qtools_Dij(0,dy)
end function dnFdzn

END MODULE modutils

