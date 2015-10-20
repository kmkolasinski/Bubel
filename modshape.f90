! ---------------------------------------------------------------------------------------
!  Created by Krzysztof Kolasinski - 2015
!  Area/Volume definitions for leads.
!
! ---------------------------------------------------------------------------------------
module modshape
use modcommons
implicit none
! ---------------------------------------------------------------------------------------
! Current supported lead areas:
! ---------------------------------------------------------------------------------------
ENUM,BIND(C)
    ENUMERATOR :: SHAPE_NONE          = 0
    ENUMERATOR :: SHAPE_RECTANGLE_XY  = 1
    ENUMERATOR :: SHAPE_CONVEX_QUAD_XY= 2
    ENUMERATOR :: SHAPE_RANGE_3D      = 3
    ENUMERATOR :: SHAPE_ATOMS_LIST    = 4
END ENUM

type qshape
    integer         :: SHAPE_TYPE
    doubleprecision :: xmin,xmax,ymin,ymax
    doubleprecision :: quad_points(4,2)
    doubleprecision :: range_direction(3),range_base_pos(3),range_n(3),range_lenght


    integer           :: no_atoms
    type(qatom),dimension(:),allocatable :: atoms_list


    contains
    procedure,pass(this) :: init
    procedure,pass(this) :: destroy_shape
    procedure,pass(this) :: copyFrom

    procedure,pass(this) :: is_inside

    procedure,pass(this) :: init_rect
    procedure,pass(this) :: init_convex_quad
    procedure,pass(this) :: init_range_3d
    procedure,pass(this) :: init_atoms_list

    procedure,pass(this) :: is_inside_rect_xy
    procedure,pass(this) :: is_inside_convex_quad_xy
    procedure,pass(this) :: is_inside_range_3d
    procedure,pass(this) :: is_inside_atom_list

    procedure,pass(this) :: flush_shape_data_to_file!(this,file_id)
end type qshape

public :: qshape
contains

! ----------------------------------------------------------
! Default 'empty' initialization procedure:
! ----------------------------------------------------------
subroutine init(this,shape_type)
    class(qshape) :: this
    integer  :: shape_type
    this%SHAPE_TYPE = shape_type
    this%no_atoms   = 0
end subroutine init


subroutine destroy_shape(this)
    class(qshape) :: this
    if(allocated(this%atoms_list)) deallocate(this%atoms_list)
    this%SHAPE_TYPE = SHAPE_NONE
    this%no_atoms   = 0
end subroutine destroy_shape

subroutine copyFrom(this,source)
    class(qshape) :: this
    type(qshape)  :: source

    this%SHAPE_TYPE = source%SHAPE_TYPE
    this%xmin = source%xmin
    this%xmax = source%xmax
    this%ymin = source%ymin
    this%ymax = source%ymax
    this%quad_points        = source%quad_points
    this%range_direction    = source%range_direction
    this%range_base_pos     = source%range_base_pos
    this%range_n            = source%range_n
    this%range_lenght       = source%range_lenght
    this%no_atoms           = source%no_atoms

    if(allocated(this%atoms_list)) deallocate(this%atoms_list)
    if(source%no_atoms > 0) then
        if(allocated(this%atoms_list)) deallocate(this%atoms_list)
        allocate(this%atoms_list(source%no_atoms))
        this%atoms_list(1:source%no_atoms) = source%atoms_list(1:source%no_atoms)
    endif

end subroutine copyFrom
! ---------------------------------------------------------------
! Check is 3D point lies in the area defined by shape.
! If YES returns true if not false
! vec(3) - double precision array containing 3D position in space
! ---------------------------------------------------------------
logical function is_inside(this,vec)
    class(qshape) :: this
    doubleprecision :: vec(3)
    is_inside =  .false.
    selectcase(this%SHAPE_TYPE)
    case(SHAPE_RECTANGLE_XY)
            is_inside = this%is_inside_rect_xy(vec)
    case(SHAPE_CONVEX_QUAD_XY)
            is_inside = this%is_inside_convex_quad_xy(vec)
    case(SHAPE_RANGE_3D)
            is_inside = this%is_inside_range_3d(vec)
    case(SHAPE_ATOMS_LIST)
            is_inside = this%is_inside_atom_list(vec)

    case default
            print*,"SYS::SHAPE::ERROR::There is no such type of shape:",this%SHAPE_TYPE
            stop -1
    endselect
end function is_inside

! ---------------------------------------------------------------
! Initialize shape as a simple 2D rectangle
! shape_type - currently it can be only SHAPE_RECTANGLE_XY
! xmin,xmax,ymin,ymax - position of left lower and right upper corner
! ---------------------------------------------------------------
subroutine init_rect(this,shape_type,xmin,xmax,ymin,ymax)
    class(qshape) :: this
    integer  :: shape_type
    doubleprecision :: xmin,xmax,ymin,ymax
    this%SHAPE_TYPE = shape_type
    this%xmin = xmin
    this%xmax = xmax
    this%ymin = ymin
    this%ymax = ymax
    this%no_atoms   = 0
    print"(A,2e12.3,A,2e12.3,A)"," SYS::SHAPE::Initializing rectagle box: x=(",xmin,xmax,"), y=(",ymin,ymax,")"
end subroutine init_rect

! ---------------------------------------------------------------
! Initialize shape as a convex quad
! quad_points(4,2) - array of corners coordinates quad_points(1,:) = (x1,y1)
!                    pay attention to keep good order of points (anticlockwise)
! ---------------------------------------------------------------
subroutine init_convex_quad(this,quad_points)
    class(qshape) :: this
    doubleprecision :: quad_points(4,2)

    this%SHAPE_TYPE  = SHAPE_CONVEX_QUAD_XY
    this%quad_points = quad_points
    this%no_atoms    = 0
    print"(A)"," SYS::SHAPE::Initializing quad shape"
end subroutine init_convex_quad


! ---------------------------------------------------------------
! Initialize shape as a 3D range: is a volume between two infinite planes
! position of firt plane is defined by base vector and normal "dir".
! next plane is shifted by vector dir.
! ---------------------------------------------------------------
subroutine init_range_3d(this,base,dir)
    class(qshape) :: this
    doubleprecision :: base(3) , dir(3)

    this%SHAPE_TYPE      = SHAPE_RANGE_3D
    this%range_base_pos  = base
    this%range_direction = dir

    this%range_lenght    = sqrt(sum(dir**2))
    this%range_n         = dir/this%range_lenght
    this%no_atoms        = 0
    print"(A)"," SYS::SHAPE::Initializing range 3d shape"

end subroutine init_range_3d


subroutine init_atoms_list(this,atoms)
    class(qshape)             :: this
    type(qatom), dimension(:) :: atoms

    this%SHAPE_TYPE = SHAPE_ATOMS_LIST
    this%no_atoms   = size(atoms)

    if(allocated(this%atoms_list)) deallocate(this%atoms_list)
    allocate(this%atoms_list(size(atoms)))

    this%atoms_list(:) = atoms(:)


    print"(A,i)"," SYS::SHAPE::Initializing atoms list of size:",this%no_atoms
end subroutine init_atoms_list
! ---------------------------------------------------------------
! Check if point vec lies in the SHAPE_RECTANGLE_XY shape
! ---------------------------------------------------------------
logical function is_inside_rect_xy(this,vec)
    class(qshape) :: this
    doubleprecision :: vec(3)
    is_inside_rect_xy = .false.
    if( this%xmin < vec(1) .and. this%xmax > vec(1) .and. &
        this%ymin < vec(2) .and. this%ymax > vec(2)  ) is_inside_rect_xy = .true.
end function is_inside_rect_xy

! ---------------------------------------------------------------
! Check if point vec lies in the SHAPE_CONVEX_QUAD_XY shape
! ---------------------------------------------------------------
logical function is_inside_convex_quad_xy(this,vec)
    class(qshape) :: this
    doubleprecision :: vec(3)
    integer :: i,j
    logical :: test
    integer,parameter :: x = 1 , y = 2
    test = .false.
    j = 4
    do i = 1 , 4
        if( (this%quad_points(i,y) > vec(y)) /= (this%quad_points(j,y) > vec(y)) .and.&
            (vec(x) < (this%quad_points(j,x) - this%quad_points(i,x))* &
            ( vec(y) - this%quad_points(i,y) )/ (this%quad_points(j,y)-this%quad_points(i,y)) + this%quad_points(i,x))) then
            test = .not. test
            endif
        j = i
    enddo
    is_inside_convex_quad_xy = test
end function is_inside_convex_quad_xy

! ---------------------------------------------------------------
! Check if point vec lies in the SHAPE_RANGE_3D shape
! ---------------------------------------------------------------
logical function is_inside_range_3d(this,vec)
    class(qshape)   :: this
    doubleprecision :: vec(3)
    ! locals
    doubleprecision :: dr(3)
    logical         :: test

    dr   = vec - this%range_base_pos
    test = (sum(this%range_n*dr) > 0) .and. (sum(this%range_n*dr) < this%range_lenght)
    is_inside_range_3d = test
end function is_inside_range_3d

! ---------------------------------------------------------------
! Check if point vec belongs to som atoms inside shape
! ---------------------------------------------------------------
logical function is_inside_atom_list(this,vec)
    class(qshape) :: this
    doubleprecision :: vec(3)
    doubleprecision :: TOLERANCE = 1.0D-5 , dist
    integer :: i
    is_inside_atom_list = .false.

    do i = 1, this%no_atoms
        if(.not.this%atoms_list(i)%bActive)cycle
        dist = sqrt(sum( ( vec - this%atoms_list(i)%atom_pos )**2 ))
        if(dist < TOLERANCE) then
            is_inside_atom_list = .true.
            exit
        endif
    enddo
end function is_inside_atom_list

! ---------------------------------------------------------------
! Flush shape data to file:
! ---------------------------------------------------------------
subroutine flush_shape_data_to_file(this,file_id)
    class(qshape) :: this
    integer :: file_id

    selectcase(this%SHAPE_TYPE)
    case(SHAPE_RECTANGLE_XY)
            write(file_id,"(A)"),"<shape_type>SHAPE_RECTANGLE_XY</shape_type>"
            write(file_id,"(A)"),"<shape_data>"
            write(file_id,"(4e20.6)"),this%xmin,this%ymin,this%xmax,this%ymax
            write(file_id,"(A)"),"</shape_data>"
    case(SHAPE_CONVEX_QUAD_XY)
            write(file_id,"(A)"),"<shape_type>SHAPE_CONVEX_QUAD_XY</shape_type>"
            write(file_id,"(A)"),"<shape_data>"
            write(file_id,"(8e20.6)"),transpose(this%quad_points)
            write(file_id,"(A)"),"</shape_data>"

    case(SHAPE_RANGE_3D)
            write(file_id,"(A)"),"<shape_type>SHAPE_RANGE_3D</shape_type>"
            write(file_id,"(A)"),"<shape_data>"
            write(file_id,"(6e20.6)"),this%range_base_pos,this%range_direction
            write(file_id,"(A)"),"</shape_data>"
    case(SHAPE_ATOMS_LIST)
            write(file_id,"(A)"),"<shape_type>SHAPE_ATOMS_LIST</shape_type>"

    case default
            print*,"SYS::SHAPE::ERROR::There is no such type of shape:",this%SHAPE_TYPE," cannot flush data to file."
            stop -1
    endselect

endsubroutine flush_shape_data_to_file
endmodule modshape
