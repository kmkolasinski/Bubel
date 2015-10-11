module modshape

ENUM,BIND(C)
    ENUMERATOR :: SHAPE_NONE          = 0
    ENUMERATOR :: SHAPE_RECTANGLE_XY  = 1
    ENUMERATOR :: SHAPE_CONVEX_QUAD_XY= 2
    ENUMERATOR :: SHAPE_RANGE_3D      = 3
END ENUM


type qshape
    integer         :: SHAPE_TYPE
    doubleprecision :: xmin,xmax,ymin,ymax
    doubleprecision :: quad_points(4,2)
    doubleprecision :: range_direction(3),range_base_pos(3),range_n(3),range_lenght
    contains

    procedure,pass(this) :: init
    procedure,pass(this) :: is_inside

    procedure,pass(this) :: init_rect
    procedure,pass(this) :: init_convex_quad
    procedure,pass(this) :: init_range_3d

    procedure,pass(this) :: is_inside_rect_xy
    procedure,pass(this) :: is_inside_convex_quad_xy
    procedure,pass(this) :: is_inside_range_3d

    procedure,pass(this) :: flush_shape_data_to_file!(this,file_id)
end type qshape

public :: qshape
contains



subroutine init(this,shape_type)
    class(qshape) :: this
    integer  :: shape_type
    this%SHAPE_TYPE = shape_type
end subroutine init


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
    case default
            print*,"SYS::SHAPE::ERROR::There is no such type of shape:",this%SHAPE_TYPE
            stop -1
    endselect

end function is_inside

subroutine init_rect(this,shape_type,xmin,xmax,ymin,ymax)
    class(qshape) :: this
    integer  :: shape_type
    doubleprecision :: xmin,xmax,ymin,ymax
    this%SHAPE_TYPE = shape_type
    this%xmin = xmin
    this%xmax = xmax
    this%ymin = ymin
    this%ymax = ymax
    print"(A,2e12.3,A,2e12.3,A)"," SYS::SHAPE::Initializing rectagle box: x=(",xmin,xmax,"), y=(",ymin,ymax,")"
end subroutine init_rect

subroutine init_convex_quad(this,quad_points)
    class(qshape) :: this
    doubleprecision :: quad_points(4,2)

    this%SHAPE_TYPE  = SHAPE_CONVEX_QUAD_XY
    this%quad_points = quad_points
    print"(A)"," SYS::SHAPE::Initializing quad shape"
end subroutine init_convex_quad

subroutine init_range_3d(this,base,dir)
    class(qshape) :: this
    doubleprecision :: base(3) , dir(3)

    this%SHAPE_TYPE      = SHAPE_RANGE_3D
    this%range_base_pos  = base
    this%range_direction = dir

    this%range_lenght    = sqrt(sum(dir**2))
    this%range_n         = dir/this%range_lenght
    print"(A)"," SYS::SHAPE::Initializing range 3d shape"

end subroutine init_range_3d

logical function is_inside_rect_xy(this,vec)
    class(qshape) :: this
    doubleprecision :: vec(3)
    is_inside_rect_xy = .false.
    if( this%xmin < vec(1) .and. this%xmax > vec(1) .and. &
        this%ymin < vec(2) .and. this%ymax > vec(2)  ) is_inside_rect_xy = .true.

end function is_inside_rect_xy


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


logical function is_inside_range_3d(this,vec)
    class(qshape) :: this
    doubleprecision :: vec(3)
    doubleprecision :: dr(3)
    logical         :: test
    dr   = vec - this%range_base_pos
    test = (sum(this%range_n*dr) > 0) .and. (sum(this%range_n*dr) < this%range_lenght)
    is_inside_range_3d = test

end function is_inside_range_3d

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
    case default
            print*,"SYS::SHAPE::ERROR::There is no such type of shape:",this%SHAPE_TYPE," cannot flush data to file."
            stop -1
    endselect

endsubroutine flush_shape_data_to_file
endmodule modshape
