module modshape

ENUM,BIND(C)
    ENUMERATOR :: SHAPE_NONE          = 0
    ENUMERATOR :: SHAPE_RECTANGLE_XY  = 1
END ENUM


type qshape
    integer         :: SHAPE_TYPE
    doubleprecision :: xmin,xmax,ymin,ymax
    contains

    procedure,pass(this) :: init
    procedure,pass(this) :: is_inside
    procedure,pass(this) :: init_rect
    procedure,pass(this) :: is_inside_rect
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

    selectcase(this%SHAPE_TYPE)
    case(SHAPE_RECTANGLE_XY)
            is_inside = this%is_inside_rect(vec)
    case default
            print*,"SYS::ERROR::There is no such type of shape:",this%SHAPE_TYPE
            stop -1
    endselect
    is_inside =  .false.

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
end subroutine init_rect

logical function is_inside_rect(this,vec)
    class(qshape) :: this
    doubleprecision :: vec(3)

    is_inside_rect = .false.
    if( this%xmin < vec(1) .and. this%xmax > vec(1) .and. &
        this%ymin < vec(2) .and. this%ymax > vec(2)  ) is_inside_rect = .true.

end function is_inside_rect

endmodule modshape
