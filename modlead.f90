! ---------------------------------------------------------------------------------------
!  Created by Krzysztof Kolasinski - 2015
!
!
! ---------------------------------------------------------------------------------------
module modlead
use modshape
use modsys
private

type qlead
    type(qshape) :: lead_shape
    contains
    procedure,pass(this) :: init_lead!()
endtype qlead

public :: qlead
contains

subroutine init_lead(this,lshape)
    class(qlead) :: this
    type(qshape) :: lshape

    doubleprecision :: vec(3)
    this%lead_shape = lshape

    !print*,"vals=", this%lead_shape%is_inside(vec)

end subroutine init_lead

endmodule modlead
