! ---------------------------------------------------------------------------------------
!  Created by Krzysztof Kolasinski - 2015
!
!
! ---------------------------------------------------------------------------------------
module modscatter
use modshape
use modsys
use modlead
use modatom
private
integer ,parameter :: QSCATTER_MAX_LEADS = 50
type qscatter

    type(qsys)   :: qsystem
    type(qlead) , dimension(QSCATTER_MAX_LEADS) :: leads
    integer      :: no_leads
    contains
endtype

public :: qscatter

contains



endmodule modscatter
