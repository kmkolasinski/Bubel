! ---------------------------------------------------------------------------------------
!  Created by Krzysztof Kolasinski - 2015
!
!
! ---------------------------------------------------------------------------------------
module modscatter
use modsys
use modlead
use modatom
use modshape
private
integer ,parameter :: QSCATTER_MAX_LEADS = 50
type qscatter
    type(qatom)       :: qatom
    type(nnb_params)  :: qnnbparam  ! auxiliary variable, can be used by user
    type(qsys)   :: qsystem
    type(qlead) , dimension(QSCATTER_MAX_LEADS) :: leads
    integer      :: no_leads
    contains
    procedure,pass(this),public :: init_system
    procedure,pass(this),public :: destroy_system
    procedure,pass(this),public :: add_lead!(this,lshape)
endtype

public :: qscatter

contains

subroutine init_system(this)
    class(qscatter) :: this
    call this%qsystem%init()
    this%no_leads = 0
end subroutine init_system

subroutine destroy_system(this)
    class(qscatter) :: this
    integer :: i
    call this%qsystem%destroy()
    do i = 1, this%no_leads
        call this%leads(i)%destroy()
    enddo
    this%no_leads = 0
end subroutine destroy_system

subroutine add_lead(this,lshape,lvec)
    class(qscatter) :: this
    type(qshape) :: lshape
    doubleprecision :: lvec(3)

    this%no_leads = this%no_leads + 1

    call this%leads(this%no_leads)%init_lead(lshape,lvec,this%qsystem%atoms)


end subroutine add_lead




endmodule modscatter
