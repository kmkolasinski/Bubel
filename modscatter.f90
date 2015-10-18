! ---------------------------------------------------------------------------------------
!  Created by Krzysztof Kolasinski - 2015
! ---------------------------------------------------------------------------------------
! Module for solving scattering problem in Quantum Mechanics problems.
! Solution for the scattergin wave function is founded with the Wave
! Function Matching approach well described in this paper:
! Scattering Matrices by Wave Function Matching: G. Brocks, V. M. Karpan et. al.
! ---------------------------------------------------------------------------------------
module modscatter
use modsys
use modlead
use modatom
use modshape
use modutils
implicit none

private
integer ,parameter :: QSCATTER_MAX_LEADS = 50
type qscatter
    type(qatom)       :: qatom
    type(nnb_params)  :: qnnbparam  ! auxiliary variable, can be used by user
    type(qsys)   :: qsystem
    type(qlead) , dimension(QSCATTER_MAX_LEADS) :: leads
    integer      :: no_leads

    complex*16,dimension(:)  , allocatable       :: MATHVALS
    integer   ,dimension(:,:), allocatable       :: ROWCOLID
    doubleprecision,allocatable                  :: densities(:,:),qauxvec(:)
    doubleprecision,allocatable                  :: Tn(:)

    integer   :: NO_NON_ZERO_VALUES,NO_VARIABLES

    contains
    procedure,pass(this),public :: init_system
    procedure,pass(this),public :: destroy_system
    procedure,pass(this),public :: add_lead!(this,lshape)
    procedure,pass(this),public :: calculate_modes!(this,Ef)
    procedure,pass(this),public :: construct_matrix!(this,Ef)
    procedure,pass(this),public :: solve!
    procedure,pass(this),public :: save_system!filename
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

    if(allocated(this%MATHVALS))    deallocate(this%MATHVALS)
    if(allocated(this%ROWCOLID))    deallocate(this%ROWCOLID)
    if(allocated(this%densities))   deallocate(this%densities)
    if(allocated(this%qauxvec))     deallocate(this%qauxvec)
    if(allocated(this%Tn))          deallocate(this%Tn)

end subroutine destroy_system

subroutine add_lead(this,lshape,lvec)
    class(qscatter) :: this
    type(qshape) :: lshape
    doubleprecision :: lvec(3)

    this%no_leads = this%no_leads + 1
    call this%leads(this%no_leads)%init_lead(lshape,lvec,this%qsystem%atoms)
end subroutine add_lead

subroutine calculate_modes(this,Ef)
    class(qscatter) :: this
    integer :: i
    doubleprecision :: Ef
    print*,"SYS::SCATTER::Calculating lead modes"
    do i = 1, this%no_leads
        print*,"-----------------------------------------"
        print*,"SYS::LEAD:",i
        print*,"-----------------------------------------"
        call this%leads(i)%calculate_modes(Ef)
    enddo
end subroutine calculate_modes

subroutine save_system(this,filename)
    class(qscatter) :: this
    character(*) :: filename
    integer,parameter :: funit = 9873721
    integer :: i
    open(unit=funit,file=filename)
    write(funit,*),"<system>"
    call this%qsystem%save_lattice(filename,ofunit=funit)
    do i = 1, this%no_leads
        call this%leads(i)%save_lead(filename,funit)
    enddo
    write(funit,*),"</system>"
    close(funit)

end subroutine save_system

subroutine construct_matrix(this,Ef)
    class(qscatter) :: this
    integer   :: i ,j, itmp ,l ,lr,lc,lid , s , no_leads , no_atoms , ts,ta,fs,fa

    doubleprecision :: Ef
    integer,allocatable :: leadFlags(:) , leadIds(:)


    no_leads = this%no_leads
    no_atoms = this%qsystem%no_atoms

    allocate(leadFlags(no_atoms))
    allocate(leadIds(this%qsystem%system_size))
    leadFlags = 0
    leadIds = 0

    ! fill leads array with flags
    do l = 1 , no_leads
        do i = 1 , this%leads(l)%no_sites

            ta = this%leads(l)%l2g(i,1) ! get atom
            ts = this%leads(l)%l2g(i,2) ! get spin
            leadFlags(ta) = l
            leadIds(this%qsystem%atoms(ta)%globalIDs(ts)) = i
        enddo
    enddo



    ! Calculate the number of non-zero elements in Hamiltonian matrix
    itmp = 0
    do i = 1, this%qsystem%no_atoms
        if(this%qsystem%atoms(i)%bActive) then
            itmp = itmp + this%qsystem%atoms(i)%no_bonds
        endif
    enddo
    this%NO_NON_ZERO_VALUES = itmp
    ! The number of unknows is taken from the last global index
    this%NO_VARIABLES       = this%qsystem%system_size

    do l = 1 , no_leads
        this%NO_NON_ZERO_VALUES = this%NO_NON_ZERO_VALUES + this%leads(l)%no_sites**2 * 2
    enddo


    print*,"SYS::SCATTERING: N=",this%NO_VARIABLES," NON ZERO:",this%NO_NON_ZERO_VALUES


    if(allocated(this%MATHVALS)) deallocate(this%MATHVALS)
    if(allocated(this%ROWCOLID)) deallocate(this%ROWCOLID)

    allocate(this%MATHVALS(this%NO_NON_ZERO_VALUES))
    allocate(this%ROWCOLID(this%NO_NON_ZERO_VALUES,2))

    ! Filling matrix and row-col array
    itmp = 0
    do i = 1 ,  this%qsystem%no_atoms

        if(this%qsystem%atoms(i)%bActive) then

        if(leadFlags(i) /= 0) then
            lid = leadFlags(i)
            ! iterate over spins
            do s = 1 , this%qsystem%atoms(i)%no_in_states
                ! get local site in lead
                lr = leadIds(this%qsystem%atoms(i)%globalIDs(s))
                ! iterate over sites in lead
                fa = this%leads(lid)%l2g(lr,1)
                fs = this%leads(lid)%l2g(lr,2)
                do lc = 1 , this%leads(lid)%no_sites

                    ta = this%leads(lid)%l2g(lc,1)
                    ts = this%leads(lid)%l2g(lc,2)
                    itmp = itmp + 1
                    if(lr == lc) then
                    this%MATHVALS(itmp)   = Ef - this%leads(lid)%SigmaMat(lr,lc)
!                    print*,"asd"
                    else
                    this%MATHVALS(itmp)   = -this%leads(lid)%SigmaMat(lr,lc)
                    endif

                    this%ROWCOLID(itmp,1) = this%qsystem%atoms(fa)%globalIDs(fs)
                    this%ROWCOLID(itmp,2) = this%qsystem%atoms(ta)%globalIDs(ts)

                    if(abs(this%MATHVALS(itmp)) < 1.0D-20) then
                    itmp = itmp - 1
                    endif

                enddo
                ! coupling matrix
                do lc = 1 , this%leads(lid)%no_sites
!                    if(lid == 1)then
                    ta = this%leads(lid)%next_l2g(lc,1)
                    ts = this%leads(lid)%next_l2g(lc,2)
                    itmp = itmp + 1
                    this%MATHVALS(itmp)   = -this%leads(lid)%valsTau(lr,lc)
                    this%ROWCOLID(itmp,1) = this%qsystem%atoms(fa)%globalIDs(fs)
                    this%ROWCOLID(itmp,2) = this%qsystem%atoms(ta)%globalIDs(ts)
!                    endif

                    if(abs(this%MATHVALS(itmp)) < 1.0D-20) then
                    itmp = itmp - 1
                    endif
                enddo

            enddo
        else

        do j = 1, this%qsystem%atoms(i)%no_bonds

            itmp = itmp + 1

            fs = this%qsystem%atoms(i)%bonds(j)%fromInnerID
            this%ROWCOLID(itmp,1) = this%qsystem%atoms(i)%globalIDs(fs)

            ta = this%qsystem%atoms(i)%bonds(j)%toAtomID
            ts = this%qsystem%atoms(i)%bonds(j)%toInnerID
            this%ROWCOLID(itmp,2) = this%qsystem%atoms(ta)%globalIDs(ts)
            if(this%ROWCOLID(itmp,1) == this%ROWCOLID(itmp,2)) then
                this%MATHVALS(itmp)   = Ef - this%qsystem%atoms(i)%bonds(j)%bondValue
            else
                this%MATHVALS(itmp)   = -this%qsystem%atoms(i)%bonds(j)%bondValue
            endif

            if(abs(this%MATHVALS(itmp)) < 1.0D-20) then
            itmp = itmp - 1
            endif
        enddo
        endif ! flag leads
        endif ! end of if active atom
    enddo
    this%NO_NON_ZERO_VALUES = itmp

    deallocate(leadFlags)
    deallocate(leadIds)
end subroutine construct_matrix


! --------------------------------------------------------------------
! Solve scattering problem for lead number leadID and given Ef
! Make sure that you run this command before calculate_modes(Ef)
! --------------------------------------------------------------------
subroutine solve(this,leadID,Ef)
    class(qscatter) :: this
    integer,intent(in)  :: leadID
    doubleprecision :: Ef
    doubleprecision :: timer_factorization,total_Tn
    integer ,allocatable   :: HBROWS(:)
    complex*16,allocatable :: phi(:)
    integer :: modin , i , j , lg , la ,ls

    ! --------------------------------------------------------------------
    !
    ! --------------------------------------------------------------------

    call this%construct_matrix(Ef)

    if(allocated(this%densities)) deallocate(this%densities)
    if(allocated(this%qauxvec))   deallocate(this%qauxvec)
    if(allocated(this%Tn))        deallocate(this%Tn)

    allocate(this%densities(this%leads(leadID)%no_in_modes,this%NO_VARIABLES))
    allocate(this%qauxvec  (this%NO_VARIABLES))
    allocate(this%Tn(this%leads(leadID)%no_in_modes))

    allocate(phi(this%NO_VARIABLES))
    allocate(HBROWS(this%NO_VARIABLES+1))

    this%densities = 0
    this%Tn        = 0

    call convert_to_HB(this%NO_NON_ZERO_VALUES,this%ROWCOLID,this%MATHVALS,HBROWS)

    timer_factorization  = get_clock()
    call solve_SSOLEQ(this%NO_VARIABLES, &
                      this%NO_NON_ZERO_VALUES,&
                      this%ROWCOLID(:,2),HBROWS,&
                      this%MATHVALS(:),phi,1)

    timer_factorization = get_clock() - timer_factorization
    print*,"SYS::SCATTERING::Factorization time:",timer_factorization

    this%Tn = 0

    ! ---------------------------------------------------------
    ! Loop over all modes
    ! ---------------------------------------------------------
    do modin = 1 , this%leads(leadID)%no_in_modes

    phi   = 0
    do i = 1 , this%leads(leadID)%no_sites
        la = this%leads(leadID)%l2g(i,1)
        ls = this%leads(leadID)%l2g(i,2)
        lg = this%qsystem%atoms(la)%globalIDs(ls);
        do j = 1 , this%leads(leadID)%no_sites
            phi(lg) = phi(lg) + this%leads(leadID)%LambdaMat(i,j)*this%leads(leadID)%modes(1,modin,j)
        enddo
    enddo

    call solve_SSOLEQ(this%NO_VARIABLES, &
                      this%NO_NON_ZERO_VALUES,&
                      this%ROWCOLID(:,2),HBROWS,&
                      this%MATHVALS(:),phi,2)

    this%densities(modin,:) = abs(phi(:))**2

    total_Tn = 0
    do i = 1 , this%no_leads
        call this%leads(i)%calculate_Tnm(this%qsystem%atoms,phi)
        ! Skip tranmission from input lead
        if( i /= leadID) then
            total_Tn = total_Tn + this%leads(i)%Tnm(1,1)
        endif
    enddo
    this%Tn(modin) = total_Tn
    enddo ! end of loop over modes

    print*,"-----------------------------------------"
    print*,"SYS::SCATTERING problem solved:"
    print*,"        T=",sum(this%Tn)
    print*,"        R=", this%leads(leadID)%no_in_modes-sum(this%Tn)
    print*,"-----------------------------------------"
    ! Free memory
    call solve_SSOLEQ(this%NO_VARIABLES, &
                      this%NO_NON_ZERO_VALUES,&
                      this%ROWCOLID(:,2),HBROWS,&
                      this%MATHVALS(:),phi,3)

    deallocate(phi)
    deallocate(HBROWS)

end subroutine solve




endmodule modscatter
