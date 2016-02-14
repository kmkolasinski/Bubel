! ---------------------------------------------------------------------------------------
!  Created by Krzysztof Kolasinski - 2015
! ---------------------------------------------------------------------------------------
! Module for solving scattering problem in Quantum Mechanics problems.
! Solution for the scattergin wave function is founded with the Wave
! Function Matching approach well described in this paper:
! Scattering Matrices by Wave Function Matching: G. Brocks, V. M. Karpan et. al.
! There is also a posibility to use QTBM to evaluate tranmission probability
! However it may be less stable.
! ---------------------------------------------------------------------------------------
module modscatter
use modcommons
use modsys
use modlead
use modshape
use modalgs

implicit none

private
integer ,parameter :: QSCATTER_MAX_LEADS = 50 ! hard coded maximum number of leads in system
type qscatter
    type(qatom)       :: qatom
    type(nnb_params)  :: qnnbparam  ! auxiliary variable, can be used by user
    type(qsys)        :: qsystem
    type(qlead) , dimension(QSCATTER_MAX_LEADS) :: leads
    integer           :: no_leads

    complex*16,dimension(:)  , allocatable       :: MATHVALS
    integer   ,dimension(:,:), allocatable       :: ROWCOLID
    doubleprecision,allocatable                  :: densities(:,:),qauxvec(:)
    complex*16     ,allocatable                  :: wavefunc(:,:)
    doubleprecision,allocatable                  :: Tn(:),Rn(:)

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
    if(allocated(this%wavefunc))    deallocate(this%wavefunc)
    if(allocated(this%qauxvec))     deallocate(this%qauxvec)
    if(allocated(this%Tn))          deallocate(this%Tn)
    if(allocated(this%Rn))          deallocate(this%Rn)

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
    if(QSYS_DEBUG_LEVEL > 0) then
    print*,"SYS::SCATTER::Calculating lead modes"
    endif
    do i = 1, this%no_leads
        if(QSYS_DEBUG_LEVEL > 0) then
        print*,"-----------------------------------------"
        print*,"SYS::LEAD:",i
        print*,"-----------------------------------------"
        endif
        call this%leads(i)%update_lead(this%qsystem%atoms)
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
    integer   :: s1,s2,ns1,ns2
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
            do j = 1 , this%qsystem%atoms(i)%no_bonds
                itmp = itmp + size(this%qsystem%atoms(i)%bonds(j)%bondMatrix)
            enddo
        endif
    enddo
    this%NO_NON_ZERO_VALUES = itmp
    ! The number of unknows is taken from the last global index
    this%NO_VARIABLES       = this%qsystem%system_size

    do l = 1 , no_leads
        this%NO_NON_ZERO_VALUES = this%NO_NON_ZERO_VALUES + this%leads(l)%no_sites**2 * 2
    enddo

    if(QSYS_DEBUG_LEVEL > 0) then
    print*,"SYS::SCATTERING: N=",this%NO_VARIABLES," NON ZERO:",this%NO_NON_ZERO_VALUES
    endif

    if(allocated(this%MATHVALS)) deallocate(this%MATHVALS)
    if(allocated(this%ROWCOLID)) deallocate(this%ROWCOLID)

    allocate(this%MATHVALS(this%NO_NON_ZERO_VALUES))
    allocate(this%ROWCOLID(this%NO_NON_ZERO_VALUES,2))
    this%ROWCOLID = 0
    this%MATHVALS = 0
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

                    this%MATHVALS(itmp)   = Ef*this%leads(lid)%valsS0(lr,lc)-this%leads(lid)%SigmaMat(lr,lc)
                    this%ROWCOLID(itmp,1) = this%qsystem%atoms(fa)%globalIDs(fs)
                    this%ROWCOLID(itmp,2) = this%qsystem%atoms(ta)%globalIDs(ts)

                    if(abs(this%MATHVALS(itmp)) < 1.0D-20) then
                    itmp = itmp - 1
                    endif
                enddo

                ! coupling matrix
                do lc = 1 , this%leads(lid)%no_sites

                    ta = this%leads(lid)%next_l2g(lc,1)
                    ts = this%leads(lid)%next_l2g(lc,2)
                    itmp = itmp + 1
                    ! ------------------------------------------------------
                    if(QSYS_SCATTERING_METHOD == QSYS_SCATTERING_WFM) then
                        this%MATHVALS(itmp)   = &
                                            Ef*this%leads(lid)%valsS1(lr,lc) -this%leads(lid)%valsTau(lr,lc)
                    ! ------------------------------------------------------
                    else if(QSYS_SCATTERING_METHOD == QSYS_SCATTERING_QTBM) then
                        this%MATHVALS(itmp)   = &
                                            Ef*this%leads(lid)%valsS1(lr,lc)       -this%leads(lid)%valsTau(lr,lc) &
                                          + conjg(Ef*this%leads(lid)%valsS1(lc,lr) -this%leads(lid)%valsTau(lc,lr))
                    endif
                    ! ------------------------------------------------------
                    this%ROWCOLID(itmp,1) = this%qsystem%atoms(fa)%globalIDs(fs)
                    this%ROWCOLID(itmp,2) = this%qsystem%atoms(ta)%globalIDs(ts)

                    if(abs(this%MATHVALS(itmp)) < 1.0D-20) then
                    itmp = itmp - 1
                    endif
                enddo ! end of lc

            enddo
        else ! end of else is lead

        ns1   = this%qsystem%atoms(i)%no_in_states
        do s1 = 1 , ns1
        do j = 1, this%qsystem%atoms(i)%no_bonds
            ns2 = size(this%qsystem%atoms(i)%bonds(j)%bondMatrix,2)
            do s2 = 1 , ns2
            itmp  = itmp + 1
            this%ROWCOLID(itmp,1) = this%qsystem%atoms(i)%globalIDs(s1)
            ta = this%qsystem%atoms(i)%bonds(j)%toAtomID
            this%ROWCOLID(itmp,2) = this%qsystem%atoms(ta)%globalIDs(s2)
            this%MATHVALS(itmp)   = this%qsystem%atoms(i)%bonds(j)%overlapMatrix(s1,s2)*Ef &
                                  - this%qsystem%atoms(i)%bonds(j)%bondMatrix(s1,s2)

            ! Remove zero elements
            if(abs(this%MATHVALS(itmp)) < 1.0D-20) then
            itmp = itmp - 1
            endif

            enddo ! end of s2
        enddo ! end of bonds in i
        enddo ! end of s1
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
    doubleprecision :: timer_factorization,total_Tn,total_Rn
    integer ,allocatable   :: HBROWS(:)
    complex*16,allocatable :: phi(:)
    complex*16             :: cval
    integer :: modin , i , j , lg , la ,ls , nrhs

    ! --------------------------------------------------------------------
    !
    ! --------------------------------------------------------------------

    call this%construct_matrix(Ef)

    if(allocated(this%densities)) deallocate(this%densities)
    if(allocated(this%wavefunc))  deallocate(this%wavefunc)
    if(allocated(this%qauxvec))   deallocate(this%qauxvec)
    if(allocated(this%Tn))        deallocate(this%Tn)
    if(allocated(this%Rn))        deallocate(this%Rn)

    allocate(this%densities(this%leads(leadID)%no_in_modes,this%NO_VARIABLES))
    allocate(this%wavefunc (this%leads(leadID)%no_in_modes,this%NO_VARIABLES))
    allocate(this%qauxvec  (this%NO_VARIABLES))
    allocate(this%Tn(this%leads(leadID)%no_in_modes))
    allocate(this%Rn(this%leads(leadID)%no_in_modes))


    this%densities = 0
    this%wavefunc  = 0
    this%Tn        = 0
    this%Rn        = 0


    ! Skip calculation for case when there is no modes in the input lead
    if(this%leads(leadID)%no_out_modes == 0 ) return;
    nrhs = this%leads(leadID)%no_in_modes
    allocate(phi(this%NO_VARIABLES*nrhs))
    allocate(HBROWS(this%NO_VARIABLES+1))


    call convert_to_HB(this%NO_NON_ZERO_VALUES,this%ROWCOLID,this%MATHVALS,HBROWS)


    timer_factorization  = get_clock()
    call zalg_PARDISO(this%NO_VARIABLES, &
                      this%NO_NON_ZERO_VALUES,&
                      this%ROWCOLID(:,2),HBROWS,&
                      this%MATHVALS(:),nrhs,phi,QSYS_LINSYS_STEP_FACTORIZE,&
                      QSYS_LINSYS_PARDISO_CMPLX_NON_SYM)

    timer_factorization = get_clock() - timer_factorization
    if(QSYS_DEBUG_LEVEL > 0) then
        print*,"SYS::SCATTERING::Factorization time:",timer_factorization
    endif

    this%Tn = 0
    this%Rn = 0
    do i = 1 , this%no_leads
        this%leads(i)%totalT = 0
    enddo

    ! ---------------------------------------------------------
    ! Loop over all modes
    ! ---------------------------------------------------------
    phi = 0
    do modin = 1 , this%leads(leadID)%no_in_modes


    do i = 1 , this%leads(leadID)%no_sites
        la = this%leads(leadID)%l2g(i,1)
        ls = this%leads(leadID)%l2g(i,2)
        lg = this%qsystem%atoms(la)%globalIDs(ls);
        do j = 1 , this%leads(leadID)%no_sites
            ! Choosing between available method
            if(QSYS_SCATTERING_METHOD == QSYS_SCATTERING_QTBM) then
            cval    = this%leads(leadID)%lambdas(M_IN,modin) - this%leads(leadID)%lambdas(M_IN,modin)**(-1)
            phi(lg+(modin-1)*this%NO_VARIABLES) = phi(lg+(modin-1)*this%NO_VARIABLES) &
                                + (this%leads(leadID)%LambdaMat(i,j) - &
                                cval*(conjg(this%leads(leadID)%valsTau(j,i)) - Ef*conjg(this%leads(leadID)%valsS1(j,i)) ) ) &
                                *this%leads(leadID)%modes(M_IN,modin,j)

            else if(QSYS_SCATTERING_METHOD == QSYS_SCATTERING_WFM) then
            cval    = this%leads(leadID)%lambdas(M_IN,modin)**(-1)
            phi(lg+(modin-1)*this%NO_VARIABLES) = phi(lg+(modin-1)*this%NO_VARIABLES) &
                              + (-this%leads(leadID)%LambdaMat(i,j) + &
                                cval*(conjg(this%leads(leadID)%valsTau(j,i)) - Ef*conjg(this%leads(leadID)%valsS1(j,i)) ) ) &
                                *this%leads(leadID)%modes(M_IN,modin,j)
            endif
        enddo ! end of j
    enddo ! end of i
    enddo ! end of modin

    call zalg_PARDISO(this%NO_VARIABLES, &
                      this%NO_NON_ZERO_VALUES,&
                      this%ROWCOLID(:,2),HBROWS,&
                      this%MATHVALS(:),nrhs,phi,QSYS_LINSYS_STEP_SOLVE,&
                      QSYS_LINSYS_PARDISO_CMPLX_NON_SYM)


    do modin = 1 , this%leads(leadID)%no_in_modes
        la = (modin-1)*this%NO_VARIABLES + 1
        ls = (modin)*this%NO_VARIABLES
        this%densities(modin,:) = abs(phi(la:ls))**2
        this%wavefunc (modin,:) = phi(la:ls)

        total_Tn = 0
        total_Rn = 0
        do i = 1 , this%no_leads
            ! Skip tranmission from input lead
            if( i /= leadID) then
                call this%leads(i)%calculate_Tnm(this%qsystem%atoms,this%wavefunc (modin,:))
                total_Tn = total_Tn  + this%leads(i)%modeT
            else
                call this%leads(i)%calculate_Tnm(this%qsystem%atoms,this%wavefunc (modin,:),modin)
                total_Rn = total_Rn  + this%leads(i)%modeT
            endif
            this%leads(i)%totalT = this%leads(i)%totalT + this%leads(i)%modeT
        enddo

        this%Tn(modin) = total_Tn
        this%Rn(modin) = total_Rn

    enddo ! end of loop over modes
    if(QSYS_DEBUG_LEVEL > 0) then
    print*,"-----------------------------------------"
    print*,"SYS::SCATTERING problem solved:"
    print*,"        T=",sum(this%Tn)
    print*,"        R=",sum(this%Rn)
    print*,"      T+R=",sum(this%Tn)+sum(this%Rn)," M=",this%leads(leadID)%no_in_modes
    print*,"-----------------------------------------"
    endif
    ! Free memory
    call zalg_PARDISO(this%NO_VARIABLES, &
                      this%NO_NON_ZERO_VALUES,&
                      this%ROWCOLID(:,2),HBROWS,&
                      this%MATHVALS(:),nrhs,phi,QSYS_LINSYS_STEP_FREE_MEMORY,&
                      QSYS_LINSYS_PARDISO_CMPLX_NON_SYM)

    deallocate(phi)

    deallocate(HBROWS)

!    allocate(phi(this%NO_VARIABLES))
!    allocate(HBROWS(this%NO_VARIABLES+1))
!
!
!    call convert_to_HB(this%NO_NON_ZERO_VALUES,this%ROWCOLID,this%MATHVALS,HBROWS)
!
!    timer_factorization  = get_clock()
!    call solve_SSOLEQ(this%NO_VARIABLES, &
!                      this%NO_NON_ZERO_VALUES,&
!                      this%ROWCOLID(:,2),HBROWS,&
!                      this%MATHVALS(:),phi,QSYS_LINSYS_STEP_FACTORIZE,&
!                      QSYS_LINSYS_PARDISO_CMPLX_NON_SYM)
!
!    timer_factorization = get_clock() - timer_factorization
!    if(QSYS_DEBUG_LEVEL > 0) then
!    print*,"SYS::SCATTERING::Factorization time:",timer_factorization
!    endif
!    this%Tn = 0
!    this%Rn = 0
!    do i = 1 , this%no_leads
!        this%leads(i)%totalT = 0
!    enddo
!
!    ! ---------------------------------------------------------
!    ! Loop over all modes
!    ! ---------------------------------------------------------
!    do modin = 1 , this%leads(leadID)%no_in_modes
!
!    phi   = 0
!    do i = 1 , this%leads(leadID)%no_sites
!        la = this%leads(leadID)%l2g(i,1)
!        ls = this%leads(leadID)%l2g(i,2)
!        lg = this%qsystem%atoms(la)%globalIDs(ls);
!        do j = 1 , this%leads(leadID)%no_sites
!            ! Choosing between available method
!            if(QSYS_SCATTERING_METHOD == QSYS_SCATTERING_QTBM) then
!            cval    = this%leads(leadID)%lambdas(M_IN,modin) - this%leads(leadID)%lambdas(M_IN,modin)**(-1)
!            phi(lg) = phi(lg) + (this%leads(leadID)%LambdaMat(i,j) - &
!                                cval*(conjg(this%leads(leadID)%valsTau(j,i)) - Ef*conjg(this%leads(leadID)%valsS1(j,i)) ) ) &
!                                *this%leads(leadID)%modes(M_IN,modin,j)
!
!            else if(QSYS_SCATTERING_METHOD == QSYS_SCATTERING_WFM) then
!            cval    = this%leads(leadID)%lambdas(M_IN,modin)**(-1)
!            phi(lg) = phi(lg) + (-this%leads(leadID)%LambdaMat(i,j) + &
!                                cval*(conjg(this%leads(leadID)%valsTau(j,i)) - Ef*conjg(this%leads(leadID)%valsS1(j,i)) ) ) &
!                                *this%leads(leadID)%modes(M_IN,modin,j)
!            endif
!        enddo
!    enddo
!
!    call solve_SSOLEQ(this%NO_VARIABLES, &
!                      this%NO_NON_ZERO_VALUES,&
!                      this%ROWCOLID(:,2),HBROWS,&
!                      this%MATHVALS(:),phi,QSYS_LINSYS_STEP_SOLVE,&
!                      QSYS_LINSYS_PARDISO_CMPLX_NON_SYM)
!
!    this%densities(modin,:) = abs(phi(:))**2
!    this%wavefunc (modin,:) = phi(:)
!
!    total_Tn = 0
!    total_Rn = 0
!    do i = 1 , this%no_leads
!        ! Skip tranmission from input lead
!        if( i /= leadID) then
!            call this%leads(i)%calculate_Tnm(this%qsystem%atoms,phi)
!            total_Tn = total_Tn  + this%leads(i)%modeT
!        else
!            call this%leads(i)%calculate_Tnm(this%qsystem%atoms,phi,modin)
!            total_Rn = total_Rn  + this%leads(i)%modeT
!        endif
!        this%leads(i)%totalT = this%leads(i)%totalT + this%leads(i)%modeT
!    enddo
!
!
!    this%Tn(modin) = total_Tn
!    this%Rn(modin) = total_Rn
!
!    enddo ! end of loop over modes
!    if(QSYS_DEBUG_LEVEL > 0) then
!    print*,"-----------------------------------------"
!    print*,"SYS::SCATTERING problem solved:"
!    print*,"        T=",sum(this%Tn)
!    print*,"        R=",sum(this%Rn)
!    print*,"      T+R=",sum(this%Tn)+sum(this%Rn)," M=",this%leads(leadID)%no_in_modes
!    print*,"-----------------------------------------"
!    endif
!    ! Free memory
!    call solve_SSOLEQ(this%NO_VARIABLES, &
!                      this%NO_NON_ZERO_VALUES,&
!                      this%ROWCOLID(:,2),HBROWS,&
!                      this%MATHVALS(:),phi,QSYS_LINSYS_STEP_FREE_MEMORY,&
!                      QSYS_LINSYS_PARDISO_CMPLX_NON_SYM)
!
!    deallocate(phi)
!    deallocate(HBROWS)

end subroutine solve


endmodule modscatter
