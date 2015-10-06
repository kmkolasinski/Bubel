! ---------------------------------------------------------------------------------------
!  Created by Krzysztof Kolasinski - 2015
!
!
! ---------------------------------------------------------------------------------------
module modlead
use modshape
use modsys
use modatom
implicit none
private
complex*16,parameter :: II = cmplx(0.0D0,1.0D0)
type qlead
    type(qshape) :: lead_shape
    complex*16,dimension(:,:) , allocatable :: valsH0 ! Hamiltonian of the lead
!    integer,   dimension(:,:) , allocatable :: rowcolH0

    complex*16,dimension(:,:) , allocatable :: valsTau ! Hamiltonian of the lead
    integer :: no_sites
  !  integer,   dimension(:,:) , allocatable :: rowcolTau
    integer,dimension(:,:),allocatable        :: l2g,next_l2g ! local id to (atom,spin)
    contains
    procedure,pass(this) :: init_lead!()
    procedure,pass(this) :: bands!(this,filename,kmin,kmax,dk,Emin,Emax)
    procedure,pass(this) :: destroy!()
endtype qlead

public :: qlead
contains

subroutine init_lead(this,lshape,all_atoms)
    class(qlead) :: this
    type(qshape) :: lshape
    type(qatom),dimension(:) :: all_atoms
    integer ,allocatable ,dimension(:) :: tmp_g2l
    integer :: i,j,b
    integer :: no_sites,system_size, atom_id,bond_id,bond_atom_id
    integer :: irow,icol,aid,gid,sid,lid
    complex :: cval
    doubleprecision,dimension(3) :: cellA_pos,cellB_pos,tmp_pos,cellBA_vec

    this%lead_shape = lshape

    system_size = 0
    do i = 1 , size(all_atoms)
        if(all_atoms(i)%bActive) then
            system_size = system_size + all_atoms(i)%no_in_states
        endif ! end of active atom
    enddo

    allocate(tmp_g2l(system_size))

    no_sites = 0
    tmp_g2l = 0
!    print*,"g2l"
    do i = 1 , size(all_atoms)
        if(all_atoms(i)%bActive) then
            if( lshape%is_inside(all_atoms(i)%atom_pos) == .true. ) then
                do j = 1 , all_atoms(i)%no_in_states
                    no_sites = no_sites + 1
                    tmp_g2l(all_atoms(i)%globalIDs(j)) = no_sites
                    cellA_pos = all_atoms(i)%atom_pos
!                    print*,all_atoms(i)%globalIDs(j),no_sites
                enddo

            endif
        endif
    enddo
!    print*,"H0 size:",no_sites
    this%no_sites = no_sites




    if(allocated(this%valsH0))deallocate(this%valsH0)
    !if(allocated(this%rowcolH0))deallocate(this%rowcolH0)

    if(allocated(this%valsTau))deallocate(this%valsTau)
    if(allocated(this%l2g))deallocate(this%l2g)
    if(allocated(this%next_l2g))deallocate(this%next_l2g)
!    if(allocated(this%l2a))deallocate(this%l2a)
    !if(allocated(this%rowcolTau))deallocate(this%rowcolTau)

    allocate(this%valsH0 (no_sites,no_sites))
    allocate(this%valsTau(no_sites,no_sites))
    allocate(this%l2g(no_sites,2))
    allocate(this%next_l2g(no_sites,2))
!    allocate(this%l2a(no_sites))
    this%valsH0   = 0
    this%valsTau  = 0
    this%l2g = 0
    this%next_l2g = 0

    no_sites = 0


!    print*,"l2g"
    do i = 1 , size(all_atoms)
        if(all_atoms(i)%bActive) then
            if( lshape%is_inside(all_atoms(i)%atom_pos) == .true. ) then
                do j = 1 , all_atoms(i)%no_in_states
                    no_sites = no_sites + 1
                    this%l2g(no_sites,1) = i
                    this%l2g(no_sites,2) = j
!                    print*,no_sites,"a=",i,"s=",j,"g=",all_atoms(i)%globalIDs(j)
                enddo
                tmp_pos = all_atoms(i)%atom_pos

                if( tmp_pos(1) < cellA_pos(1) .or. &
                    tmp_pos(2) < cellA_pos(2) .or. &
                    tmp_pos(3) < cellA_pos(3)) then
                    cellA_pos = tmp_pos
                endif
            endif
        endif
    enddo

    ! initializing lowest atom position in next cell
    do i = 1 , no_sites
        atom_id = this%l2g(i,1) ! get atom index
!        print*,"First atom:",i," id=",atom_id," no_bonds=",all_atoms(atom_id)%no_bonds
        do b = 1 , all_atoms(atom_id)%no_bonds
!            print*,"    Testint bond=",b

            bond_atom_id = all_atoms(atom_id)%bonds(b)%toAtomID
            bond_id      = all_atoms(bond_atom_id)%globalIDs(1)
!            print*,"    Connection with atom:",bond_atom_id, " local index:",tmp_g2l(bond_id)
            if( tmp_g2l(bond_id) == 0 ) then
                cellB_pos = all_atoms(bond_atom_id)%atom_pos
!                print*,"Lowest atom in next cell:",cellB_pos

            endif
        enddo
    enddo

    ! initializing lowest atom position in next cell
    do i = 1 , no_sites
        atom_id = this%l2g(i,1) ! get atom index
        do b = 1 , all_atoms(atom_id)%no_bonds
            bond_atom_id = all_atoms(atom_id)%bonds(b)%toAtomID
            bond_id      = all_atoms(bond_atom_id)%globalIDs(1)
            if( tmp_g2l(bond_id) == 0 ) then
                tmp_pos = all_atoms(bond_atom_id)%atom_pos

                if( tmp_pos(1) < cellB_pos(1) .or. &
                    tmp_pos(2) < cellB_pos(2) .or. &
                    tmp_pos(3) < cellB_pos(3)) then
                    cellB_pos = tmp_pos
                endif
            endif
        enddo
    enddo

!    print*,"Lowest atom in zero cell:",cellA_pos
!    print*,"Lowest atom in next cell:",cellB_pos
    cellBA_vec = cellB_pos - cellA_pos

    do i = 1 , no_sites
        atom_id = this%l2g(i,1) ! get atom index
        irow = i
!        print*,"site:",i," -> ",atom_id
        do b = 1 , all_atoms(atom_id)%no_bonds
            if( all_atoms(atom_id)%bonds(b)%fromInnerID == this%l2g(i,2) ) then

            bond_atom_id = all_atoms(atom_id)%bonds(b)%toAtomID
            bond_id      = all_atoms(bond_atom_id)%globalIDs(1)
            if( tmp_g2l(bond_id) == 0 ) then
!                print*,"Sasiednia komurka"
                ! we search for an atom in the main uint cell with the same position
                do j = 1 , no_sites
                    cellA_pos = all_atoms(this%l2g(j,1))%atom_pos + cellBA_vec
                    cellB_pos = all_atoms(bond_atom_id)%atom_pos
                    if( sqrt(sum( cellB_pos - cellA_pos )**2) < 1.0E-20 ) then
                        exit
                    endif
                enddo
!                print*,"atom:",bond_atom_id," to ten sam co=",this%l2g(j,1)

                aid = this%l2g(j,1)
                sid = all_atoms(atom_id)%bonds(b)%toInnerID

                gid  = all_atoms(aid)%globalIDs(sid)
                lid  = tmp_g2l(gid)
                irow = i
                icol = lid
                cval = all_atoms(atom_id)%bonds(b)%bondValue
!                print*,irow,icol," h=",cval
                this%valsTau(irow,icol) = cval

                this%next_l2g(i,1) = bond_atom_id
                this%next_l2g(i,2) = sid
            else
!                print*,"at l=",atom_id," checking bound ",b," with atom:",bond_atom_id
                aid = bond_atom_id
                sid = all_atoms(atom_id)%bonds(b)%toInnerID
                gid = all_atoms(aid)%globalIDs(sid)
                lid = tmp_g2l(gid)
                irow = i
                icol = lid
                cval = all_atoms(atom_id)%bonds(b)%bondValue
!                print*,irow,icol,b," h=",cval
                this%valsH0(irow,icol) = cval
            endif
            endif
        enddo

    enddo


    deallocate(tmp_g2l)

end subroutine init_lead


subroutine destroy(this)
    class(qlead) :: this
    if(allocated(this%valsH0))deallocate(this%valsH0)
!    if(allocated(this%rowcolH0))deallocate(this%rowcolH0)
    if(allocated(this%valsTau))deallocate(this%valsTau)
    if(allocated(this%l2g))deallocate(this%l2g)
    if(allocated(this%next_l2g))deallocate(this%next_l2g)
!    if(allocated(this%l2a))deallocate(this%l2a)
!    if(allocated(this%rowcolTau))deallocate(this%rowcolTau)
end subroutine destroy


subroutine bands(this,filename,kmin,kmax,dk,Emin,Emax)
    class(qlead) :: this
    double precision :: kmin,kmax,dk,Emin,Emax
    character(*)     :: filename


    complex*16 :: kvec
    doubleprecision :: skank
    integer :: i,j,N


    ! -------------------------------------------------
    !                    LAPACK
    ! -------------------------------------------------
    !     .. Local Scalars ..
    INTEGER          INFO, LWORK, LRWORK, LIWORK, IL, IU, M , LWMAX
    DOUBLE PRECISION ABSTOL, VL, VU

    !     .. Local Arrays ..
    INTEGER,allocatable,dimension(:)          :: ISUPPZ, IWORK
    DOUBLE PRECISION,allocatable,dimension(:) :: EVALS( : ), RWORK( : )
    COMPLEX*16,allocatable,dimension(:)       :: Z( :, : ), WORK( : )
    COMPLEX*16,allocatable,dimension(:,:)     :: Hamiltonian



    N = this%no_sites
    LWMAX = N*50

    allocate(ISUPPZ( 2*N ))
    allocate(EVALS ( N ))
    allocate(Z( N, N ))

    allocate(IWORK ( LWMAX ))
    allocate(RWORK ( LWMAX ))
    allocate(WORK  ( LWMAX ))
    allocate(Hamiltonian(N,N))

    !
    !     Query the optimal workspace.
    !
    LWORK  = -1
    LRWORK = -1
    LIWORK = -1
    !
    ABSTOL = -1.0
    VL     = Emin
    VU     = Emax




    do i = 1 , N
    do j = 1 , N
            Hamiltonian(i,j) = this%valsH0(i,j) + this%valsTau(i,j) + conjg(this%valsTau(j,i))
    enddo
    enddo

    open(unit = 33345, file= "Hamiltonian.txt" )
    do i = 1 , N
    write(33345,"(5000f10.4)"),dble(Hamiltonian(i,:))
    enddo
    close(33345)


    CALL ZHEEVR( "N", 'Values', 'Lower', N, Hamiltonian, N, VL, VU, IL,&
                &             IU, ABSTOL, M, EVALS, Z, N, ISUPPZ, WORK, LWORK, RWORK,&
                &             LRWORK, IWORK, LIWORK, INFO )

    if( INFO /= 0 ) then
        print*,"  SPINMODPOP: Problem zle zdefiniowany INFO:",INFO
        stop
    endif
    LWORK  = MIN( LWMAX, INT( WORK( 1 ) ) )
    LRWORK = MIN( LWMAX, INT( RWORK( 1 ) ) )
    LIWORK = MIN( LWMAX, IWORK( 1 ) )
    deallocate( WORK)
    deallocate(RWORK)
    deallocate(IWORK)

    allocate(IWORK ( LIWORK ))
    allocate(RWORK ( LRWORK ))
    allocate(WORK  ( LWORK  ))


!    print*,"LWORK=" ,LWORK
!    print*,"LRWORK=",LRWORK
!    print*,"LIWORK=",LIWORK
!    print*,"N=",N


    open(unit = 782321, file= filename )
    do skank = kmin , kmax , dk

        kvec = II*skank
        do i = 1 , N
        do j = 1 , N
                Hamiltonian(i,j) = this%valsH0(i,j) + this%valsTau(i,j)*exp(kvec) + conjg(this%valsTau(j,i)*exp(kvec))
        enddo
        enddo
        ABSTOL = -1.0
        !     Set VL, VU to compute eigenvalues in half-open (VL,VU] interval
        VL = Emin
        VU = Emax
        !
        !     Solve eigenproblem.
        !
        EVALS = 0
        M  = 0
        CALL ZHEEVR( "V", 'Values', 'Lower', N, Hamiltonian, N, VL, VU, IL,&
        &             IU, ABSTOL, M, EVALS, Z, N, ISUPPZ, WORK, LWORK, RWORK,&
        &             LRWORK, IWORK, LIWORK, INFO )

        if(M > 0) then
            write(782321,"(1000e20.10)"),skank,EVALS(1:M)
        endif

    enddo
    close(782321)

    deallocate(ISUPPZ)
    deallocate(EVALS)
    deallocate(Z)

    deallocate(IWORK)
    deallocate(RWORK)
    deallocate(WORK)
    deallocate(Hamiltonian)

end subroutine bands


endmodule modlead
