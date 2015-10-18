! ---------------------------------------------------------------------------------------
!  Created by Krzysztof Kolasinski - 2015
!
!
! ---------------------------------------------------------------------------------------

! -----------------------------------------------
! Stucture which holds the unit cell of periotic
! structure. In calculated hamiltonian of that
! cell and coupling matrix Tau to next unit cell.
! Those information can be used to generate the
! band structure and for scattering problems.
! -----------------------------------------------
module modlead
use modshape
use modsys
use modatom
use modutils
implicit none
private
complex*16,parameter :: II = cmplx(0.0D0,1.0D0)
integer ,parameter :: M_IN = 1 , M_OUT = 2
type qlead
    type(qshape) :: lead_shape                         ! contains the area/volume in which are
                                                       ! located atoms which belong to the lead.
    doubleprecision :: lead_vector(3)                  ! direction of propagation of the lead in the periodic structur

    complex*16,dimension(:,:) , allocatable :: valsH0  ! Hamiltonian of the lead - unit cell
    complex*16,dimension(:,:) , allocatable :: valsTau ! Coupling matrix to the next unit cell

    integer :: no_sites                                ! number of unknowns (number of sites in the lead,
                                                       ! including spin states)

    integer :: no_in_modes                             ! when modes are calculated this number contains no. of incoming modes
    integer :: no_out_modes                            ! same but for outgoing modes
    integer :: no_in_em                                ! number of "incoming" evanescent modes
    integer :: no_out_em
    ! ----------------------------------------------------------------------
    ! Modes parameters:
    ! According to this paper:http://www.psi-k.org/newsletters/News_80/Highlight_80.pdf
    ! see equation: (53) Arrays below are filled after calling calculate_modes function
    ! ----------------------------------------------------------------------
    complex*16,allocatable :: lambdas(:,:)             ! lambda(M_IN,:) => \lambda_+ and lambda(M_OUT,:) => \lambda_-
    complex*16,allocatable :: modes (:,:,:)            ! modes({M_IN,M_OUT},M,:) N-th mode vector
    complex*16,allocatable :: SigmaMat(:,:)            ! See Eq. (68) here we added H0 to this matrix
    complex*16,allocatable :: LambdaMat(:,:)           ! See Eq. (67) is the definition of the Qm matrix.
    doubleprecision,allocatable :: Tnm(:,:)
    doubleprecision,allocatable :: currents(:,:)       ! Fluxes carried by mode M currents({M_IN,M_OUT},M) - note modes are sorted
                                                       ! from largest current to smallest
    complex*16,allocatable      :: UTildeDagger(:,:)   ! Dual Vectors matrix see Eq. (55)


    integer,dimension(:,:),allocatable        ::      l2g ! mapping from local id in lead to global id (:,1)  = atom_ID , (:,2) = spin_ID
    integer,dimension(:,:),allocatable        :: next_l2g ! the same as above but mapping to the atoms in the next unit cell
    logical :: LEAD_BAD_NEARST = .false.
    contains
    procedure,pass(this) :: init_lead!()
    procedure,pass(this) :: bands!(this,filename,kmin,kmax,dk,Emin,Emax)
    procedure,pass(this) :: destroy!()
    procedure,pass(this) :: save_lead!(this,filename,funit)
    procedure,pass(this) :: calculate_modes!(this,Ef)
    procedure,pass(this) :: calculate_Tnm!(this,all_atoms,n,phi)

endtype qlead

public :: qlead , M_IN , M_OUT
contains


! ------------------------------------------------------------------------
! Initialize lead structure - which is a unit cell in periodic lattice.
! lshape - contains information of the volume in which all the atoms which
!          belong to the lead are located.
! lvec   - contains vector (x,y,z) of the translational symmetry of the lead
! all_atoms - list of all atoms in the system.
! Based on those both information lead will generate Hamiltonian matrix
! for atoms located in this lead as well as coupling matrix to the next
! unit cell.
! Note: Unit cell must lie at the edge of the system such that there will
! be only one unit cell ajacent to main unit cell of the lead.
! ------------------------------------------------------------------------
subroutine init_lead(this,lshape,lvec,all_atoms)
    class(qlead) :: this
    type(qshape) :: lshape
    type(qatom),dimension(:) :: all_atoms
    doubleprecision :: lvec(3)
    ! local variables
    integer ,allocatable ,dimension(:) :: tmp_g2l ! contains mapping between global ID of the state (atom+spin)
                                                  ! to local id of that site (atom,spin) in unit cell
    integer :: i,j,k,b
    integer :: no_sites,system_size, atom_id,bond_atom_id,bond_id
    integer :: irow,icol,aid,gid,sid,lid
    integer :: no_atoms
    complex :: cval
    integer,allocatable :: next_cell_atoms(:)
    doubleprecision,dimension(3) :: cellA_pos,cellB_pos,tmp_pos,cellBA_vec
    doubleprecision :: minimum_distance,dist

    ! copy shape for future possible usage
    this%lead_shape  = lshape
    this%lead_vector = lvec

    print*,"SYS::LEAD::Initializing lead unit cell parameters"
    ! calculate the number of states in whole system: sum over active atoms and their spins
    system_size = 0
    do i = 1 , size(all_atoms)
        if(all_atoms(i)%bActive) then
            system_size = system_size + all_atoms(i)%no_in_states
        endif ! end of active atom
    enddo

    allocate(tmp_g2l(system_size))
    ! -----------------------------------------------------------
    ! Find all the atoms which are located in the lead area/volume
    ! and creat mapping between its globalID and localID
    ! -----------------------------------------------------------
    no_sites = 0
    no_atoms = 0
    tmp_g2l  = 0
    do i = 1 , size(all_atoms)
        if(all_atoms(i)%bActive) then
            if( lshape%is_inside(all_atoms(i)%atom_pos) == .true. ) then
                no_atoms = no_atoms + 1
                do j = 1 , all_atoms(i)%no_in_states ! loop over spins of that atom

                    no_sites = no_sites + 1
                    tmp_g2l(all_atoms(i)%globalIDs(j)) = no_sites ! creating mappin globalID => local id in lead
                    cellA_pos = all_atoms(i)%atom_pos

                enddo
            endif ! end of if this atom is located in the lead
        endif ! end of if active atom
    enddo ! end of do all_atoms

    allocate(next_cell_atoms(no_atoms))
    this%no_sites = no_sites
    print*,"SYS::LEAD::Number of states in lead=",no_sites," with n=",no_atoms," atoms"


    if(allocated(this%valsH0))  deallocate(this%valsH0)
    if(allocated(this%valsTau)) deallocate(this%valsTau)
    if(allocated(this%l2g))     deallocate(this%l2g)
    if(allocated(this%next_l2g))deallocate(this%next_l2g)

    allocate(this%valsH0    (no_sites,no_sites))
    allocate(this%valsTau   (no_sites,no_sites))
    allocate(this%l2g       (no_sites,2))
    allocate(this%next_l2g  (no_sites,2))

    this%valsH0   = 0
    this%valsTau  = 0
    this%l2g      = 0
    this%next_l2g = 0



    ! -----------------------------------------------------------
    ! Now we want to find some unique atom in lead -> has lowest
    ! corrdinates (X_min,Y_min,Z_min). Then we will find the same
    ! atom but in next unit cell. Thus we will be able to calculate
    ! translation vector, between lead and next unit cell.
    ! -----------------------------------------------------------

    no_sites = 0
    no_atoms = 0
    do i = 1 , size(all_atoms) ! Search in all active atoms
        if(all_atoms(i)%bActive) then
            ! if atom is in the unit cell
            tmp_pos = all_atoms(i)%atom_pos-lvec
            if( lshape%is_inside(tmp_pos) == .true. ) then
                no_atoms = no_atoms + 1

                if( no_atoms > this%no_sites ) then
                    print*,"---------------------------------------------------------------------------"
                    print*,"SYS::LEAD::WARNING::Unit cell parameters (shape area?) has not "
                    print*,"           been chosen properly. Plot lead data to see where "
                    print*,"           is the problem."
                    print*,"---------------------------------------------------------------------------"
                    this%LEAD_BAD_NEARST = .true.
                    return
                endif

                next_cell_atoms(no_atoms) = i
!                print*,"next cell i=",i
            endif



            if( lshape%is_inside(all_atoms(i)%atom_pos) == .true. ) then
                ! create mapping between localID (l2g array) and the atom (j) and its spin (j)

                do j = 1 , all_atoms(i)%no_in_states
                    no_sites = no_sites + 1
                    this%l2g(no_sites,1) = i ! } mapping
                    this%l2g(no_sites,2) = j ! }
!                    print*,no_sites,i,j
                enddo

                ! Find atom with the lowest coordinates in that lead:
                tmp_pos = all_atoms(i)%atom_pos
                if( tmp_pos(1) < cellA_pos(1) .or. &
                    tmp_pos(2) < cellA_pos(2) .or. &
                    tmp_pos(3) < cellA_pos(3)) then
                    cellA_pos = tmp_pos
                endif
            endif ! end of if atom is in lead
        endif ! end of if active atom
    enddo ! end of do all_atoms

    ! -----------------------------------------------------
    ! Searching for the minium distance between atoms. Stage 1
    ! -----------------------------------------------------
    do i = 1 , no_sites
        atom_id = this%l2g(i,1) ! get atom index , this%l2g(i,2) tells what is the spinID of that atom
        ! iterate over all bonds in that atom (A)
        do b = 1 , all_atoms(atom_id)%no_bonds
        bond_atom_id = all_atoms(atom_id)%bonds(b)%toAtomID
        if(bond_atom_id /= atom_id)then
            minimum_distance = sqrt(sum((all_atoms(atom_id)%atom_pos-all_atoms(bond_atom_id)%atom_pos)**2))
        endif
        enddo
    enddo

    ! Searching for the minium distance between atoms. Stage 2
    do i = 1 , no_sites
        atom_id = this%l2g(i,1) ! get atom index , this%l2g(i,2) tells what is the spinID of that atom
        ! iterate over all bonds in that atom (A)
        do b = 1 , all_atoms(atom_id)%no_bonds
        bond_atom_id = all_atoms(atom_id)%bonds(b)%toAtomID
        if(bond_atom_id /= atom_id)then
            dist = sqrt(sum((all_atoms(atom_id)%atom_pos-all_atoms(bond_atom_id)%atom_pos)**2))
            if(dist < minimum_distance ) minimum_distance = dist
        endif
        enddo
    enddo
!    print*,"minum distance=",minimum_distance

    no_atoms = 0
    do i = 1 , size(all_atoms) ! Search in all active atoms
        if(all_atoms(i)%bActive) then
            ! if atom is in the unit cell
            tmp_pos = all_atoms(i)%atom_pos - lvec
            if( lshape%is_inside(tmp_pos) == .true. ) then
                ! search for atom with the same position
                do j = 1 , no_sites

                    cellA_pos = all_atoms(this%l2g(j,1))%atom_pos
!                    print"(2i4,5f8.4)",i,j,sqrt(sum( tmp_pos - cellA_pos )**2),tmp_pos
                    if( sqrt(sum( (tmp_pos - cellA_pos)**2 )) < minimum_distance*1.0D-5 ) then
                        exit
                    endif

                enddo
!                print*,"i=",i,"j=",j
                if( j > no_sites ) then
                    print"(A)",              " SYS::LEAD::ERROR::The translation"
                    print"(A,i9,A,3e12.4,A)","                   at atom:",i," with position r=(",all_atoms(i)%atom_pos,")"
                    print"(A,3e12.4,A)"     ,"                   with translation vector: v=(",lvec,")"
                    print"(A)",              "                   is impossible by applying lead translation vector. Maybe some of the "
                    print"(A)",              "                   atoms have not been inlcuded in the lead. Check lead area and T vector."
                    stop -1
                endif
                ! j keeps local ID of the same atom in unit
                do k = 1 , all_atoms(this%l2g(j,1))%no_in_states
!                    print*,"j=",j,"spin k=",k
                    b = tmp_g2l(all_atoms(this%l2g(j,1))%globalIDs(k))
                    this%next_l2g(b,1) = i
                    this%next_l2g(b,2) = this%l2g(b,2)
!                    print*,b,this%next_l2g(b,1),this%next_l2g(b,2)
!                    print*,"atom a=",i," z b=",j," ",this%l2g(j,1)
                enddo

            endif

        endif ! end of if active atom
    enddo ! end of do all_atoms


    cellBA_vec = lvec
    print"(A,3f12.4,A)"," SYS::LEAD::Lead translation vector XYZ=(",cellBA_vec,")"

    ! -----------------------------------------------------------
    ! Creating the Hamiltonian of the Lead and Coupling matrix
    ! -----------------------------------------------------------
    do i = 1 , no_sites
        atom_id = this%l2g(i,1) ! get atom index , this%l2g(i,2) tells what is the spinID of that atom
        ! iterate over all bonds in that atom (A)
        do b = 1 , all_atoms(atom_id)%no_bonds
            ! Filter out bond from atom A but which do not have spin this%l2g(i,1)
            if( all_atoms(atom_id)%bonds(b)%fromInnerID == this%l2g(i,2) ) then
            ! get ids of associated atom
            bond_atom_id = all_atoms(atom_id)%bonds(b)%toAtomID
            bond_id      = all_atoms(bond_atom_id)%globalIDs(1)
            ! if tmp_g2l(bond_id) == 0 it means that this atom with id = bond_atom_id
            ! is not in the lead, so it has to be located in the next unit cell.
            if( tmp_g2l(bond_id) == 0 ) then
                ! now we search for an atom in the main uint cell with the same position
                ! what assiociated atom in the next lead.

                do j = 1 , no_sites
                    cellA_pos = all_atoms(this%l2g(j,1))%atom_pos + cellBA_vec
                    cellB_pos = all_atoms(bond_atom_id)%atom_pos
                    if( sqrt(sum( (cellB_pos - cellA_pos )**2)) < minimum_distance*1.0D-5 ) then
                        exit
                    endif
                enddo
                if( j > no_sites ) then
                    print"(A)",              " SYS::LEAD::ERROR::The translation"
                    print"(A,i9,A,3e12.4,A)","                   from atom:",atom_id," at position r=(",all_atoms(atom_id)%atom_pos,")"
                    print"(A,i9,A,3e12.4,A)","                   to atom  :",bond_atom_id," at position r=(",all_atoms(bond_atom_id)%atom_pos,")"
                    print"(A)",              "                   is impossible by applying lead translation vector. Maybe some of the "
                    print"(A)",              "                   atoms have not been inlcuded in the lead. Check lead area and T vector."
                    cycle
                    !stop -1
                endif
                ! now "j" contains of localID of the state in the lead which is equivalent
                ! to same state but in the next unit cell.

                aid  = this%l2g(j,1) ! get globalID of that atom
                sid  = all_atoms(atom_id)%bonds(b)%toInnerID ! and its spin
                gid  = all_atoms(aid)%globalIDs(sid) ! get global ID of site (atom+spin)
                lid  = tmp_g2l(gid) ! remap it to local ID in the lead
                irow = i
                icol = lid
                cval = all_atoms(atom_id)%bonds(b)%bondValue

                ! fill coupling matrix
                this%valsTau(irow,icol) = cval
!                print*,i,lid,cval


            else ! in case of coupling occures within lead

                aid = bond_atom_id ! id of connected atom
                sid = all_atoms(atom_id)%bonds(b)%toInnerID ! and its spin ID
                gid = all_atoms(aid)%globalIDs(sid) ! convert to global state id (atom+spin)
                lid = tmp_g2l(gid) ! remap to local ID in lead
                irow = i
                icol = lid
                cval = all_atoms(atom_id)%bonds(b)%bondValue
                this%valsH0(irow,icol) = cval
!                print*,
!                print*,irow,icol,atom_id,aid,sqrt(sum((all_atoms(atom_id)%atom_pos-all_atoms(aid)%atom_pos)**2 ))
            endif ! end of else if "atom is in the lead or outside"
            endif ! end of if bond comes from the atom with the same spin what site has
        enddo ! end of do over bond in site

    enddo ! end of do over sites in lead
    ! deallocate unnecessary array
    deallocate(tmp_g2l)
    deallocate(next_cell_atoms)

end subroutine init_lead

! ------------------------------------------------------------------------
! Free allocated memory
! ------------------------------------------------------------------------
subroutine destroy(this)
    class(qlead) :: this
    print*,"SYS::LEAD::freeing memory"
    if(allocated(this%valsH0))  deallocate(this%valsH0)
    if(allocated(this%valsTau)) deallocate(this%valsTau)
    if(allocated(this%l2g))     deallocate(this%l2g)
    if(allocated(this%next_l2g))deallocate(this%next_l2g)

    if(allocated(this%modes))        deallocate(this%modes)
    if(allocated(this%lambdas))      deallocate(this%lambdas)
    if(allocated(this%SigmaMat))     deallocate(this%SigmaMat)
    if(allocated(this%LambdaMat))    deallocate(this%LambdaMat)
    if(allocated(this%UTildeDagger)) deallocate(this%UTildeDagger)
    if(allocated(this%Tnm))          deallocate(this%Tnm)
    if(allocated(this%currents))     deallocate(this%currents)

    this%no_sites    = 0
    this%no_in_modes = 0
    this%no_out_em   = 0
    this%no_out_modes= 0
    this%no_in_em    = 0
end subroutine destroy

! ------------------------------------------------------------------------
! Calculate the band structure for lead. Lead has to be initialized before.
! filename      - the name of the file with generated data (k_vector,Energy1,Energy2,...)
! kmin,kmax,dk  - the range of the perfomed scan with step dk. k is dimensionless.
! Emin,Emax     - search for the eigenvalues which lie in that energy window.
! ------------------------------------------------------------------------
subroutine bands(this,filename,kmin,kmax,dk,Emin,Emax)
    class(qlead) :: this
    double precision :: kmin,kmax,dk,Emin,Emax
    character(*)     :: filename

    ! Local variables
    complex*16      :: kvec
    doubleprecision :: skank
    integer         :: i,j,N


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

    print*,"SYS::LEAD::Generating band structure data to file:",filename
    ! The size of the system
    N = this%no_sites


    ! Initialize lapack and allocate arrays
    LWMAX = N*24

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


    ! Perform scan
    open(unit = 782321, file= filename )
    do skank = kmin , kmax + dk/2 , dk

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


! ------------------------------------------------------------------------
! Calculate the modes for a given Fermi energy. Those modes can be used
! to build scattering boundary condidionts in Quantum transport methods.
! Ef - value of Fermi level energy
! ------------------------------------------------------------------------
subroutine calculate_modes(this,Ef)
    class(qlead)    :: this
    doubleprecision :: Ef
    integer :: N
    complex*16  , allocatable , dimension(:,:)    :: MA,MB,Z
    complex*16, allocatable , dimension(:,:)      :: Mdiag,d,c
    complex*16, allocatable , dimension(:,:,:)    :: blochF

    integer :: k,i,j,no_in,no_out,no_e_in,no_e_out

    INTEGER                                      :: LDVL, LDVR , LWMAX , LWORK , INFO
    COMPLEX*16 , dimension(:) ,     allocatable  :: ALPHA , BETA , WORK
    double precision, dimension(:), allocatable  :: RWORK

    doubleprecision :: tmpc,current,time,dval1,dval2
    COMPLEX*16 :: DUMMY(1,1),lambda


    time = get_clock()
    N = this%no_sites

    print*,"SYS::LEAD::Finding lead modes for N=",N
    ! ---------------------------------------
    ! Memoru allocations
    ! ---------------------------------------
    allocate(MA (2*N,2*N))
    allocate(MB (2*N,2*N))
    allocate(Z  (2*N,2*N))
    allocate(Mdiag(N,N))
    allocate(blochF(2,N,N))
    allocate(d(2*N,N))
    allocate(c(2*N,N))


    if(allocated(this%modes))     deallocate(this%modes)
    if(allocated(this%lambdas))  deallocate(this%lambdas)
    if(allocated(this%SigmaMat)) deallocate(this%SigmaMat)
    if(allocated(this%LambdaMat)) deallocate(this%LambdaMat)
    if(allocated(this%Tnm))      deallocate(this%Tnm)
    if(allocated(this%currents))      deallocate(this%currents)
    if(allocated(this%UTildeDagger)) deallocate(this%UTildeDagger)

    allocate(this%modes(2,N,N))
    allocate(this%lambdas(2,N))
    allocate(this%SigmaMat(N,N))
    allocate(this%LambdaMat(N,N))
    allocate(this%currents(2,N))
    allocate(this%UTildeDagger(N,N))
    this%currents = 0
    this%UTildeDagger = 0
    this%modes     = 0
    this%lambdas  = 0
    this%SigmaMat = 0
    this%LambdaMat = 0

    ! -------------------------------------------------------
    ! Creation of the Generalize eigenvalue problem Eq. (52)
    ! -------------------------------------------------------
    do i = 1 , N
    do j = 1 , N
        Mdiag(i,j) =  conjg(this%valsTau(j,i)) ! Dag of Tau
    enddo
    enddo
    ! Filling matrices
    MA  = 0
    MB  = 0
    !   (
    !   (   0       1   )
    !   (  -t^*   E-H   )
    !   (
    MA(N+1:2*N,1:N) = -Mdiag ! diag contains hermitian cojugate of Tau

    ! filling diag with diagonal matrix
    Mdiag = 0
    do i = 1 , N
        Mdiag(i,i) =  1
    enddo
    MA(1:N,N+1:2*N)     =  Mdiag
    MA(N+1:2*N,N+1:2*N) =  Mdiag*Ef - this%valsH0


    MB(1:N,1:N)         = Mdiag
    MB(N+1:2*N,N+1:2*N) = this%valsTau

    ! Ustalenie parametrow LAPACKA
    LWMAX = 20 * N
    LDVL  = 2  * N
    LDVR  = 2  * N
    ! Alokacja macierzy LAPACKA
    allocate(ALPHA(2*N))
    allocate(BETA(2*N))
    allocate(RWORK(8*N))
    allocate(WORK(LWMAX))


    LWORK = -1
    ! Initalization
    CALL ZGGEV("N","N", 2*N, MA, 2*N, MB,2*N, ALPHA,BETA, &
                DUMMY, 1, Z, LDVR, WORK, LWORK, RWORK, INFO )
    LWORK = MIN(LWMAX, INT( WORK(1)))
    ! Solving GGEV problem
    CALL ZGGEV("N","V", 2*N, MA, 2*N, MB , 2*N, ALPHA,BETA, &
                DUMMY, 1, Z, LDVR, WORK, LWORK, RWORK, INFO )
    ! Checking solution
    if( INFO /= 0 ) then
        print*,"SYS::LEAD::Cannot solve generalized eigenvalue problem for eigenmodes: ZGGEV info:",INFO
        stop
    endif

    ! -------------------------------------------------------
    ! Calculating the number of modes
    ! -------------------------------------------------------
    this%no_in_modes  = 0
    this%no_out_modes = 0
    this%no_in_em     = 0
    this%no_out_em    = 0

    do i = 1 , 2*N
        if(abs(Beta(i))>1e-16) then !
            lambda= (ALPHA(i)/BETA(i))
!            print"(i,2f10.5,A,f10.6)",i,lambda,"abs=",abs(lambda)
            c(i,:) =  Z(1:N,i)
            d(i,:) =  Z(N+1:2*N,i)
            ! Normalize vectors
            c(i,:) = c(i,:)/sqrt(sum(abs(c(i,:))**2))
            d(i,:) = d(i,:)/sqrt(sum(abs(d(i,:))**2))

            if(  abs(abs(lambda) - 1.0) < 1.0E-6 ) then ! check if propagating mode
                current = (mode_current(N,c(i,:),d(i,:),this%valsTau))
                tmpc    = current
                ! Normalize to current = 1
                c(i,:)  = c(i,:)/sqrt(abs(current))
                d(i,:)  = d(i,:)/sqrt(abs(current))
                ! Sort between incoming and outgoing
                if(current > 0) then
                    this%no_in_modes  = this%no_in_modes + 1
                    this%currents(M_IN,this%no_in_modes) = tmpc
                else
                    this%no_out_modes = this%no_out_modes + 1
                    this%currents(M_OUT,this%no_out_modes) = tmpc
                endif
            ! Evanescent modes filtering
            else if( abs(lambda) > 1.0 ) then
                this%no_out_em = this%no_out_em + 1

            else if( abs(lambda) < 1.0 .and. abs(lambda)> 1.D-16 ) then
                this%no_in_em  = this%no_in_em + 1
            else ! Strange case when lambda = 0 "standing mode" we assume that means
                 ! this case belongs to both evanescent modes
                this%no_in_em  = this%no_in_em  + 1
                this%no_out_em = this%no_out_em + 1
            endif ! end of filtering
        endif ! end of beta > 0
    enddo

    print*,"SYS::LEAD::Lead stats:"
    print*,"           No. incoming modes:",this%no_in_modes
    print*,"           No. outgoing modes:",this%no_out_modes
    print*,"           No. in. evan.modes:",this%no_in_em
    print*,"           No. out.evan.modes:",this%no_out_em


    ! Allocate T-Matrix
    allocate(this%Tnm(this%no_out_modes,this%no_out_modes))
    this%Tnm      = 0

    ! -------------------------------------------------------
    ! Filling arrays
    ! -------------------------------------------------------
    no_in    = 0
    no_e_in  = 0
    no_out   = 0
    no_e_out = 0
    do i = 1 , 2*N
        if(abs(Beta(i))>1e-16) then
            lambda= (ALPHA(i)/BETA(i))
            if(  abs(abs(lambda) - 1.0) < 1.0E-6 ) then
                current = mode_current(N,c(i,:),d(i,:),this%valsTau)
                if(current > 0) then
                    no_in  = no_in + 1
                    this%modes   (M_IN,no_in,:)  = c(i,:)
                    this%lambdas (M_IN,no_in)    = lambda
                else
                    no_out = no_out + 1
                    this%modes   (M_OUT,no_out,:) = c(i,:)
                    this%lambdas (M_OUT,no_out)   = lambda
                endif
            else if( abs(lambda) > 1.0 ) then
                    no_e_out = no_e_out + 1
                    this%modes  (M_OUT,this%no_out_modes+no_e_out,:) = c(i,:)
                    this%lambdas(M_OUT,this%no_out_modes+no_e_out)   = lambda
            else if( abs(lambda) < 1.0 .and. abs(lambda)> 1.D-16 ) then
                    no_e_in = no_e_in + 1
                    this%modes  (M_IN,this%no_in_modes+no_e_in,:) = c(i,:)
                    this%lambdas(M_IN,this%no_in_modes+no_e_in)   = lambda
            else
                    no_e_out = no_e_out + 1
                    this%modes  (M_OUT,this%no_out_modes+no_e_out,:) = c(i,:)
                    this%lambdas(M_OUT,this%no_out_modes+no_e_out)   = 1.0D20
                    no_e_in = no_e_in + 1
                    this%modes  (M_IN,this%no_in_modes+no_e_in,:)    = c(i,:)
                    this%lambdas(M_IN,this%no_in_modes+no_e_in)      = 1.0D20
                    !print*,"Problematic case! Lambda == 0 "
            endif
        endif
    enddo

    deallocate(ALPHA)
    deallocate(BETA)
    deallocate(RWORK)
    deallocate(WORK)
    deallocate(MA)
    deallocate(MB)
    deallocate(Z)
    deallocate(d)
    deallocate(c)


    ! Sorting propagating modes by the current amplitde
    call sort_modes(this%no_in_modes ,this%modes(M_IN,:,:) ,this%lambdas(M_IN,:),this%currents(M_IN,:))
    call sort_modes(this%no_out_modes,this%modes(M_OUT,:,:),this%lambdas(M_OUT,:),this%currents(M_OUT,:))
    print*,"-------------------------------------------------------"
    print*," K vec.  :    In     |       Out  |          Fluxes   "
    print*,"-------------------------------------------------------"
    do i = 1 , this%no_in_modes
        dval1 = log(this%lambdas(M_IN ,i))/II
        dval2 = log(this%lambdas(M_OUT,i))/II
        print"(A,i4,A,f10.4,A,1f10.4,A,2f10.4)","   K[",i,"]:",dval1," | ",dval2 ," | " , this%currents(:,i)
    enddo
    print*,"-------------------------------------------------------"


    ! Construction of the F^-1+ matrix Eq. (57)
    blochF = 0
    Mdiag(:,:) = this%modes(M_IN,:,:)
    call inverse_matrix(N,Mdiag)
    do k = 1, N
        do i = 1, N
        do j = 1, N
        blochF(M_IN,i,j) = blochF(M_IN,i,j) + this%lambdas(M_IN,k)**(-1)  * this%modes(M_IN,k,i) * Mdiag(j,k)
        enddo
        enddo
    enddo
    ! Construction of the F^-1+ matrix Eq. (57)
    Mdiag(:,:) = this%modes(M_OUT,:,:)
    call inverse_matrix(N,Mdiag)
    this%UTildeDagger = Mdiag
    do k = 1, N
        do i = 1, N
        do j = 1, N
        blochF(M_OUT,i,j) = blochF(M_OUT,i,j) + this%lambdas(M_OUT,k)**(-1)  * this%modes(M_OUT,k,i) * Mdiag(j,k)
        enddo
        enddo
    enddo
    ! Calculating of SigmaMatrix
    do i = 1 , N
    do j = 1 , N
        Mdiag(i,j) =  conjg(this%valsTau(j,i)) ! Dag of Tau
    enddo
    enddo
    do i = 1 , N
    do j = 1 , N
        do k = 1 , N
            this%SigmaMat(i,j) = this%SigmaMat(i,j) + Mdiag(i,k)*blochF(M_OUT,k,j)
        enddo
    enddo

    enddo
    ! add to sigma H0 internal hamiltonian
    this%SigmaMat = this%SigmaMat + this%valsH0

    ! Lambda matrix calculation:
    do i = 1 , N
    do j = 1 , N
        do k = 1 , N
            this%LambdaMat(i,j) = this%LambdaMat(i,j) + Mdiag(i,k)*(blochF(M_IN,k,j)-blochF(M_OUT,k,j))
        enddo
    enddo
    enddo

    print*,"SYS::LEAD::Modes calculated in time:", get_clock() - time , "[s]"
    deallocate(Mdiag)
    deallocate(blochF)


end subroutine calculate_modes

! ------------------------------------------------------------------------
! Calculate tranimission amplitude in the output lead
! all_atoms - reference to all atoms in the system
! phi       - calculated scattering wave function
! The result is holded: this%Tnm(1,1) element. For each incoming mode
! this number should be <= 1
! ------------------------------------------------------------------------
subroutine calculate_Tnm(this,all_atoms,phi)
    class(qlead)              :: this
    type(qatom),dimension(:)  :: all_atoms
    complex*16 ,dimension(:)  :: phi
    ! local variables
    complex*16 ,allocatable , dimension(:) :: leadPhi
    complex :: tmpT
    integer :: i,la,ls,lg

    allocate(leadPhi(this%no_sites))
    do i = 1 , this%no_sites
        la = this%l2g(i,1)
        ls = this%l2g(i,2)
        lg = all_atoms(la)%globalIDs(ls)
        leadPhi(i) = phi(lg)
    enddo
    this%Tnm = 0
    do i = 1 , this%no_out_modes
        tmpT = sum(this%UTildeDagger(:,i) * leadPhi(:))
        this%Tnm(i,1) =  abs(tmpT)**2
    enddo
    tmpT = sum(this%Tnm(:,1))
    this%Tnm(1,1) = tmpT
    deallocate(leadPhi)
endsubroutine calculate_Tnm

doubleprecision function mode_current(N,c,d,tau)
    integer :: N
    complex*16 , dimension(:)   :: c,d
    complex*16 , dimension(:,:) :: tau
    integer :: i,j
    complex*16 :: current,tmp

    current=0
    do i = 1 , N
        tmp = 0
    do j = 1 , N
        tmp = tmp + conjg(tau(j,i)) * c(j)
    enddo
        current = current +  conjg(d(i)) * tmp
    enddo
    mode_current = 2 * Imag(current)
end function mode_current

doubleprecision function mode_current_lambda(N,c,lambda,tau)
    integer :: N
    complex*16 , dimension(:)   :: c
    complex*16 :: lambda
    complex*16 , dimension(:,:) :: tau
    integer :: i,j
    complex*16 :: current,tmp

    current=0
    do i = 1 , N
        tmp = 0
    do j = 1 , N
        tmp = tmp + (tau(i,j)) * c(j)
    enddo
        current = current +  lambda*conjg(c(i)) * tmp
    enddo
    mode_current_lambda = -2 * IMAG(current)

end function mode_current_lambda

subroutine inverse_matrix(N,A)
  integer :: N
  complex*16,dimension(:,:):: A
  complex*16,allocatable,dimension(:)  :: WORK
  integer,allocatable,dimension (:)    :: IPIV
  integer info,error

  allocate(WORK(N),IPIV(N),stat=error)
  if (error.ne.0)then
    print *,"SYS::LEAD::ZGETRF::error:not enough memory"
    stop
  end if
  call ZGETRF(N,N,A,N,IPIV,info)
  if(info .eq. 0) then
!    write(*,*)"succeded"
  else
    write(*,*)"SYS::LEAD::ZGETRF::failed with info:",info
  end if
  call ZGETRI(N,A,N,IPIV,WORK,N,info)
  if(info .eq. 0) then
!    write(*,*)"succeded"
  else
   write(*,*)"SYS::LEAD::ZGETRI::failed with info:",info
  end if
  deallocate(IPIV,WORK,stat=error)
  if (error.ne.0)then
    print *,"SYS::LEAD::ZGETRF::error:fail to release"
    stop
  end if
end subroutine inverse_matrix


subroutine sort_modes(N,vectors,lambdas,currents)
    complex*16,dimension(:,:) :: vectors
    complex*16,dimension(:)   :: lambdas
    doubleprecision,dimension(:)   :: currents
    integer :: N
    complex*16,dimension(:),allocatable  :: tmpvec
    complex*16 :: tmpval
    doubleprecision :: dval , mindval , tmpcurr
    integer :: i , j  , nvec , imin

    nvec = size(vectors(1,:))
    allocate(tmpvec(nvec))
    do i = 1 , n
        imin   = i
        do j = i+1 , n
           dval     = abs(currents(j))
           mindval  = abs(currents(imin))
           if( dval > mindval ) imin = j
        enddo

        tmpcurr        = currents(i)
        currents(i)    = currents(imin)
        currents(imin) = tmpcurr

        tmpval        = lambdas(i)
        lambdas(i)    = lambdas(imin)
        lambdas(imin) = tmpval

        tmpvec          = vectors(i,:)
        vectors(i,:)    = vectors(imin,:)
        vectors(imin,:) = tmpvec
    enddo

    deallocate(tmpvec)
end subroutine sort_modes

subroutine save_lead(this,filename,ofunit)
    class(qlead) :: this
    character(*) :: filename


    integer,optional :: ofunit
    integer         :: funit = 5437629
    integer         :: i,j,id_atom_a,id_spin_a,id_atom_b,id_spin_b
    doubleprecision :: max_abs_matrix_element,normalized_value
    doubleprecision :: CUTOFF_LEVEL


    CUTOFF_LEVEL = 1.0D-10
    print*,"SYS::LEAD::Writing lead data to file:",filename


    if(present(ofunit)) then
        funit = ofunit
    else
        open(unit=funit,file=filename)
    endif


    write(funit,"(A)"),"<lead>"
    call this%lead_shape%flush_shape_data_to_file(funit)
    write(funit,"(A)"),"<vector>"
    write(funit,"(3e20.6)"),this%lead_vector
    write(funit,"(A)"),"</vector>"

    write(funit,"(A)"),"<atoms>"
    do i = 1 , this%no_sites
        id_atom_a = this%l2g(i,1)
        if(this%l2g(i,2) == 1) then
            write(funit,*),"<d>",id_atom_a,"</d>"
        endif
    enddo
    write(funit,"(A)"),"</atoms>"

    write(funit,"(A)"),"<next_atoms>"
    do i = 1 , this%no_sites
        id_atom_b = this%next_l2g(i,1)
!        print*,i,id_atom_b,this%next_l2g(i,2)
        if(this%next_l2g(i,2) == 1) then
            write(funit,*),"<d>",id_atom_b,"</d>"
        endif
    enddo
    write(funit,"(A)"),"</next_atoms>"

    ! ----------------------------------------------------------------------------
    ! COUPLING TO THE NEXT CELL DATA
    ! ----------------------------------------------------------------------------
    write(funit,"(A)"),"<lead_coupling>"
    max_abs_matrix_element = 0
    do i = 1 , this%no_sites
        do j = 1 , this%no_sites
            if( abs(this%valsTau(i,j)) > max_abs_matrix_element ) max_abs_matrix_element = abs(this%valsTau(i,j))
        enddo
    enddo
    do i = 1 , this%no_sites
        id_atom_a = this%l2g(i,1)
        id_spin_a = this%l2g(i,2)
        do j = 1 , this%no_sites
        if(this%next_l2g(j,1) == 0) then
            exit
        endif
        id_atom_b = this%next_l2g(j,1)
        id_spin_b = this%next_l2g(j,2)
        normalized_value = abs(this%valsTau(i,j))/max_abs_matrix_element
        if( normalized_value > CUTOFF_LEVEL ) then
            write(funit,"(A,4i10,A)"),"   <d>",id_atom_a,id_atom_b,id_spin_a,id_spin_b,"</d>"
        endif
        enddo
    enddo
    write(funit,"(A)"),"</lead_coupling>"

    write(funit,"(A)"),"<inner_coupling>"
    max_abs_matrix_element = 0
    do i = 1 , this%no_sites
        do j = i , this%no_sites
            if( abs(this%valsH0(i,j)) > max_abs_matrix_element ) max_abs_matrix_element = abs(this%valsH0(i,j))
        enddo
    enddo
    do i = 1 , this%no_sites
        id_atom_a = this%l2g(i,1)
        id_spin_a = this%l2g(i,2)
        do j = i , this%no_sites
        id_atom_b = this%l2g(j,1)
        id_spin_b = this%l2g(j,2)
        normalized_value = abs(this%valsH0(i,j))/max_abs_matrix_element
        if( normalized_value > CUTOFF_LEVEL ) then
            write(funit,"(A,4i10,A)"),"   <d>",id_atom_a,id_atom_b,id_spin_a,id_spin_b,"</d>"
        endif
        enddo

    enddo
    write(funit,"(A)"),"</inner_coupling>"


    write(funit,"(A)"),"</lead>"
    if(.not. present(ofunit)) close(funit)


endsubroutine save_lead

endmodule modlead
