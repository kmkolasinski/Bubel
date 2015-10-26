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
use modcommons
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
    complex*16,dimension(:,:) , allocatable :: valsS0  ! Overlaps in the lead in most cases this is diagonal matrix
    complex*16,dimension(:,:) , allocatable :: valsS1  ! Overlaps in the next unit cell in most cases this is "zero matrix"

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
    integer :: i,j,k,b,s1,s2,ns1,ns2
    integer :: no_sites,system_size, atom_id,bond_atom_id,bond_id
    integer :: irow,icol,aid,gid,sid,lid
    integer :: no_atoms
    complex :: cval
    integer,allocatable :: next_cell_atoms(:)
    doubleprecision,dimension(3) :: cellA_pos,cellB_pos,tmp_pos,cellBA_vec
    doubleprecision :: minimum_distance,dist

    ! copy shape for future possible usage
    call this%lead_shape%copyFrom(lshape)


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
    if(allocated(this%valsS0))  deallocate(this%valsS0)
    if(allocated(this%valsS1))  deallocate(this%valsS1)
    if(allocated(this%l2g))     deallocate(this%l2g)
    if(allocated(this%next_l2g))deallocate(this%next_l2g)

    allocate(this%valsH0    (no_sites,no_sites))
    allocate(this%valsTau   (no_sites,no_sites))
    allocate(this%valsS0    (no_sites,no_sites))
    allocate(this%valsS1    (no_sites,no_sites))
    allocate(this%l2g       (no_sites,2))
    allocate(this%next_l2g  (no_sites,2))

    this%valsH0   = 0
    this%valsTau  = 0
    this%valsS0   = 0
    this%valsS1   = 0
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
    print*,"minum distance=",minimum_distance

    no_atoms = 0
    do i = 1 , size(all_atoms) ! Search in all active atoms
        if(all_atoms(i)%bActive) then
            ! if atom is in the unit cell
            tmp_pos = all_atoms(i)%atom_pos - lvec
            if( lshape%is_inside(tmp_pos) == .true. ) then
                ! search for atom with the same position
                do j = 1 , no_sites

                    cellA_pos = all_atoms(this%l2g(j,1))%atom_pos
!                    print"(2i4,5f12.4)",i,j,sqrt(sum( (tmp_pos - cellA_pos )**2)),tmp_pos
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
            ns1 = size(all_atoms(atom_id)%bonds(b)%bondMatrix,1)
            ns2 = size(all_atoms(atom_id)%bonds(b)%bondMatrix,2)

            do s1 = 1 , ns1
            do s2 = 1 , ns2

            ! Filter out bond from atom A but which do not have spin this%l2g(i,1)
            if( s1 == this%l2g(i,2) ) then
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
                sid  = s2            ! and its spin
                gid  = all_atoms(aid)%globalIDs(sid) ! get global ID of site (atom+spin)
                lid  = tmp_g2l(gid) ! remap it to local ID in the lead
                irow = i
                icol = lid
                cval = all_atoms(atom_id)%bonds(b)%bondMatrix(s1,s2)

                ! fill coupling matrix
                this%valsTau(irow,icol) = cval
                this%valsS1 (irow,icol) = all_atoms(atom_id)%bonds(b)%overlapMatrix(s1,s2)

            else ! in case of coupling occures within lead

                aid  = bond_atom_id ! id of connected atom
                sid  = s2 ! and its spin ID
                gid  = all_atoms(aid)%globalIDs(sid) ! convert to global state id (atom+spin)
                lid  = tmp_g2l(gid) ! remap to local ID in lead
                irow = i
                icol = lid
                cval = all_atoms(atom_id)%bonds(b)%bondMatrix(s1,s2)
                this%valsH0(irow,icol) = cval
                this%valsS0(irow,icol) = all_atoms(atom_id)%bonds(b)%overlapMatrix(s1,s2)

!               print*,irow,icol,atom_id,aid,sqrt(sum((all_atoms(atom_id)%atom_pos-all_atoms(aid)%atom_pos)**2 ))
            endif ! end of else if "atom is in the lead or outside"
            endif ! end of if bond comes from the atom with the same spin what site has

            enddo ! end of s1
            enddo ! end of s2

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
    if(allocated(this%valsS0))  deallocate(this%valsS0)
    if(allocated(this%valsS1))  deallocate(this%valsS1)
    if(allocated(this%l2g))     deallocate(this%l2g)
    if(allocated(this%next_l2g))deallocate(this%next_l2g)

    if(allocated(this%modes))        deallocate(this%modes)
    if(allocated(this%lambdas))      deallocate(this%lambdas)
    if(allocated(this%SigmaMat))     deallocate(this%SigmaMat)
    if(allocated(this%LambdaMat))    deallocate(this%LambdaMat)
    if(allocated(this%UTildeDagger)) deallocate(this%UTildeDagger)
    if(allocated(this%Tnm))          deallocate(this%Tnm)
    if(allocated(this%currents))     deallocate(this%currents)

    call this%lead_shape%destroy_shape()

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
    INTEGER,allocatable,dimension(:)          :: ISUPPZ, IWORK , IFAIL
    DOUBLE PRECISION,allocatable,dimension(:) :: EVALS( : ), RWORK( : )
    COMPLEX*16,allocatable,dimension(:)       :: Z( :, : ), WORK( : )
    COMPLEX*16,allocatable,dimension(:,:)     :: Hamiltonian,Overlaps
    logical :: bIsIdentity
    doubleprecision :: absSum
    print*,"SYS::LEAD::Generating band structure data to file:",filename
    ! The size of the system
    N = this%no_sites


    ! Initialize lapack and allocate arrays
    LWMAX = N*24

    allocate(ISUPPZ( 2*N ))
    allocate(EVALS ( N ))
    allocate(Z( N, N ))
    allocate(IFAIL( 2*N ))

    allocate(IWORK ( LWMAX ))
    allocate(RWORK ( LWMAX ))
    allocate(WORK  ( LWMAX ))
    allocate(Hamiltonian(N,N))
    allocate(Overlaps(N,N))
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

    bIsIdentity = .true.
    ! Checking if overlap matrix is idenity matrix
    absSum      = sum(abs(this%valsS1))
    if( absSum > 1.0D-6 ) bIsIdentity = .false.
    do i = 1 , N
        if(.not. bIsIdentity) exit
    do j = 1 , N
        if(i == j .and. abs(this%valsS0(i,j)-1.0) > 1.0D-6 ) then
            bIsIdentity = .false.
            exit
        else if( i /= j .and. abs(this%valsS0(i,j)) > 1.0D-6  ) then
            bIsIdentity = .false.
            exit
        endif
    enddo
    enddo
    ! Perform scan
    open(unit = 782321, file= filename )
    do skank = kmin , kmax + dk/2 , dk

        kvec = II*skank
        do i = 1 , N
        do j = 1 , N
                Hamiltonian(i,j) = this%valsH0(i,j) + this%valsTau(i,j)*exp(kvec) + conjg(this%valsTau(j,i)*exp(kvec))
        enddo
        enddo
        Overlaps = 0
        do i = 1 , N
        do j = 1 , N
                Overlaps(i,j)    = this%valsS0(i,j) + this%valsS1(i,j)*exp(kvec) + conjg(this%valsS1(j,i)*exp(kvec))
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

        if(bIsIdentity) then
        CALL ZHEEVR( "V", 'Values', 'Lower', N, Hamiltonian, N, VL, VU, IL,&
        &             IU, ABSTOL, M, EVALS, Z, N, ISUPPZ, WORK, LWORK, RWORK,&
        &             LRWORK, IWORK, LIWORK, INFO )
        else
        CALL ZHEGVX(1,'Vectors','Values in range','Upper',N,Hamiltonian,N,Overlaps,N,VL,VU,IL, &
                      IU,ABSTOL,M,EVALS,Z,N,WORK,LWORK,RWORK, &
                      IWORK,IFAIL,INFO)
        endif
        if(M > 0) then
            write(782321,"(1000e20.10)"),skank,EVALS(1:M)
        endif

    enddo
    close(782321)

    deallocate(ISUPPZ)
    deallocate(EVALS)
    deallocate(Z)
    deallocate(IFAIL)

    deallocate(IWORK)
    deallocate(RWORK)
    deallocate(WORK)
    deallocate(Hamiltonian)
    deallocate(Overlaps)

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

    integer :: k,p,q,i,j,no_in,no_out,no_e_in,no_e_out,no_modes

    INTEGER                                      :: LDVL, LDVR , LWMAX , LWORK , INFO
    COMPLEX*16 , dimension(:) ,     allocatable  :: ALPHA , BETA , WORK
    double precision, dimension(:), allocatable  :: RWORK

    doubleprecision :: tmpc,current,time,dval1,dval2,time_calc
    COMPLEX*16 :: DUMMY(1,1),lambda,one,zero


    time = get_clock()
    N = this%no_sites
    if(QSYS_DEBUG_LEVEL > 0) then
    print*,"SYS::LEAD::Finding lead modes for N=",N
    endif
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

    ! Setting the parameters for LAPACK
    LWMAX = 40 * N
    LDVL  = 2  * N
    LDVR  = 2  * N

    allocate(ALPHA(2*N))
    allocate(BETA(2*N))
    allocate(RWORK(8*N))
    allocate(WORK(LWMAX))
    ! -------------------------------------------------------
    ! Creation of the Generalized eigenvalue problem Eq. (52)
    ! -------------------------------------------------------
!    Mdiag = -this%valsH0
!    do i = 1 , N
!!        print"(50f9.4)", dble(this%valsS0(i,:))
!        Mdiag(i,i) =   Mdiag(i,i) + Ef
!    enddo
!    open(unit=1112,file="tau.dat")
!    open(unit=1113,file="h0.dat")
!    do i = 1 , N
!        write(1112,"(500e20.8)"),this%valsTau(i,:)
!        write(1113,"(500e20.8)"),Mdiag(i,:)
!    enddo
!    close(1112)
!    close(1113)


!    call solve_GGEV(this%valsS0*Ef - this%valsH0,this%valsTau)

888 if(QSYS_USE_ZGGEV_TO_FIND_MODES) then
    do i = 1 , N
!        print"(30f7.3)",dble(this%valsTau(:,i))
    do j = 1 , N
        Mdiag(i,j) =  conjg(this%valsTau(j,i)) - Ef*conjg(this%valsS1(j,i)) ! Dag of Tau
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
    MA(N+1:2*N,N+1:2*N) =  this%valsS0*Ef - this%valsH0

    MB(1:N,1:N)         = Mdiag
    MB(N+1:2*N,N+1:2*N) = this%valsTau  - this%valsS1*Ef


!    call try_ZGGEVX(2*N,MA,MB,Z,ALPHA,BETA)

    LWORK = -1
    ! Initalization
    CALL ZGGEV("N","N", 2*N, MA, 2*N, MB,2*N, ALPHA,BETA, &
                DUMMY, 1, Z, LDVR, WORK, LWORK, RWORK, INFO )
    LWORK = MIN(LWMAX, INT( WORK(1)))
    ! Solving GGEV problem
    CALL ZGGEV("N","V", 2*N, MA, 2*N, MB , 2*N, ALPHA,BETA, &
                DUMMY, 1, Z, LDVR, WORK, LWORK, RWORK, INFO )
!     Checking solution
    if( INFO /= 0 ) then
        print*,"SYS::LEAD::Cannot solve generalized eigenvalue problem for eigenmodes: ZGGEV info:",INFO
        stop
    endif

    ! ------------------------------------------------------------
    else ! Convert Generalized eigenvalue problem to standard one
    ! ------------------------------------------------------------
    do i = 1 , N
    do j = 1 , N
        Mdiag(i,j) =  conjg(this%valsTau(j,i)) - Ef*conjg(this%valsS1(j,i)) ! Dag of Tau
    enddo
    enddo
    ! Filling matrices
    MA  = 0
    MB  = 0
    ! Invert block NxN matrix
    blochF(1,:,:) = this%valsTau  - this%valsS1*Ef
    call inverse_matrix(N,blochF(1,:,:))



    if(B_SINGULAR_MATRIX) then
        print*,"==============================================================================="
        print*,"SYS::ERROR::Lead B matrix is singular, trying to find eigen modes "
        print*,"            modes with ZGGEV solver."
        print*,"            Set QSYS_USE_ZGGEV_TO_FIND_MODES = .true. to disable this message."
        print*,"==============================================================================="
        QSYS_USE_ZGGEV_TO_FIND_MODES = .true.
        goto 888
    endif

   one  = -1
   zero = 0


   ! ------------------------------------------------------
   ! Simple apporach
   ! ------------------------------------------------------
   call ZGEMM( 'N', 'N', N, N, N, one ,blochF(1,:,:)  , N, &
                        Mdiag, N, zero,MA(N+1:2*N,1:N), N )
    ! filling diag with diagonal matrix
    Mdiag = 0
    do i = 1 , N
        Mdiag(i,i) =  1
    enddo
    MA(1:N,N+1:2*N)     =  Mdiag

    one = +1
    blochF(2,:,:) = this%valsS0*Ef - this%valsH0
    call ZGEMM( 'N', 'N', N, N, N, one ,blochF(1,:,:)  , N, &
                        blochF(2,:,:), N, zero,MA(N+1:2*N,N+1:2*N), N )

    LWORK = -1
    CALL ZGEEV( 'Not Left', 'Vectors',2*N, MA, 2*N, ALPHA, DUMMY, 1,&
                Z, LDVR, WORK, LWORK, RWORK, INFO )
    LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
    CALL ZGEEV( 'Not Left', 'Vectors', 2*N, MA, 2*N, ALPHA, DUMMY, 1, &
                Z, LDVR, WORK, LWORK, RWORK, INFO )

   ! ------------------------------------------------------
   ! Basis apporach
   ! ------------------------------------------------------
!   call ZGEMM( 'N', 'N', N, N, N, one ,blochF(1,:,:)  , N, &
!                        Mdiag, N, zero,MA(N+1:2*N,1:N), N )
!
!
!    MA(N+1:2*N,1:N) = -Mdiag
!    ! filling diag with diagonal matrix
!    MA(1:N,N+1:2*N)     =  blochF(1,:,:)
!
!    one = +1
!    blochF(2,:,:) = this%valsS0*Ef - this%valsH0
!    call ZGEMM( 'N', 'N', N, N, N, one ,blochF(2,:,:)  , N, &
!                        blochF(1,:,:), N, zero,MA(N+1:2*N,N+1:2*N), N )
!
!    LWORK = -1
!    CALL ZGEEV( 'Not Left', 'Vectors',2*N, MA, 2*N, ALPHA, DUMMY, 1,&
!                Z, LDVR, WORK, LWORK, RWORK, INFO )
!    LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!    CALL ZGEEV( 'Not Left', 'Vectors', 2*N, MA, 2*N, ALPHA, DUMMY, 1, &
!                Z, LDVR, WORK, LWORK, RWORK, INFO )

   ! ------------------------------------------------------
   ! Further
   ! ------------------------------------------------------
    ! Checking solution
    if( INFO /= 0 ) then
        print*,"SYS::LEAD::Cannot solve generalized eigenvalue problem for eigenmodes: ZGEEV info:",INFO
        stop
    endif
    BETA  = 1.0
    endif ! else choose solver SGGEV


    ! -------------------------------------------------------
    ! Calculating the number of modes
    ! -------------------------------------------------------
    this%no_in_modes  = 0
    this%no_out_modes = 0
    this%no_in_em     = 0
    this%no_out_em    = 0

! ----------------------------------------------------------------------------
! Solve problem using Wave function matching approach
! ----------------------------------------------------------------------------
    if(QSYS_SCATTERING_METHOD == QSYS_SCATTERING_WFM ) then

    do i = 1 , 2*N
        if(abs(Beta(i))>1e-16) then !
            lambda= (ALPHA(i)/BETA(i))
!            print"(i,4f16.5,A,f16.6)",i,abs(BETA(i)),abs(ALPHA(i)),lambda,"abs=",abs(lambda)
            c(i,:) =  Z(1:N,i)
            d(i,:) =  Z(N+1:2*N,i)
!            one = 1
!            call ZGEMV ( 'N', N, N, one, blochF(1,:,:), N,d(i,:),1, zero,Mdiag(1,:) , 1 )
!            d(i,:) = Mdiag(1,:)
            ! Normalize vectors
            c(i,:) = c(i,:)/sqrt(sum(abs(c(i,:))**2))
            d(i,:) = d(i,:)/sqrt(sum(abs(d(i,:))**2))
!            print*,i,sqrt(sum(abs(c(i,:))**2)),sqrt(sum(abs(d(i,:))**2))
            if(  abs(abs(lambda) - 1.0) < 1.0E-6 ) then ! check if propagating mode
                current = (mode_current(N,c(i,:),d(i,:),this%valsTau))
!                current =  mode_current_lambda(N,c(i,:),lambda,this%valsTau)
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
            else ! Strange case when lambda = 0 "standing mode" we assume that
                 ! this case belongs to both evanescent modes
                this%no_in_em  = this%no_in_em  + 1
                this%no_out_em = this%no_out_em + 1
            endif ! end of filtering
        endif ! end of beta > 0
    enddo
    if(QSYS_DEBUG_LEVEL > 0) then
    print*,"SYS::LEAD::Lead stats:"
    print*,"           No. incoming modes:",this%no_in_modes
    print*,"           No. outgoing modes:",this%no_out_modes
    print*,"           No. in. evan.modes:",this%no_in_em
    print*,"           No. out.evan.modes:",this%no_out_em
    endif

    ! checksum
    no_in = this%no_in_modes - this%no_out_modes + this%no_in_em - this%no_out_em

    if( no_in /= 0 ) then
        print*,"SYS::LEAD::ERROR::The total number of propagating modes and evanescent ones"
        print*,"           does not sum up to same number in for IN/OUT modes. The difference"
        print*,"           is following:",no_in," which means that your system is probably"
        print*,"           higlhy degenerated (symmetry reasons or badly formulated problem)."
        print*,"           The program will stop!"
        print*,"           You can try to remove the degeneracy by adding small perturbation"
        print*,"           to your system."
        stop -1
    endif


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
    k = 0
    do i = 1 , 2*N
        if(abs(Beta(i))>1e-16) then
            lambda= (ALPHA(i)/BETA(i))
            if(  abs(abs(lambda) - 1.0) < 1.0E-6 ) then
!                current = mode_current(N,c(i,:),d(i,:),this%valsTau)
                current =  mode_current_lambda(N,c(i,:),lambda,this%valsTau)
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
                    this%lambdas(M_OUT,this%no_out_modes+no_e_out)   = 1.0D10
                    no_e_in = no_e_in + 1
                    this%modes  (M_IN,this%no_in_modes+no_e_in,:)    = c(i,:)
                    this%lambdas(M_IN,this%no_in_modes+no_e_in)      = 1.0D10
                    k = k +1
!                    print*,"Problematic case! Lambda == 0 ",k
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

!    open(unit=111,file="rho.dat")
!    do i = 1 , N
!        write(111,"(500e20.6)"),abs(this%modes(M_IN,1:this%no_in_modes,i))**2,abs(this%modes(M_OUT,1:this%no_out_modes,i))**2
!    enddo
!    close(111)

    ! Sorting propagating modes by the current amplitde
    call sort_modes(this%no_in_modes ,this%modes(M_IN,:,:) ,this%lambdas(M_IN,:),this%currents(M_IN,:))
    call sort_modes(this%no_out_modes,this%modes(M_OUT,:,:),this%lambdas(M_OUT,:),this%currents(M_OUT,:))
    if(QSYS_DEBUG_LEVEL > 0) then
    print*,"-------------------------------------------------------"
    print*," K vec.  :    In     |       Out  |          Fluxes   "
    print*,"-------------------------------------------------------"
    do i = 1 , this%no_in_modes
        dval1 = log(this%lambdas(M_IN ,i))/II
        dval2 = log(this%lambdas(M_OUT,i))/II
        print"(A,i4,A,f10.4,A,1f10.4,A,2f12.8)","   K[",i,"]:",dval1," | ",dval2 ," | " , this%currents(:,i)
!        print"(A,i4,A,2f10.4,A,2f10.4,A,2f10.4)","   K[",i,"]:",this%lambdas(M_IN ,i)," | ",this%lambdas(M_OUT ,i) ," | " , this%currents(:,i)
    enddo
    print*,"-------------------------------------------------------"
    endif
    no_modes = N

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
        Mdiag(i,j) =  conjg(this%valsTau(j,i)) - Ef * conjg(this%valsS1(j,i)) ! Dag of Tau
    enddo
    enddo
    one  = 1.0
    zero = 0.0
    call ZGEMM( 'N', 'N', N, N, N, one ,Mdiag  , N, &
                        blochF(M_OUT,:,:), N, zero,this%SigmaMat, N )
    ! add to sigma H0 internal hamiltonian
    this%SigmaMat = this%SigmaMat + this%valsH0
    ! Lambda matrix calculation:
    blochF(M_OUT,:,:) = blochF(M_IN,:,:)-blochF(M_OUT,:,:)
    call ZGEMM( 'N', 'N', N, N, N, one ,Mdiag  , N, &
                        blochF(M_OUT,:,:), N, zero,this%LambdaMat, N )


! ----------------------------------------------------------------------------
! Solve problem using Quantum transmissing boundary method
! ----------------------------------------------------------------------------
    else if(QSYS_SCATTERING_METHOD == QSYS_SCATTERING_QTBM ) then

    do i = 1 , 2*N
        if(abs(Beta(i))>1e-16) then !
            lambda= (ALPHA(i)/BETA(i))
!            print"(i,4f16.5,A,f16.6)",i,abs(BETA(i)),abs(ALPHA(i)),lambda,"abs=",abs(lambda)
            c(i,:) =  Z(1:N,i)
            d(i,:) =  Z(N+1:2*N,i)
            ! Normalize vectors
            c(i,:) = c(i,:)/sqrt(sum(abs(c(i,:))**2))
            d(i,:) = d(i,:)/sqrt(sum(abs(d(i,:))**2))
!            print*,i,sqrt(sum(abs(c(i,:))**2)),sqrt(sum(abs(d(i,:))**2))
            if(  abs(abs(lambda) - 1.0) < 1.0E-6 ) then ! check if propagating mode
!                current = (mode_current(N,c(i,:),d(i,:),this%valsTau))
                current =  mode_current_lambda(N,c(i,:),lambda,this%valsTau)
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
            else ! Strange case when lambda = 0 "standing mode" we assume that
                 ! this case belongs to both evanescent modes
                !this%no_in_em  = this%no_in_em  + 1
                !this%no_out_em = this%no_out_em + 1
            endif ! end of filtering
        endif ! end of beta > 0
    enddo
    if(QSYS_DEBUG_LEVEL > 0) then
    print*,"SYS::LEAD::Lead stats:"
    print*,"           No. incoming modes:",this%no_in_modes
    print*,"           No. outgoing modes:",this%no_out_modes
    print*,"           No. in. evan.modes:",this%no_in_em
    print*,"           No. out.evan.modes:",this%no_out_em
    endif

!    ! checksum
!    no_in = this%no_in_modes - this%no_out_modes + this%no_in_em - this%no_out_em
!
!    if( no_in /= 0 ) then
!        print*,"SYS::LEAD::ERROR::The total number of propagating modes and evanescent ones"
!        print*,"           does not sum up to same number in for IN/OUT modes. The difference"
!        print*,"           is following:",no_in," which means that your system is probably"
!        print*,"           higlhy degenerated (symmetry reasons or badly formulated problem)."
!        print*,"           The program will stop!"
!        print*,"           You can try to remove the degeneracy by adding small perturbation"
!        print*,"           to your system."
!        stop -1
!    endif


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
    this%modes = 0
    do i = 1 , 2*N
        if(abs(Beta(i))>1e-16) then
            lambda= (ALPHA(i)/BETA(i))
            if(  abs(abs(lambda) - 1.0) < 1.0E-6 ) then
!                current = mode_current(N,c(i,:),d(i,:),this%valsTau)
                current =  mode_current_lambda(N,c(i,:),lambda,this%valsTau)
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
!                    no_e_out = no_e_out + 1
!                    this%modes  (M_OUT,this%no_out_modes+no_e_out,:) = c(i,:)
!                    this%lambdas(M_OUT,this%no_out_modes+no_e_out)   = 1.0D20
!                    no_e_in = no_e_in + 1
!                    this%modes  (M_IN,this%no_in_modes+no_e_in,:)    = c(i,:)
!                    this%lambdas(M_IN,this%no_in_modes+no_e_in)      = 1.0D20
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
    if(QSYS_DEBUG_LEVEL > 0) then
    print*,"-------------------------------------------------------"
    print*," K vec.  :    In     |       Out  |          Fluxes   "
    print*,"-------------------------------------------------------"
    do i = 1 , this%no_in_modes
        dval1 = log(this%lambdas(M_IN ,i))/II
        dval2 = log(this%lambdas(M_OUT,i))/II
        print"(A,i4,A,f10.4,A,1f10.4,A,2f10.4)","   K[",i,"]:",dval1," | ",dval2 ," | " , this%currents(:,i)
!        print"(A,i4,A,2f10.4,A,2f10.4,A,2f10.4)","   K[",i,"]:",this%lambdas(M_IN ,i)," | ",this%lambdas(M_OUT ,i) ," | " , this%currents(:,i)
    enddo
    print*,"-------------------------------------------------------"
    endif



    ! Calculate overlap matrix
    no_modes = this%no_out_modes + QSYS_SCATTERING_QTBM_NO_EVAN

    if(QSYS_SCATTERING_QTBM_NO_EVAN == QSYS_SCATTERING_QTBM_TAKE_ALL_EVAN) &
                no_modes = this%no_out_modes + this%no_out_em



    allocate(MA (no_modes,no_modes))
    allocate(MB (N,N))
    MA = 0
    MB = 0
    do i = 1 , no_modes
    do j = 1 , no_modes
        MA(i,j) = sum(conjg(this%modes(M_OUT,i,:))*this%modes(M_OUT,j,:))
    enddo
    enddo
    if(no_modes > 0) call inverse_matrix(no_modes,MA)

    blochF(M_OUT,:,:) = 0
    do i = 1 , N
    do j = 1 , N

    do k = 1 , no_modes
    do p = 1 , no_modes
!        MB(j,i) = MB(j,i) + (-this%lambdas(M_OUT,k)**(-1))*this%modes(M_OUT,k,j)*MA(k,p)*conjg(this%modes(M_OUT,p,i))
        MB(i,j) = MB(i,j) + (this%lambdas(M_OUT,k)-this%lambdas(M_OUT,k)**(-1)) * &
                            this%modes(M_OUT,k,i)*MA(k,p)*conjg(this%modes(M_OUT,p,j))
    enddo
    enddo
    enddo

    enddo

    blochF(M_OUT,:,:) =  MB

    ! Calculating of SigmaMatrix
    do i = 1 , N
    do j = 1 , N
        Mdiag(i,j) = conjg(this%valsTau(j,i)) - Ef * conjg(this%valsS1(j,i)) ! Dag of Tau
    enddo
    enddo

    one  = 1.0
    zero = 0.0
    call ZGEMM( 'N', 'N', N, N, N, one ,Mdiag  , N, &
                         MB, N, zero,this%SigmaMat, N )
    ! add to sigma H0 internal hamiltonian
    this%LambdaMat =  this%SigmaMat
    this%SigmaMat  = -this%SigmaMat + this%valsH0

    deallocate(this%UTildeDagger)
    allocate(this%UTildeDagger(no_modes,no_modes))
    this%UTildeDagger = MA

    deallocate(MA)
    deallocate(MB)

    ! ----------------------------------------------------------------------------
    ! Other possibilities ?
    ! ----------------------------------------------------------------------------
    else ! else of scattering method
        print*,"SYS::ERROR::Scattering problem method is not define: QSYS_SCATTERING_METHOD=",QSYS_SCATTERING_METHOD
        print*,"            QSYS_SCATTERING_METHOD=",QSYS_SCATTERING_METHOD
        stop
    endif

    if(QSYS_DEBUG_LEVEL > 0) then
    print*,"SYS::LEAD::Modes calculated in time:", get_clock() - time , "[s]"
    endif
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
    complex*16 ,allocatable , dimension(:) :: leadPhi,tmpVec
    complex*16 :: tmpT
    integer :: i,la,ls,lg,no_modes

    allocate(leadPhi(this%no_sites))
    do i = 1 , this%no_sites
        la = this%l2g(i,1)
        ls = this%l2g(i,2)
        lg = all_atoms(la)%globalIDs(ls)
        leadPhi(i) = phi(lg)
    enddo
    this%Tnm = 0
    if(QSYS_SCATTERING_METHOD == QSYS_SCATTERING_WFM) then
        do i = 1 , this%no_out_modes
            tmpT = sum(this%UTildeDagger(:,i) * leadPhi(:))
            this%Tnm(i,1) =  abs(tmpT)**2
        enddo
    else if(QSYS_SCATTERING_METHOD == QSYS_SCATTERING_QTBM) then

        no_modes = this%no_out_modes + QSYS_SCATTERING_QTBM_NO_EVAN
        if(QSYS_SCATTERING_QTBM_NO_EVAN == QSYS_SCATTERING_QTBM_TAKE_ALL_EVAN) &
               no_modes = this%no_out_modes + this%no_out_em

        allocate(tmpVec(no_modes))
        do i = 1 , no_modes
            tmpVec(i) = sum(conjg(this%modes(M_OUT,i,:))*leadPhi)
        enddo
        do i = 1 , this%no_out_modes
            tmpT = sum(this%UTildeDagger(i,:) * tmpVec(:))
            this%Tnm(i,1) =  abs(tmpT)**2
        enddo
        deallocate(tmpVec)
    endif

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


  B_SINGULAR_MATRIX = .false.
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
    write(*,*)"           It seems your matrix is singular check your code"
    B_SINGULAR_MATRIX = .true.
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


!subroutine try_ZGGEVX(N,A,B,Z,ALPHA,BETA)
!
!      integer    :: N
!      complex*16 :: A(:,:),B(:,:),Z(:,:),ALPHA(:),BETA(:)
!!  .. "Local Scalars" ..
!      INTEGER :: I, J
!      doubleprecision :: ABNRM, BBNRM
!!  .. "Local Arrays" ..
!      CHARACTER(LEN=*), PARAMETER :: FMT = '(4(1X,1H(,F7.3,1H,,F7.3,1H):))'
!
!
!        INTEGER                                   :: LDVL, LDVR , LWMAX , LWORK , INFO , LDA, LDB , ILO , IHI
!        COMPLEX*16 , dimension(:) ,  allocatable  ::  WORK
!
!        logical,allocatable :: BWORK(:)
!        integer,allocatable :: IWORK(:)
!        COMPLEX*16 :: DUMMY(1,1)
!        doubleprecision, ALLOCATABLE :: LSCALE(:), RSCALE(:), RCONDE(:), RCONDV(:) , RWORK(:)
!
!
!    ! Setting the parameters for LAPACK
!    LWMAX = 4*N**2 + 4*N
!    LDVL  = N
!    LDVR  = N
!    LDA = N
!    LDB = N
!
!    allocate(RWORK(8*N))
!    allocate(WORK(LWMAX))
!    allocate(IWORK(N+2))
!    allocate(BWORK(N))
!    ALLOCATE( LSCALE(N), RSCALE(N), RCONDE(N), RCONDV(N))
!
!    LWORK = LWMAX
!    IWORK = 6*N
!
!
!!  .. "Executable Statements" ..
!      WRITE (*,*) 'GGEVX Example Program Results'
!
!
!
!
!!      INTERFACE
!      call ZGGEVX('B','N','V','B',N,A,LDA,B,LDB,ALPHA,   &
!     &                  BETA,DUMMY,1,Z,LDVR,ILO,IHI,LSCALE,RSCALE,     &
!     &                  ABNRM,BBNRM,RCONDE,RCONDV,WORK,LWORK,RWORK,     &
!     &                  IWORK,BWORK,INFO)
!!      CHARACTER          BALANC,JOBVL,JOBVR,SENSE
!!      INTEGER            IHI,ILO,INFO,LDA,LDB,LDVL,LDVR,LWORK,N
!!      DOUBLE PRECISION   ABNRM,BBNRM
!!      LOGICAL            BWORK(*)
!!      INTEGER            IWORK(*)
!!      DOUBLE PRECISION   LSCALE(*),RCONDE(*),RCONDV(*),RSCALE(*),       &
!!     &                   RWORK(*)
!!      COMPLEX*16         A(LDA,*),ALPHA(*),B(LDB,*),BETA(*),VL(LDVL,*), &
!!     &                   VR(LDVR,*),WORK(*)
!!      END
!!      END INTERFACE
!
!      WRITE(*,*)
!      WRITE(*,*)'LSCALE : '
!      DO I=1,N
!         WRITE(*,'(F9.5)') LSCALE(I)
!      ENDDO
!      WRITE(*,*)
!      WRITE(*,*)'RSCALE : '
!      DO I=1,N
!         WRITE(*,'(F9.5)') RSCALE(I)
!      ENDDO
!      WRITE(*,*)
!      WRITE(*,*)'ABNRM = ', ABNRM
!      WRITE(*,*)
!      WRITE(*,*)'BBNRM = ', BBNRM
!      WRITE(*,*)
!      WRITE(*,*)'RCONDE : '
!      DO I=1,N
!         WRITE(*,'(F9.5)') RCONDE(I)
!      ENDDO
!      WRITE(*,*)
!      WRITE(*,*)'RCONDV : '
!      DO I=1,N
!         WRITE(*,'(F9.5)') RCONDV(I)
!      ENDDO
!
!
!      DEALLOCATE(LSCALE,BWORK,IWORK,WORK,RWORK)
!      DEALLOCATE(RSCALE, RCONDE, RCONDV)
!
!end subroutine try_ZGGEVX


subroutine ZSVD(N,A,U,S,VT)
      integer :: N
      complex*16 , dimension(:,:) :: A
      complex*16 ,allocatable , dimension(:,:) :: U,VT
      doubleprecision,allocatable,dimension(:) :: S

!*     .. Parameters ..
      INTEGER          LDA, LDU, LDVT
      INTEGER          LWMAX
      INTEGER          INFO, LWORK
      complex*16,allocatable,dimension(:) :: WORK
      doubleprecision,allocatable,dimension(:) :: RWORK

      LDA   = N
      LDU   = N
      LDVT  = N
      LWORK = -1
      LWMAX = 40*N

      allocate(U( LDU, N ), VT( LDVT, N ), S( N ),WORK(LWMAX) , RWORK(5*N))


      CALL ZGESVD( 'All', 'All', N, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK , RWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

      CALL ZGESVD( 'All', 'All', N, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK , RWORK, INFO )

      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm computing SVD failed to converge.'
         STOP
      END IF


end subroutine ZSVD

subroutine solve_GGEV(H0,Tau)
    complex*16 :: H0(:,:) , Tau(:,:)
    integer :: N,M
    complex*16 , allocatable  :: U(:,:),VT(:,:) , V(:,:) , tU(:,:),tVT(:,:)
    doubleprecision , allocatable  :: S(:) , tS(:)
    doubleprecision :: rcond,eps
    integer :: i,j
    logical :: bStabilize
    N = size(H0,1)
    print*,"Perform SVD",N

    call ZSVD(N,Tau,tU,tS,tVT)
    M = 0
    do i = 1 , N
       if(abs(tS(i)) > 1.0d-20 ) M = M+1
    enddo


    allocate(S(M),U(N,M),Vt(M,N),V(N,M))
    S = tS(1:M)
    print*,S,M
    U  = tU(1:N,1:M)
    Vt = tVt(1:M,1:N)
    V  = conjg(transpose(Vt))

    do i = 1 , N
    do j = 1 , M
        U(i,j) =   U(i,j) * sqrt( S(j) )
        V(i,j) =   V(i,j) * sqrt( S(j) )
    enddo
    enddo
    rcond = cond_number(N,tau)

    eps        = 1.0D-10
    bStabilize = .false.

    if( rcond < eps  ) bStabilize = .true.

    if(bStabilize) then
        do i = 1 , N
        do j = 1 , N
            tU(i,j) = sum(U(i,:)*conjg(U(j,:))) + sum(V(i,:)*conjg(V(j,:)))
        enddo
        enddo
        tVt = H0 + II * tU
        rcond = cond_number(N,tVt)
        if( rcond < eps ) then
            print*,"SYS::ERROR::Hoping matrix baldy defined"
            stop -1
        endif

    endif

end subroutine solve_GGEV


doubleprecision function cond_number(N,A) result(rval)
  integer :: N
  complex*16,dimension(:,:):: A
  complex*16,allocatable,dimension(:)       :: WORK
  doubleprecision,allocatable,dimension(:)  :: RWORK
  integer,allocatable,dimension (:)    :: IPIV
  integer info,error
  doubleprecision :: anorm,norm
  external zlange
  doubleprecision zlange
  B_SINGULAR_MATRIX = .false.
  allocate(WORK(2*N),IPIV(N),RWORK(2*N),stat=error)
  if (error.ne.0)then
    print *,"SYS::LEAD::ZGETRF::error:not enough memory"
    stop
  end if

  anorm = ZLANGE( '1', N, N, A, N, WORK )

  call ZGETRF(N,N,A,N,IPIV,info)
  call zgecon( '1', N, A, N, anorm, norm, work, rwork, info )
  rval = norm


  deallocate(IPIV,WORK,RWORK,stat=error)
  if (error.ne.0)then
    print *,"SYS::LEAD::ZGETRF::error:fail to release"
    stop
  end if
end function cond_number


!
!subroutine feast_dense()
!!subroutine feast_dense(N,A,B,lambdas,Z)
!
!
!
!
!!!!!!!!!!!!!!!!!! Feast declaration variable
!  integer,dimension(64) :: feastparam
!  integer :: loop
!  character(len=1) :: UPLO='F' ! 'L' or 'U' also fine
!
!!!!!!!!!!!!!!!!!! Matrix declaration variable
!  character(len=100) :: name
!  integer :: n,nnz
!  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: A
!
!!!!!!!!!!!!!!!!!! Contour
!  integer :: ccN
!  complex(kind=kind(1.0d0)),dimension(:),allocatable :: Zne, Wne, Zedge
!  integer, dimension(:), allocatable :: Nedge, Tedge
!
!!!!!!!!!!!!!!!!!! Others
!  integer :: t1,t2,tim
!  integer :: i,j,k
!  double precision :: rea,img
!
!!!!!!!!!!!!!!!!!! FEAST
!  integer :: M0,M,info
!  complex(kind=kind(1.0d0)) :: Emid
!  double precision :: r, epsout
!  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: XR ! eigenvectors
!  complex(kind=kind(1.0d0)),dimension(:),allocatable :: E ! eigenvalues
!  double precision,dimension(:),allocatable :: resr ! residual
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!Read Coordinate format and convert to dense format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  name='../../system4'
!
!  open(10,file=trim(name),status='old')
!  read(10,*) n,nnz
!  allocate(A(1:n,1:n))
!  A(1:N,1:N)=(0.0d0,0.0d0)
!  do k=1,nnz
!     read(10,*) i,j,rea,img
!     A(i,j)=rea*(1.0d0,0.0d0)+img*(0.0d0,1.0d0)
!  enddo
!  close(10)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  print *,'dense matrix -system4- size',n
!
!  call system_clock(t1,tim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! FEAST in dense format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  call feastinit(feastparam)
!  feastparam(1)=1
!  M0=50 !! M0>=M
!
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ! Create Custom Contour
!  ccN = 3     !! number of pieces that make up contour
!  allocate(Zedge(1:ccN))
!  allocate(Nedge(1:ccN))
!  allocate(Tedge(1:ccN))
!  !!! Example contour - triangle
!  Zedge = (/(0.10d0,0.410d0),(4.2d0,0.41d0),(4.2d0,-8.3d0)/)
!  Tedge(:) = (/0,0,0/)
!  Nedge(:) = (/6,6,18/)
!  !! Note: user must specify total # of contour points and edit feastparam(8)
!  feastparam(8) = sum(Nedge(1:ccN))
!  allocate(Zne(1:feastparam(8))) !! Contains the complex valued contour points
!  allocate(Wne(1:feastparam(8))) !! Contains the complex valued integrations weights
!
!  !! Fill Zne/Wne
!  print *, 'Enter FEAST'
!  call zfeast_customcontour(feastparam(8),ccN,Nedge,Tedge,Zedge,Zne,Wne)
!  print *,'---- Printing Countour Nodes ----'
!  do i=1,feastparam(8)
!     write(*,*)i,dble(Zne(i)),aimag(Zne(i))
!  enddo
!  print *,'---- Printing Contour Weights ----'
!  do i=1,feastparam(8)
!     write(*,*)i,dble(Wne(i)),aimag(Wne(i))
!  enddo
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!! ALLOCATE VARIABLE
!  allocate(E(1:M0))     ! Eigenvalue
!  allocate(XR(1:n,1:M0)) ! Right Eigenvectors ( XL = CONJG(XR) )
!  allocate(resr(1:M0))   ! Residual (if needed)
!
!!!!!!!!!!!!!  FEAST
!  print *, 'Enter FEAST'
!  call zfeast_syevx(UPLO,N,A,N,feastparam,epsout,loop,Emid,r,M0,E,XR,M,resr,info,Zne,Wne)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! POST-PROCESSING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  call system_clock(t2,tim)
!  print *,'FEAST OUTPUT INFO',info
!  if (info==0) then
!     print *,'*************************************************'
!     print *,'************** REPORT ***************************'
!     print *,'*************************************************'
!     print *,'SIMULATION TIME',(t2-t1)*1.0d0/tim
!     print *,'# mode found/subspace',M,M0
!     print *,'# iterations',loop
!     print *,'TRACE',sum(E(1:M))
!     print *,'Relative error on the Trace',epsout
!     print *,'Eigenvalues/Residuals'
!     do i=1,M
!        print *,i,E(i),resr(i)
!     enddo
!  endif
!
!!    complex*16,dimension(:,:) :: A,B,Z
!!    complex*16,dimension(:)   :: lambdas
!!    integer                   :: N
!!
!!
!!    doubleprecision :: emin , emax
!!    integer :: fpm(128),m0
!!    doubleprecision :: epsout
!!
!!    emin = 0.0
!!    emax = 1.0
!!    call feastinit(fpm)
!!    fpm(1)=1                 ! do not show any information
!!    fpm(2)=8                 ! number of contours
!!    fpm(3)=12                ! exponent of the error
!!    fpm(4)=10                ! maximum number of iteration
!!    fpm(5)=0                 ! we start with default vector (if 1 then with provided)
!!    fpm(6)=0                 ! convergece criterium with value of residuum (0 albo 1)
!!
!!
!!    m0 = 10
!!
!!
!!    call zfeast_hegv('F', N, A, N, B, N, fpm, epsout, loop, emin, emax, m0, lambdas, Z, m, res, info)
!
!
!!    class(qsys) :: sys
!!    doubleprecision   :: pEmin, pEmax
!!    integer           :: NoStates
!!    integer,optional  :: no_feast_contours,print_info,pmaks_iter
!!
!!
!!
!!    integer :: i,j,info,itmp,nw,M0,loop,no_evals,ta,ts,ns1,ns2,s1,s2
!!    integer :: no_contours,display_info,maks_iter
!!    doubleprecision :: epsout
!!    doubleprecision :: Emin, Emax
!!    doubleprecision :: time_start
!!
!!    integer,allocatable                          :: HBROWS(:)
!!    complex*16,dimension(:,:), allocatable       :: EVectors
!!    complex*16,dimension(:)  , allocatable       :: MATHVALS
!!    integer   ,dimension(:,:), allocatable       :: ROWCOLID
!!    double precision,dimension(:), allocatable   :: Evalues,Rerrors
!!
!!    integer :: NO_NON_ZERO_VALUES , NO_VARIABLES
!!
!!    if(sys%bOverlapMatrixEnabled) then
!!        print*,"==============================================================================="
!!        print*,"SYS::ERROR::Eigenvalue problem not supported for system with "
!!        print*,"            Overlap matrix different than identity matrix."
!!        print*,"==============================================================================="
!!        stop -1
!!    endif
!!    time_start = get_clock()
!!    ! Przejscie do jednostek donorowych
!!    Emin = pEmin
!!    Emax = pEmax
!!
!!
!!
!!    ! setting the default parameters
!!    if(.not. present(no_feast_contours)) then
!!        no_contours = 8
!!    else
!!        no_contours = no_feast_contours
!!    endif
!!    if(.not. present(print_info)) then
!!        display_info = 0
!!    else
!!        display_info = print_info
!!    endif
!!    if(.not. present(pmaks_iter)) then
!!        maks_iter = 20
!!    else
!!        maks_iter = pmaks_iter
!!    endif
!!
!!
!!
!!    call feastinit(fpm)
!!    fpm(1)=display_info      ! do not show any information
!!    fpm(2)=no_contours       ! number of contours
!!    fpm(3)=12                ! exponent of the error
!!    fpm(4)=maks_iter         ! maximum number of iteration
!!    fpm(5)=0                 ! we start with default vector (if 1 then with provided)
!!    fpm(6)=0                 ! convergece criterium with value of residuum (0 albo 1)
!!
!!    ! Calculate the number of non-zero elements in Hamiltonian matrix
!!    itmp = 0
!!    do i = 1, sys%no_atoms
!!        if(sys%atoms(i)%bActive) then
!!            do j = 1, sys%atoms(i)%no_bonds
!!            ns1 = size(sys%atoms(i)%bonds(j)%bondMatrix,1)
!!            ns2 = size(sys%atoms(i)%bonds(j)%bondMatrix,2)
!!!            if( abs(sys%atoms(i)%bonds(j)%bondValue) > QSYS_COUPLING_CUTOFF ) &
!!            itmp = itmp + ns1*ns2
!!            enddo
!!        endif
!!    enddo
!!    NO_NON_ZERO_VALUES = itmp
!!    ! The number of unknows is taken from the last global index
!!    NO_VARIABLES       = sys%system_size
!!    if(QSYS_DEBUG_LEVEL > 0) then
!!    print*,"SYS::Calulating eigenvalue problem using FEAST solver"
!!    endif
!!    allocate(MATHVALS(NO_NON_ZERO_VALUES))
!!    allocate(ROWCOLID(NO_NON_ZERO_VALUES,2))
!!
!!    ! Filling matrix and row-col array
!!    itmp = 0
!!    do i = 1 ,  sys%no_atoms
!!
!!        if(sys%atoms(i)%bActive) then
!!
!!        ns1   = sys%atoms(i)%no_in_states
!!        do s1 = 1 , ns1
!!
!!        do j  = 1, sys%atoms(i)%no_bonds
!!        ns2 = size(sys%atoms(i)%bonds(j)%bondMatrix,2)
!!
!!
!!        do s2 = 1 , ns2
!!            itmp = itmp + 1
!!            MATHVALS(itmp)   = sys%atoms(i)%bonds(j)%bondMatrix(s1,s2)
!!            ROWCOLID(itmp,1) = sys%atoms(i)%globalIDs(s1)
!!
!!            ta = sys%atoms(i)%bonds(j)%toAtomID
!!            ROWCOLID(itmp,2) = sys%atoms(ta)%globalIDs(s2)
!!
!!!            print"(3i4,2f10.4)",itmp,ROWCOLID(itmp,1),ROWCOLID(itmp,2),MATHVALS(itmp)
!!
!!        enddo
!!
!!        ! remove zero ellements
!!!        if( abs(MATHVALS(itmp)) < QSYS_COUPLING_CUTOFF )  itmp = itmp - 1
!!
!!        enddo
!!        enddo
!!        endif ! end of if active atom
!!    enddo
!!
!!    ! auxiliary variable
!!    nw   = NO_NON_ZERO_VALUES
!!
!!    if(display_info==1) then
!!        print*,"--------------------------------------------------"
!!        print*,"SYS::EIGENVALUE PROBLEM INPUTS"
!!        print*,"--------------------------------------------------"
!!        print*,"SYS::FEAST::Number of variables   ",NO_VARIABLES
!!        print*,"SYS::FEAST::Energy range from     ",Emin," to ",Emax
!!        print*,"SYS::FEAST::Number of contours    ",no_contours
!!        print*,"SYS::FEAST::Max. no. of iters     ",maks_iter
!!        print*,"SYS::FEAST::Expected no. of states",NoStates
!!    endif
!!
!!
!!    ! -----------------------------------------------------------
!!    allocate(HBROWS(NO_VARIABLES+1))
!!    call convert_to_HB(NO_NON_ZERO_VALUES,ROWCOLID,MATHVALS,HBROWS)
!!
!!    ! zgadujemy liczbe stanow
!!    M0  = NoStates
!!    allocate(EVectors(NO_VARIABLES,M0))
!!    allocate(Evalues(M0))
!!    allocate(Rerrors(M0))
!!
!!
!!
!!    call zfeast_hcsrev('F',&               ! - 'F' oznacza ze podawana jest pelna macierz
!!                          NO_VARIABLES,&   ! - rozmiar problemu (ile wezlow z flaga B_NORMAL)
!!                          MATHVALS(1:nw),& ! - kolejne nie zerowe wartosci w macierzy H
!!                          HBROWS,&         ! - numeracja wierszy (rodzaj zapisu macierzy rzakidch)
!!                          ROWCOLID(1:nw,2),& ! - indeksy kolumn odpowiadaja tablicy wartosci CMATA
!!                          fpm,&            ! - wektor z konfiguracja procedury
!!                          epsout,&         ! - Residuum wyjsciowe
!!                          loop, &          ! - Koncowa liczba iteracji
!!                          Emin,&           ! - Minimalna energia przeszukiwania
!!                          Emax,&           ! - Maksymalna energia
!!                          M0,&             ! - Spodziewana liczba modow w zakresie (Emin,Emax)
!!                          Evalues,&        ! - Wektor z otrzymanymi wartosciami wlasnymi
!!                          EVectors,&       ! - Macierz z wektorami (kolejne kolumny odpowiadaja kolejnym wartoscia z tablicy Evalues)
!!                          no_evals,&       ! - Liczba otrzymanych wartosci z przedziale (Emin,Emax)
!!                          Rerrors,&        ! - Wektor z bledami dla kolejnych wartosci wlasnych
!!                          info)            ! - Ewentualne informacje o bledach
!!
!!
!!
!!        if(display_info==1) then
!!            print*,"SYS::FEAST::Output error       ",  epsout
!!            print*,"SYS::FEAST::No. interations    ",  loop
!!            print*,"SYS::FEAST::No. states         ",  no_evals
!!            print*,"SYS::FEAST::Ouput info value   ",  info
!!            print*,"SYS::FEAST::Calulation time [s]",  get_clock() - time_start
!!            print*,"--------------------------------------------------"
!!        endif
!!
!!
!!
!!        sys%no_eigenvalues = no_evals
!!
!!        ! ----------------------------------------------------------------------------------
!!        ! Obsluga bledow:
!!        ! ----------------------------------------------------------------------------------
!!        selectcase(info)
!!        case( 202 )
!!            print*," Error : Problem with size of the system n (n≤0) "
!!            stop
!!        case( 201 )
!!            print*," Error : Problem with size of initial subspace m0 (m0≤0 or m0>n) "
!!            stop
!!        case( 200 )
!!            print*," Error : Problem with emin,emax (emin≥emax) "
!!            stop
!!        case(100:199)
!!            print"(A,I4,A)"," Error : Problem with ",info-100,"-th value of the input Extended Eigensolver parameter (fpm(i)). Only the parameters in use are checked. "
!!            sys%no_eigenvalues = 0
!!            stop
!!        case( 4 )
!!            print*," Warning : Successful return of only the computed subspace after call withfpm(14) = 1 "
!!            sys%no_eigenvalues = 0
!!
!!        case( 3 )
!!            print*," Warning : Size of the subspace m0 is too small (m0<m) "
!!            sys%no_eigenvalues = 0
!!
!!        case( 2 )
!!            print*," Warning : No Convergence (number of iteration loops >fpm(4))"
!!            sys%no_eigenvalues = 0
!!        case( 1 )
!!            print*," Warning : No eigenvalue found in the search interval. See remark below for further details. "
!!            sys%no_eigenvalues = 0
!!        case( 0 )
!!            print*,               "---------------------------------------------"
!!            print"(A,i12)",       " SYS::FEAST:: No. states  :",sys%no_eigenvalues
!!            print"(A,f12.3)",     "              Time [s]    :",get_clock() - time_start
!!            print"(A,e12.4)",     "              Error       :",epsout
!!            print"(A,i12)",       "              No. iters   :",loop
!!            print*,               "---------------------------------------------"
!!        case( -1 )
!!            print*," Error : Internal error for allocation memory. "
!!            stop
!!        case( -2 )
!!            print*," Error : Internal error of the inner system solver. Possible reasons: not enough memory for inner linear system solver or inconsistent input. "
!!            stop
!!        case( -3 )
!!            print*," Error : Internal error of the reduced eigenvalue solver Possible cause: matrix B may not be positive definite. It can be checked with LAPACK routines, if necessary."
!!            stop
!!        case(-199:-100)
!!            print"(A,I4,A)"," Error : Problem with the ",-info-100,"-th argument of the Extended Eigensolver interface. "
!!            stop
!!        endselect
!!
!!        ! -----------------------------------------------------------------
!!        ! Kopiowanie wynikow do odpowiednich tablic
!!        ! -----------------------------------------------------------------
!!
!!
!!        if(allocated(sys%eigenvals)) deallocate(sys%eigenvals)
!!        if(allocated(sys%eigenvecs)) deallocate(sys%eigenvecs)
!!
!!
!!        if(sys%no_eigenvalues > 0 ) then
!!            allocate(sys%eigenvecs(NO_VARIABLES,sys%no_eigenvalues))
!!
!!            sys%eigenvecs(:,1:sys%no_eigenvalues) = EVectors(:,1:sys%no_eigenvalues)
!!            allocate(sys%eigenvals(sys%no_eigenvalues))
!!            sys%eigenvals(1:sys%no_eigenvalues) = Evalues(1:sys%no_eigenvalues)
!!
!!        endif
!!
!!
!!        deallocate(MATHVALS)
!!        deallocate(ROWCOLID)
!!        deallocate(HBROWS)
!!        deallocate(EVectors)
!!        deallocate(Evalues)
!!        deallocate(Rerrors)
!!        if(QSYS_DEBUG_LEVEL > 0) then
!!        print*,"SYS::Eigenvalues calculated. Found:",sys%no_eigenvalues
!!        endif
!end subroutine feast_dense
!

endmodule modlead
