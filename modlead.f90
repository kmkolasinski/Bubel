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
include 'lapack.f90'
module modlead
use modshape
use modsys
use modcommons
use modutils
use mkl95_lapack

implicit none
private
complex*16,parameter :: II = cmplx(0.0D0,1.0D0)
integer ,parameter :: M_IN = 1 , M_OUT = 2

! create sparse version of tau matrix
complex*16, allocatable :: sparse_tau_vals(:)
integer,   allocatable  :: sparse_tau_rcid(:,:)
integer                 :: sparse_tau_novals

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
    integer :: lead_type
    complex*16      :: pseudo_lead_phase
    doubleprecision :: modeT  ! total transmission for a given incoming mode
    doubleprecision :: totalT ! total transmission summed over all incoming modes. For input lead this variable
                              ! contains the reflection probability
    ! ----------------------------------------------------------------------
    ! Modes parameters:
    ! According to this paper:http://www.psi-k.org/newsletters/News_80/Highlight_80.pdf
    ! see equation: (53) Arrays below are filled after calling calculate_modes function
    ! ----------------------------------------------------------------------
    complex*16,allocatable :: lambdas(:,:)             ! lambda(M_IN,:) => \lambda_+ and lambda(M_OUT,:) => \lambda_-
    complex*16,allocatable :: modes (:,:,:)            ! modes({M_IN,M_OUT},M,:) N-th mode vector
    complex*16,allocatable :: SigmaMat(:,:)            ! See Eq. (68) here we added H0 to this matrix
    complex*16,allocatable :: LambdaMat(:,:)           ! See Eq. (67) is the definition of the Qm matrix.
    doubleprecision,allocatable :: Tnm(:,:),Rnm(:,:)
    doubleprecision,allocatable :: currents(:,:)       ! Fluxes carried by mode M currents({M_IN,M_OUT},M) - note modes are sorted
                                                       ! from largest current to smallest
    complex*16,allocatable      :: UTildeDagger(:,:)   ! Dual Vectors matrix see Eq. (55)


    integer,dimension(:,:),allocatable        ::      l2g ! mapping from local id in lead to global id (:,1)  = atom_ID , (:,2) = spin_ID
    integer,dimension(:,:),allocatable        :: next_l2g ! the same as above but mapping to the atoms in the next unit cell
    integer,dimension(:,:,:),allocatable      :: next_update_map ! contains all indices to update lead matrices fastly
    integer,dimension(:,:,:),allocatable      :: update_map

    logical :: LEAD_BAD_NEARST       = .false.
    logical :: bUseShurDecomposition = .false.
    contains
    procedure,pass(this) :: init_lead!()
    procedure,pass(this) :: update_lead!(this,all_atoms)
    procedure,pass(this) :: bands!(this,filename,kmin,kmax,dk,Emin,Emax)
    procedure,pass(this) :: destroy!()
    procedure,pass(this) :: save_lead!(this,filename,funit)
    procedure,pass(this) :: calculate_modes!(this,Ef)
    procedure,pass(this) :: calculate_Tnm!(this,all_atoms,n,phi)
    procedure,pass(this) :: diagonalize_currents!(this,all_atoms,n,phi)
    procedure,pass(this) :: calc_average_spins !(this,dir,spins) - calculate averago polariaztions in XYZ for systems with s=1/2

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
    this%lead_type   =  QSYS_LEAD_TYPE_NORMAL
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
    if(allocated(this%next_update_map))deallocate(this%next_update_map)
    if(allocated(this%update_map))deallocate(this%update_map)

    allocate(this%valsH0    (no_sites,no_sites))
    allocate(this%valsTau   (no_sites,no_sites))
    allocate(this%valsS0    (no_sites,no_sites))
    allocate(this%valsS1    (no_sites,no_sites))
    allocate(this%l2g       (no_sites,2))
    allocate(this%next_l2g  (no_sites,2))
    allocate(this%next_update_map(no_sites,no_sites,4))
    allocate(this%update_map(no_sites,no_sites,4))

    this%valsH0   = 0
    this%valsTau  = 0
    this%valsS0   = 0
    this%valsS1   = 0
    this%l2g      = 0
    this%next_l2g = 0

    this%next_update_map = -1 ! set null poiners
    this%update_map      = -1 ! set null poiners
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
            endif

            if( lshape%is_inside(all_atoms(i)%atom_pos) == .true. ) then
                ! create mapping between localID (l2g array) and the atom (j) and its spin (j)
                do j = 1 , all_atoms(i)%no_in_states
                    no_sites = no_sites + 1
                    this%l2g(no_sites,1) = i ! } mapping
                    this%l2g(no_sites,2) = j ! }
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

    no_atoms = 0
    do i = 1 , size(all_atoms) ! Search in all active atoms
        if(all_atoms(i)%bActive) then
            ! if atom is in the unit cell
            tmp_pos = all_atoms(i)%atom_pos - lvec
            if( lshape%is_inside(tmp_pos) == .true. ) then
                ! search for atom with the same position
                do j = 1 , no_sites
                    cellA_pos = all_atoms(this%l2g(j,1))%atom_pos
                    if( sqrt(sum( (tmp_pos - cellA_pos)**2 )) < minimum_distance*1.0D-5 ) then
                        exit
                    endif
                enddo
                if( j > no_sites ) then
                    print"(A)",              " SYS::LEAD::ERROR::The translation"
                    print"(A,i9,A,3e12.4,A)","                   at atom:",i," with position r=(",all_atoms(i)%atom_pos,")"
                    print"(A,3e12.4,A)"     ,"                   with translation vector: v=(",lvec,")"
                    print"(A)",              "                   is impossible by applying lead translation vector. Maybe some of the "
                    print"(A)",              "                   atoms have not been inlcuded in the lead. Check lead area and T vector."
!                    stop -1
                    cycle
                endif
                ! j keeps local ID of the same atom in unit
                do k = 1 , all_atoms(this%l2g(j,1))%no_in_states
                    b = tmp_g2l(all_atoms(this%l2g(j,1))%globalIDs(k))
                    this%next_l2g(b,1) = i
                    this%next_l2g(b,2) = this%l2g(b,2)
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
                this%next_update_map(irow,icol,:) = (/ atom_id , b , s1 ,s2 /)
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
                this%update_map(irow,icol,:) = (/ atom_id , b , s1 ,s2 /)
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




subroutine update_lead(this,all_atoms)
    class(qlead) :: this
    type(qatom),dimension(:) :: all_atoms
    integer    :: i,j,a,b,s1,s2,no_sites

    no_sites = this%no_sites

    this%valsH0  = 0
    this%valsTau = 0
    this%valsS1  = 0
    this%valsS0  = 0
    do i = 1 , no_sites
    do j = 1 , no_sites
        if(this%update_map(i,j,1) /= -1)then
            a  = this%update_map(i,j,1)
            b  = this%update_map(i,j,2)
            s1 = this%update_map(i,j,3)
            s2 = this%update_map(i,j,4)
            this%valsH0(i,j) = all_atoms(a)%bonds(b)%bondMatrix   (s1,s2)
            this%valsS0(i,j) = all_atoms(a)%bonds(b)%overlapMatrix(s1,s2)
        endif
        if(this%next_update_map(i,j,1) /= -1)then
            a  = this%next_update_map(i,j,1)
            b  = this%next_update_map(i,j,2)
            s1 = this%next_update_map(i,j,3)
            s2 = this%next_update_map(i,j,4)
            this%valsTau(i,j) = all_atoms(a)%bonds(b)%bondMatrix   (s1,s2)
            this%valsS1 (i,j) = all_atoms(a)%bonds(b)%overlapMatrix(s1,s2)
        endif
    enddo ! end loop j
    enddo ! end loop i

end subroutine update_lead

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
    if(allocated(this%next_update_map))deallocate(this%next_update_map)
    if(allocated(this%update_map))     deallocate(this%update_map)

    if(allocated(this%modes))        deallocate(this%modes)
    if(allocated(this%lambdas))      deallocate(this%lambdas)
    if(allocated(this%SigmaMat))     deallocate(this%SigmaMat)
    if(allocated(this%LambdaMat))    deallocate(this%LambdaMat)
    if(allocated(this%UTildeDagger)) deallocate(this%UTildeDagger)
    if(allocated(this%Tnm))          deallocate(this%Tnm)
    if(allocated(this%Rnm))          deallocate(this%Rnm)
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
    complex*16  , allocatable , dimension(:,:)    :: MA,MB,Z,Qshur,Zshur,Sshur,Pshur
    complex*16, allocatable , dimension(:,:)      :: Mdiag,d,c,Z11,Z21
    complex*16, allocatable , dimension(:,:,:)    :: blochF

    integer :: k,p,i,j,no_in,no_out,no_e_in,no_e_out,no_modes,no_inf_modes

    INTEGER                                      :: LDVL, LDVR , LWMAX , LWORK , INFO
    COMPLEX*16 , dimension(:) ,     allocatable  :: ALPHA , BETA , WORK
    double precision, dimension(:), allocatable  :: RWORK
    logical, dimension(:), allocatable    :: iselect

    doubleprecision :: tmpc,current,time,dval1,dval2
    COMPLEX*16 :: DUMMY(1,1),lambda,one,zero
    logical :: bShurDecompositionForStandEP


    bShurDecompositionForStandEP = .false.
    this%bUseShurDecomposition   = .false.

    time = get_clock()
    N = this%no_sites
    if(QSYS_DEBUG_LEVEL > 0) then
    print*,"SYS::LEAD::Finding lead modes for N=",N
    endif


    ! Create sparse version of tau matrix: since current has to
    ! be counted many times it is good to use the fact that tau
    ! is usually very sparse matrix
    if(allocated(sparse_tau_vals))deallocate(sparse_tau_vals)
    if(allocated(sparse_tau_rcid))deallocate(sparse_tau_rcid)
    allocate(sparse_tau_vals(N*N))
    allocate(sparse_tau_rcid(N*N,2))
    sparse_tau_novals = 0
    do i = 1 , N
    do j = 1 , N
        lambda = this%valsTau(i,j)
        if(abs(lambda) > qsys_double_error ) then
            sparse_tau_novals = sparse_tau_novals + 1
            k = sparse_tau_novals
            sparse_tau_vals(k)   = lambda
            sparse_tau_rcid(k,1) = i
            sparse_tau_rcid(k,2) = j
        endif
    enddo
    enddo

    ! ---------------------------------------
    ! Memory allocations
    ! ---------------------------------------
    if(allocated(this%modes))       deallocate(this%modes)
    if(allocated(this%lambdas))     deallocate(this%lambdas)
    if(allocated(this%SigmaMat))    deallocate(this%SigmaMat)
    if(allocated(this%LambdaMat))   deallocate(this%LambdaMat)
    if(allocated(this%Tnm))         deallocate(this%Tnm)
    if(allocated(this%Rnm))         deallocate(this%Rnm)
    if(allocated(this%currents))    deallocate(this%currents)
    if(allocated(this%UTildeDagger))deallocate(this%UTildeDagger)

    allocate(this%modes     (2,N,N))
    allocate(this%lambdas   (2,N))
    allocate(this%SigmaMat  (N,N))
    allocate(this%LambdaMat (N,N))
    allocate(this%currents  (2,N))
    allocate(this%UTildeDagger(N,N))

    this%currents       = 0
    this%UTildeDagger   = 0
    this%modes          = 0
    this%lambdas        = 0
    this%SigmaMat       = 0
    this%LambdaMat      = 0


    if(this%lead_type == QSYS_LEAD_TYPE_PSEUDO_TRANSPARENT) then
        allocate(Mdiag  (N,N))
        allocate(this%Tnm(1,1))
        do i = 1 , N
        do j = 1 , N
            Mdiag(i,j) =  conjg(this%valsTau(j,i)) - Ef * conjg(this%valsS1(j,i)) ! Dag of Tau
        enddo
        enddo

        Mdiag = this%pseudo_lead_phase*Mdiag
        this%SigmaMat = Mdiag + this%valsH0
        deallocate(Mdiag)
        return
    endif


    ! Setting the parameters for LAPACK
    LWMAX = 40 * N
    LDVL  = 2  * N
    LDVR  = 2  * N

    allocate(ALPHA  (2*N))
    allocate(BETA   (2*N))
    allocate(RWORK  (8*N))
    allocate(WORK   (LWMAX))
    allocate(MA     (2*N,2*N))
    allocate(MB     (2*N,2*N))
    allocate(Z      (2*N,2*N))
    allocate(Mdiag  (N,N))
    allocate(blochF (2,N,N))
    allocate(d      (2*N,N))
    allocate(c      (2*N,N))
    allocate(Z11 (N,N))
    allocate(Z21 (N,N))
    allocate(iselect(2*N))

    ! -------------------------------------------------------
    ! Creation of the Generalized eigenvalue problem Eq. (52)
    ! Use this method when coupling matrix is singular.
    ! -------------------------------------------------------
888 if(QSYS_USE_ZGGEV_TO_FIND_MODES) then

    allocate(Qshur      (2*N,2*N))
    allocate(Zshur      (2*N,2*N))
    allocate(Sshur      (2*N,2*N))
    allocate(Pshur      (2*N,2*N))

    do i = 1 , N
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
    MA(1:N,N+1:2*N)     = Mdiag
    MA(N+1:2*N,N+1:2*N) = this%valsS0*Ef - this%valsH0
    MB(1:N,1:N)         = Mdiag
    MB(N+1:2*N,N+1:2*N) = this%valsTau  - this%valsS1*Ef

    Sshur = MA
    Pshur = MB
    ! ---------------------------------------------------------------
    ! Shur decompostion
    ! ---------------------------------------------------------------
    this%bUseShurDecomposition = .true.
    if(QSYS_FORCE_SCHUR_DECOMPOSITION) then
        call gges(Sshur, Pshur, alpha, beta , vsl=Qshur,vsr=Zshur,info=info)
        if( INFO /= 0 ) then
            print*,"SYS::LEAD::Shur decomposition failed: gges info:",INFO
            stop
        endif
        ! Calculate eigen vectors of translation operator
        Z  = Zshur
        call tgevc(Sshur, Pshur ,howmny = 'B' ,vr=Z ,info=info)
        if( INFO /= 0 ) then
            print*,"SYS::LEAD::Shur decomposition failed: tgevc info:",INFO
            stop
        endif
    else ! end of else force shur
        call GGEV(MA, MB, alpha, beta , vr = Z,info=info)
        if( INFO /= 0 ) then
            print*,"SYS::LEAD::Cannot solve generalized eigenvalue problem for eigenmodes: ZGGEV info:",INFO
            stop
        endif
    endif

         ! ------------------------------------------------------------
    else ! Try to convert Generalized eigenvalue problem to standard one
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
    Z11 = this%valsTau  - this%valsS1*Ef
    call inverse_matrix(N,Z11)
    blochF(1,:,:) = Z11

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
   zero =  0
   ! ------------------------------------------------------
   ! Simple apporach
   ! ------------------------------------------------------
   call ZGEMM( 'N', 'N', N, N, N, one ,Z11 , N , &
                        Mdiag, N, zero,Z21 , N )
    MA(N+1:2*N,1:N) = Z21
    ! filling diag with diagonal matrix
    Mdiag = 0
    do i = 1 , N
        Mdiag(i,i) =  1
    enddo
    MA(1:N,N+1:2*N)     =  Mdiag


    one = +1
    Z21 = this%valsS0*Ef - this%valsH0

    call ZGEMM( 'N', 'N', N, N, N, one ,Z11  , N, &
                        Z21, N, zero,Mdiag, N )
    MA(N+1:2*N,N+1:2*N) = Mdiag

    ! Solve eigen problem normally or use Schur decomposition
    ! for more stable results
    CONDA = 1.0/cond_number(N,blochF(1,:,:))
    if(CONDA*qsys_double_error < 1.0D-6 ) then
        LWORK = -1
        CALL ZGEEV( 'Not Left', 'Vectors',2*N, MA, 2*N, ALPHA, DUMMY, 1,&
                    Z, LDVR, WORK, LWORK, RWORK, INFO )
        LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
        CALL ZGEEV( 'Not Left', 'Vectors', 2*N, MA, 2*N, ALPHA, DUMMY, 1, &
                    Z, LDVR, WORK, LWORK, RWORK, INFO )
        ! Checking solution
        if( INFO /= 0 ) then
            print*,"SYS::LEAD::Cannot solve generalized eigenvalue problem for eigenmodes: ZGEEV info:",INFO
            stop
        endif
    else ! use shur factorization to calculate lead self energy
        this%bUseShurDecomposition   = .true.
        bShurDecompositionForStandEP = .true.
        allocate(Zshur      (2*N,2*N))
        allocate(Sshur      (2*N,2*N))
        Sshur = MA
        call gees(Sshur, ALPHA ,vs=Zshur,info=info)
        ! Checking solution
        if( INFO /= 0 ) then
            print*,"SYS::LEAD::Shur decomposition failed: gees info:",INFO
            stop
        endif
        Z      = Zshur
        call trevc(Sshur, howmny='B',vr=Z ,info=info)
        if( INFO /= 0 ) then
            print*,"SYS::LEAD::Shur decomposition failed: trevc info:",INFO
            stop
        endif
    ! ---------------------------------------------------
    endif ! choose between shur or normal diagonalization

    BETA  = 1.0
    endif ! else choose solver SGGEV
    ! -------------------------------------------------------
    ! Calculating the number of modes
    ! -------------------------------------------------------
    this%no_in_modes  = 0
    this%no_out_modes = 0
    this%no_in_em     = 0
    this%no_out_em    = 0

    ! --------------------------------------------------------
    ! Solve problem using Wave function matching approach
    ! --------------------------------------------------------
    if(QSYS_SCATTERING_METHOD == QSYS_SCATTERING_WFM ) then
    iselect = .false. ! reordeing by Shur decomposition
    no_inf_modes = 0
    do i = 1 , 2*N

        lambda = (ALPHA(i)/BETA(i))
        c(i,:) = Z(1:N,i)
        d(i,:) = Z(N+1:2*N,i)

        if(abs(Beta(i))>1e-16) then !
            ! Normalization make sense only for valid states ( lambda != +inf )
            c(i,:) = c(i,:)/sqrt(sum(abs(c(i,:))**2))
            d(i,:) = d(i,:)/sqrt(sum(abs(d(i,:))**2))
            ! Normalize vectors
            if(  abs(abs(lambda) - 1.0) < 1.0E-6 ) then ! check if propagating mode
                current = (mode_current_sparse(c(i,:),lambda,sparse_tau_vals,sparse_tau_rcid,sparse_tau_novals))
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
                    iselect(i) = .true.
                endif
            ! Evanescent modes filtering
            else if( abs(lambda) > 1.0 ) then
                this%no_out_em = this%no_out_em + 1
                iselect(i) = .true.
            else if( abs(lambda) < 1.0 .and. abs(lambda)> 1.D-16 ) then
                this%no_in_em  = this%no_in_em + 1
            else ! Strange case when lambda = 0 "standing mode" we assume that
                 ! this case belongs to both evanescent modes
                this%no_in_em  = this%no_in_em  + 1
            endif ! end of filtering
        else  ! else of beta > 0
            this%no_out_em = this%no_out_em + 1
            iselect(i) = .true.
            no_inf_modes = no_inf_modes + 1
        endif ! end of beta > 0
    enddo


    ! Allocate T-Matrix
    allocate(this%Tnm(this%no_out_modes,this%no_out_modes))
    allocate(this%Rnm(this%no_out_modes,this%no_out_modes))
    this%Tnm      = 0
    this%Rnm      = 0
    ! -------------------------------------------------------
    ! Filling arrays: modes, lambdas...
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
                current = (mode_current_sparse(c(i,:),lambda,sparse_tau_vals,sparse_tau_rcid,sparse_tau_novals))
                if(current > 0) then
                    no_in  = no_in + 1
                    this%modes   (M_IN,no_in,:)  = c(i,:)
                    this%lambdas (M_IN,no_in)    = lambda
                else
                    no_out = no_out + 1
                    this%modes   (M_OUT,no_out,:) = c(i,:)
                    this%lambdas (M_OUT,no_out)   = lambda

                    Z11(:,no_out)         = Z(1:N,i)
                    Z21(:,no_out)         = Z(1+N:2*N,i)
                    blochF(M_IN,1,no_out) = lambda
                endif
            else if( abs(lambda) > 1.0 ) then
                    no_e_out = no_e_out + 1
                    this%modes  (M_OUT,this%no_out_modes+no_e_out,:) = c(i,:)
                    this%lambdas(M_OUT,this%no_out_modes+no_e_out)   = lambda

                    Z11(:,this%no_out_modes+no_e_out) = Z(1:N,i)
                    Z21(:,this%no_out_modes+no_e_out) = Z(1+N:2*N,i)
                    blochF(M_IN,1,this%no_out_modes+no_e_out) = lambda

            else if( abs(lambda) < 1.0 .and. abs(lambda)> 1.D-16 ) then
                    no_e_in = no_e_in + 1
                    this%modes  (M_IN,this%no_in_modes+no_e_in,:)    = c(i,:)
                    this%lambdas(M_IN,this%no_in_modes+no_e_in)      = lambda
            else
                    no_e_in = no_e_in + 1
                    this%modes  (M_IN,this%no_in_modes+no_e_in,:)    = c(i,:)
                    this%lambdas(M_IN,this%no_in_modes+no_e_in)      = 0.0D0
            endif
        else ! else beta > 0
                    no_e_out = no_e_out + 1
                    this%modes  (M_OUT,this%no_out_modes+no_e_out,:) = c(i,:)
                    this%lambdas(M_OUT,this%no_out_modes+no_e_out)   = 1.0D50

                    Z11(:,this%no_out_modes+no_e_out)         = Z(1:N,i)
                    Z21(:,this%no_out_modes+no_e_out)         = Z(1+N:2*N,i)
                    blochF(M_IN,1,this%no_out_modes+no_e_out) = 1.0D50
        endif
    enddo

    ! Sorting propagating modes by the current amplitde
    ! and diagonalizing current operator for degenerated states
    call sort_modes(this%no_in_modes ,this%modes(M_IN,:,:) ,this%lambdas(M_IN,:),this%currents(M_IN,:))
    lambda = this%lambdas(M_IN,1)
    j      = 1
    do i = 2 ,this%no_in_modes
        if( abs(lambda -this%lambdas(M_IN,i)) > 1.0D-10  ) then
            call this%diagonalize_currents(M_IN,j,i-1)
            j = i
            lambda = this%lambdas(M_IN,i)
        endif
    enddo
    call this%diagonalize_currents(M_IN,j,i-1)

    ! Sorting propagating modes by the current amplitde
    ! and diagonalizing current operator for degenerated states
    call sort_modes(this%no_out_modes,this%modes(M_OUT,:,:),this%lambdas(M_OUT,:),this%currents(M_OUT,:))
    lambda = this%lambdas(M_OUT,1)
    j      = 1
    do i = 2 ,this%no_out_modes
        if( abs(lambda - this%lambdas(M_OUT,i)) > 1.0D-10  ) then
            call this%diagonalize_currents(M_OUT,j,i-1)
            j = i
            lambda = this%lambdas(M_OUT,i)
        endif
    enddo
    call this%diagonalize_currents(M_OUT,j,i-1)


    if(.not. this%bUseShurDecomposition) then
        ! Construction of the F^-1+ matrix Eq. (57)
        Mdiag(:,:) = this%modes(M_OUT,:,:)
        call inverse_matrix(N,Mdiag)
        this%UTildeDagger  = Mdiag
        blochF(M_OUT,:,:)  = 0
        ! Calculate inverse of lambda once
        do k = 1 , N
            Z11(1,k) = this%lambdas(M_OUT,k)**(-1)
        enddo
        do i = 1, N
        do j = 1, N
             Mdiag(i,j) = Mdiag(i,j)*Z11(1,j)
        enddo
        enddo
        one  = 1.0
        zero = 0.0
        ! Use blas for matrix multiplication

        Z11 = this%modes(M_OUT,:,:)
        call ZGEMM('T','T',N,N,N,one,Z11,N,Mdiag,N,zero,Z21,N)
        blochF(M_OUT,:,:) = Z21
    else
        ! Reorded Shur Unitary matrix - Z

        if(.not. QSYS_FORCE_SCHUR_DECOMPOSITION) then
            ! Use singular value decomposition to inverse matrices
            ! and remove vectors with lambda = infty
            ! Z21   - U
            ! RWORK - S
            ! Mdiag - VT

            ALPHA(1:N)        = blochF(M_IN,1,:)
!            blochF(M_OUT,:,:) = Z11
            call ZSVD(N,Z11,Z21,RWORK,Mdiag)

            blochF(M_OUT,:,:) = 0
            do i = 1, N
            do j = 1, N
                do p = 1 , N - no_inf_modes
                    blochF(M_OUT,i,j) = blochF(M_OUT,i,j) + Mdiag(i,p)*conjg(Mdiag(j,p))/ALPHA(p)
                enddo
            enddo
            enddo
            do i = 1, N - no_inf_modes
            do j = 1, N - no_inf_modes
                blochF(M_OUT,i,j) = blochF(M_OUT,i,j)*RWORK(i)/RWORK(j)
            enddo
            enddo
            one  = 1.0
            zero = 0.0
            Z11  = 0.0
            Mdiag = blochF(M_OUT,:,:)
            call ZGEMM( 'N', 'N', N,N - no_inf_modes,N - no_inf_modes, one , Z21,N, Mdiag, N, zero,Z11,N)
            blochF(M_OUT,:,:) = 0
            call ZGEMM( 'N', 'C', N,N,N - no_inf_modes, one , Z11,N,Z21,N, zero,Mdiag,N)
            blochF(M_OUT,:,:) = Mdiag
        else ! force shur method

            Z  = Zshur
            if(bShurDecompositionForStandEP) then
                call trsen(Sshur, select=iselect,q=Z ,info=info)
            else
                call tgsen(a=Sshur,b=Pshur, select=iselect,z=Z ,info=info)
            endif
            if( INFO /= 0 ) then
                print*,"SYS::LEAD::Shur decomposition failed: tgsen/trsen info:",INFO
                stop
            endif
            ! Take block elements - Michael Wimmer work
            Z11(:,:) = Z(1:N,1:N)
            Z21(:,:) = Z(N+1:2*N,1:N)
            ! Try to inverse - if this will fail, use QTBM
            call inverse_matrix(N,Z21)
            blochF(M_OUT,:,:) = 0
            do i = 1, N
            do j = 1, N
                do k = 1, N
                    blochF(M_OUT,i,j) = blochF(M_OUT,i,j) +   Z11(i,k) * Z21(k,j)
                enddo
            enddo
            enddo
        endif
    endif


    ! Calculating of SigmaMatrix
    do i = 1 , N
    do j = 1 , N
        Mdiag(i,j) =  conjg(this%valsTau(j,i)) - Ef * conjg(this%valsS1(j,i)) ! Dag of Tau
    enddo
    enddo
    one  = 1.0
    zero = 0.0
    Z11 = blochF(M_OUT,:,:)
    call ZGEMM( 'N', 'N', N, N, N, one ,Mdiag  , N, Z11, N, zero,this%SigmaMat, N )

    ! add to sigma H0 internal hamiltonian
    this%SigmaMat = this%SigmaMat + this%valsH0

    ! Lambda matrix calculation:
    call ZGEMM( 'N', 'N', N, N, N, one ,Mdiag  , N, Z11, N, zero,this%LambdaMat, N )


    if(this%bUseShurDecomposition) then

        ! -------------------------------------------------
        ! Calculate overlap matrix
        ! -------------------------------------------------
        no_modes = this%no_out_modes + this%no_out_em - no_inf_modes

        deallocate(MA)
        allocate(MA (no_modes,N))
        deallocate(this%UTildeDagger)
        allocate(this%UTildeDagger(no_modes,no_modes))


        MA  = conjg((this%modes(M_OUT,1:no_modes,:)))
        Z11 = this%modes(M_OUT,:,:)
        call ZGEMM( 'N', 'T', no_modes,no_modes,N,one, &
                    MA, no_modes, Z11, N, zero,this%UTildeDagger, no_modes )


!        MA = transpose(this%modes(M_OUT,:,:))
!        call ZGEMM( 'C', 'T', no_modes,no_modes,no_modes, one , &
!                    MA, no_modes, this%modes(M_OUT,:,:), no_modes, zero,this%UTildeDagger, no_modes )
        if(no_modes > 0) call inverse_matrix(no_modes,this%UTildeDagger)


        if(allocated(Qshur))deallocate(Qshur)
        if(allocated(Zshur))deallocate(Zshur)
        if(allocated(Sshur))deallocate(Sshur)
        if(allocated(Pshur))deallocate(Pshur)

    endif


! ----------------------------------------------------------------------------
! Solve problem using Quantum transmissing boundary method
! ----------------------------------------------------------------------------
    else if(QSYS_SCATTERING_METHOD == QSYS_SCATTERING_QTBM ) then
    no_inf_modes = 0
    one          = 1.0
    zero         = 0.0


    do i = 1 , 2*N

        if(abs(Beta(i))>1e-16) then !
            lambda= (ALPHA(i)/BETA(i))
            c(i,:) =  Z(1:N,i)
            d(i,:) =  Z(N+1:2*N,i)
            ! Normalize vectors
            c(i,:) = c(i,:)/sqrt(sum(abs(c(i,:))**2))
            d(i,:) = d(i,:)/sqrt(sum(abs(d(i,:))**2))
            if(  abs(abs(lambda) - 1.0) < 1.0E-6 ) then ! check if propagating mode
                current = (mode_current_sparse(c(i,:),lambda,sparse_tau_vals,sparse_tau_rcid,sparse_tau_novals))
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
                ! Save for sorting reasons
            else if( abs(lambda) < 1.0 .and. abs(lambda)> 1.D-16 ) then
                this%no_in_em  = this%no_in_em + 1
            endif ! end of filtering
        else ! calculate number of states with lambda = infty
            no_inf_modes = no_inf_modes + 1
        endif ! end of beta > 0
    enddo



    ! Allocate T-Matrix
    allocate(this%Tnm(this%no_out_modes,this%no_out_modes))
    allocate(this%Rnm(this%no_out_modes,this%no_out_modes))
    this%Tnm      = 0
    this%Rnm      = 0

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
                current = (mode_current_sparse(c(i,:),lambda,sparse_tau_vals,sparse_tau_rcid,sparse_tau_novals))
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
                    this%currents(M_OUT,this%no_out_modes+no_e_out)  = abs(lambda)
            else if( abs(lambda) < 1.0 .and. abs(lambda)> 1.D-16 ) then
                    no_e_in = no_e_in + 1
                    this%modes  (M_IN,this%no_in_modes+no_e_in,:) = c(i,:)
                    this%lambdas(M_IN,this%no_in_modes+no_e_in)   = lambda
                    this%currents(M_IN,this%no_in_modes+no_e_in)  = abs(lambda)
            endif
        endif
    enddo


    ! Sorting propagating modes by the current amplitde
    ! and diagonalizing current operator for degenerated states
    call sort_modes(this%no_in_modes ,this%modes(M_IN,:,:) ,this%lambdas(M_IN,:),this%currents(M_IN,:))
    lambda = this%lambdas(M_IN,1)
    j      = 1
    do i = 2 ,this%no_in_modes
        if( abs(lambda -this%lambdas(M_IN,i)) > 1.0D-10  ) then
            call this%diagonalize_currents(M_IN,j,i-1)
            j = i
            lambda = this%lambdas(M_IN,i)
        endif
    enddo
    call this%diagonalize_currents(M_IN,j,i-1)

    ! Sorting propagating modes by the current amplitde
    ! and diagonalizing current operator for degenerated states
    call sort_modes(this%no_out_modes,this%modes(M_OUT,:,:),this%lambdas(M_OUT,:),this%currents(M_OUT,:))


    lambda = this%lambdas(M_OUT,1)
    j      = 1
    do i = 2 ,this%no_out_modes
        if( abs(lambda - this%lambdas(M_OUT,i)) > 1.0D-10  ) then
            call this%diagonalize_currents(M_OUT,j,i-1)
            j = i
            lambda = this%lambdas(M_OUT,i)
        endif
    enddo
    call this%diagonalize_currents(M_OUT,j,i-1)


    call sort_modes(this%no_out_em,&
        this%modes   (M_OUT,this%no_out_modes+1:this%no_out_modes+this%no_out_em,:),&
        this%lambdas (M_OUT,this%no_out_modes+1:this%no_out_modes+this%no_out_em),  &
        this%currents(M_OUT,this%no_out_modes+1:this%no_out_modes+this%no_out_em), bInverse=.false.)

    ! Calculate overlap matrix of size (no_modes)x(no_modes)
    ! note that in systems like graphene we avoid overlaps of
    ! states with |lambda| = +inf.
    no_modes = this%no_out_modes + QSYS_SCATTERING_QTBM_NO_EVAN
    if(QSYS_SCATTERING_QTBM_NO_EVAN == QSYS_SCATTERING_QTBM_TAKE_ALL_EVAN) &
                no_modes = this%no_out_modes + this%no_out_em


    deallocate(MA,MB,this%UTildeDagger)
    allocate  (MA (no_modes,no_modes))
    allocate  (MB (no_modes,N))
    allocate  (this%UTildeDagger(no_modes,no_modes))

    MB   = conjg((this%modes(M_OUT,1:no_modes,:)))
    ! Calculate overlap matrix and its inverse
    ! Skp = < Uk- | Up- >
    call ZGEMM( 'N', 'T', no_modes,no_modes,N,one, &
                MB, no_modes, this%modes(M_OUT,:,:), N, zero,this%UTildeDagger,no_modes )
    MA  = this%UTildeDagger
    print*,"SYS::QTBM overlap matrix size:",no_modes," cond(S)=",cond_SVD(no_modes,MA)
    if(no_modes > 0) call inverse_matrix(no_modes,this%UTildeDagger)

    MA   = 0
    do i = 1 , no_modes
        Z21(i,:) = this%modes(M_OUT,i,:)
        Z21(i,:) = Z21(i,:)/sqrt(sum(abs(Z21(i,:))**2))
    enddo
    MB   = conjg(Z21)
    blochF(M_OUT,:,:) = 0

    ! Calculate overlap matrix and its inverse
    ! Skp = < Uk- | Up- >
    call ZGEMM( 'N', 'T', no_modes,no_modes,N,one, &
                MB, no_modes, Z21, N, zero,MA,no_modes )
    MB(1:no_modes,1:no_modes) = MA
    if(no_modes > 0) call inverse_matrix(no_modes,MA)

    deallocate(MB)
    allocate(MB (N,N))
    MB = 0

    ! Premultiply Inv(S) matrix by DGamma_k=(lambda_k - lambda_k^-1) :
    ! Tilde(S) = Inv(S)*DGamma
    ! Where DGamma is diagonal matrix, this can be done faster without blas
    do k = 1 , no_modes
    do p = 1 , no_modes
        MA(k,p) = MA(k,p)*(this%lambdas(M_OUT,k)-this%lambdas(M_OUT,k)**(-1))
    enddo
    enddo
    ! Calculate two matrix multiplication with blas
    ! Q = U * Tilde(S)
    MB = 0
    call ZGEMM( 'T', 'N', N,no_modes,no_modes,one,Z21,N, MA, no_modes, zero,MB,N)
    ! F = Q * U^dagger
    Z11 = conjg(Z21) ! We don't need to transpose because matrix U
                                       ! is already stored in that form
    call ZGEMM( 'N', 'N', N,N,no_modes,one,MB,N,Z11,N,zero,blochF(M_OUT,:,:),N)

    ! Calculating of SigmaMatrix
    do i = 1 , N
    do j = 1 , N
        Mdiag(i,j) = conjg(this%valsTau(j,i)) - Ef * conjg(this%valsS1(j,i)) ! Dag of Tau
    enddo
    enddo
    ! Perform multiplication with blas
    call ZGEMM('N','N',N,N,N,one,Mdiag, N, blochF(M_OUT,:,:), N,zero,this%SigmaMat, N )
    ! add to sigma H0 internal hamiltonian
    this%LambdaMat =  this%SigmaMat
    this%SigmaMat  = -this%SigmaMat + this%valsH0


    ! ----------------------------------------------------------------------------
    ! Other possibilities ?
    ! ----------------------------------------------------------------------------
    else ! else of scattering method
        print*,"SYS::ERROR::Scattering problem method is not define: QSYS_SCATTERING_METHOD=",QSYS_SCATTERING_METHOD
        print*,"            QSYS_SCATTERING_METHOD=",QSYS_SCATTERING_METHOD
        stop
    endif

    if(QSYS_DEBUG_LEVEL > 0) then
        print*,"SYS::LEAD::Lead stats:"
        print*,"           No. incoming modes:",this%no_in_modes
        print*,"           No. outgoing modes:",this%no_out_modes
        print*,"           No. in. evan.modes:",this%no_in_em
        print*,"           No. out.evan.modes:",this%no_out_em
        print*,"           No. inf modes     :",no_inf_modes
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

    if(QSYS_DEBUG_LEVEL > 0) then
    print*,"-------------------------------------------------------"
    print*," K vec.  :    In     |       Out  |          Fluxes   "
    print*,"-------------------------------------------------------"
    do i = 1 , this%no_in_modes
        dval1 = log(this%lambdas(M_IN ,i))/II
        dval2 = log(this%lambdas(M_OUT,i))/II
        print"(A,i4,A,f10.4,A,1f10.4,A,2f12.8)","   K[",i,"]:",dval1," | ",dval2 ," | " , this%currents(:,i)
    enddo
    print*,"-------------------------------------------------------"
    endif



    if(allocated(iselect))        deallocate(iselect)
    if(allocated(ALPHA))          deallocate(ALPHA)
    if(allocated(BETA))           deallocate(BETA)
    if(allocated(RWORK))          deallocate(RWORK)
    if(allocated(WORK))           deallocate(WORK)
    if(allocated(MA))             deallocate(MA)
    if(allocated(MB))             deallocate(MB)
    if(allocated(Z21))            deallocate(Z21)
    if(allocated(Z11))            deallocate(Z11)
    if(allocated(Z))              deallocate(Z)
    if(allocated(d))              deallocate(d)
    if(allocated(c))              deallocate(c)
    if(allocated(Mdiag))          deallocate(Mdiag)
    if(allocated(blochF))         deallocate(blochF)
    if(allocated(sparse_tau_vals))deallocate(sparse_tau_vals)
    if(allocated(sparse_tau_rcid))deallocate(sparse_tau_rcid)


    if(QSYS_DEBUG_LEVEL > 0) then
    print*,"SYS::LEAD::Modes calculated in time:", get_clock() - time , "[s]"
    endif
end subroutine calculate_modes

! ------------------------------------------------------------------------
! Calculate tranimission amplitude in the output lead
! all_atoms - reference to all atoms in the system
! phi       - calculated scattering wave function
! The result is holded: this%Tnm(1,1) element. For each incoming mode
! this number should be <= 1
! ------------------------------------------------------------------------
subroutine calculate_Tnm(this,all_atoms,phi,inputmode)
    class(qlead)              :: this
    type(qatom),dimension(:)  :: all_atoms
    complex*16 ,dimension(:)  :: phi
    integer,optional :: inputmode
    ! local variables
    complex*16 ,allocatable , dimension(:) :: leadPhi,tmpVec
    complex*16 :: tmpT
    logical :: bInputLead
    integer :: i,la,ls,lg,no_modes

    bInputLead = .false.
    if(present(inputmode)) bInputLead = .true.

    if(this%lead_type == QSYS_LEAD_TYPE_PSEUDO_TRANSPARENT) then
        this%Tnm = 0
        return
    endif

    allocate(leadPhi(this%no_sites))
    do i = 1 , this%no_sites
        la = this%l2g(i,1)
        ls = this%l2g(i,2)
        lg = all_atoms(la)%globalIDs(ls)
        if(bInputLead) then
            leadPhi(i) = phi(lg) - this%modes(M_IN,inputmode,i)
        else
            leadPhi(i) = phi(lg)
        endif
    enddo



    if(QSYS_SCATTERING_METHOD == QSYS_SCATTERING_WFM) then
        ! Use standard equation when Schur vector are not used
        if(this%bUseShurDecomposition) then
            no_modes = size(this%UTildeDagger,1)
            allocate(tmpVec(no_modes))

            do i = 1 , no_modes
                tmpVec(i) = sum(conjg(this%modes(M_OUT,i,:))*leadPhi)
            enddo
            do i = 1 , this%no_out_modes
                tmpT = sum(this%UTildeDagger(i,:) * tmpVec(:))
                this%Tnm(i,1) =  abs(tmpT)**2
            enddo
            deallocate(tmpVec)
        else
            do i = 1 , this%no_out_modes
                tmpT = sum(this%UTildeDagger(:,i) * leadPhi(:))
                this%Tnm(i,1) =  abs(tmpT)**2
            enddo
        endif
    else if(QSYS_SCATTERING_METHOD == QSYS_SCATTERING_QTBM) then

        no_modes = size(this%UTildeDagger,1)
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

    this%modeT = 0
    if(this%no_out_modes > 0) this%modeT = sum(this%Tnm(:,1))
    deallocate(leadPhi)
endsubroutine calculate_Tnm

! ------------------------------------------------------------------------
! Calculate flux current for a givem mode c with lambda using the fact
! that the coupling matrix is usually sparse. Performs sparse multiplication.
! Returns the value of the current carried by this mode.
! The current is calculated with general formula:
! j =  - 2 * Imag( lambda c(:)^dag * Tau * c(:) )
! ------------------------------------------------------------------------
doubleprecision function mode_current_sparse(c,lambda,svals,rowcols,nvals)

    complex*16 , dimension(:)   :: c
    complex*16                  :: lambda
    complex*16 , dimension(:)   :: svals
    integer    , dimension(:,:) :: rowcols
    integer                     :: nvals
    integer :: k,i,j
    complex*16 :: current
    current=0
    do k = 1 , nvals
        i = rowcols(k,1)
        j = rowcols(k,2)
        current = current + conjg(c(i)) * svals(k) * c(j)
    enddo
    mode_current_sparse = -2 * IMAG(lambda*current)
end function mode_current_sparse

! ------------------------------------------------------------------------
! Sort propagating modes in ascending order of the current amplitude
! ------------------------------------------------------------------------
subroutine sort_modes(N,vectors,lambdas,currents,bInverse)
    complex*16,dimension(:,:) :: vectors
    complex*16,dimension(:)   :: lambdas
    doubleprecision,dimension(:)   :: currents
    logical,optional :: bInverse
    integer :: N
    complex*16,dimension(:),allocatable  :: tmpvec
    complex*16 :: tmpval
    doubleprecision :: dval , mindval , tmpcurr
    integer :: i , j  , nvec , imin
    logical :: bInv

    nvec = size(vectors(1,:))
    allocate(tmpvec(nvec))

    bInv = .false.
    if(present(bInverse)) then
    if(bInverse == .true.) bInv = .true.
    endif


    do i = 1 , n
        imin   = i
        if(bInv) then ! sort order
            do j = i+1 , n
               dval     = abs(currents(j))
               mindval  = abs(currents(imin))
               if( dval < mindval ) imin = j
            enddo
        else
            do j = i+1 , n
               dval     = abs(currents(j))
               mindval  = abs(currents(imin))
               if( dval > mindval ) imin = j
            enddo
        endif

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

! ------------------------------------------------------------------------
! Save lead data for further debuging with the Viewer.
! ------------------------------------------------------------------------
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
! -------------------------------------------------
! Diagonalize current operator for degenerate states
! for more details see Micheal Wimmer book.
! -------------------------------------------------
subroutine diagonalize_currents(this,mdir,ifrom,ito)
    class(qlead) :: this

    integer :: mdir,ifrom,ito,N,i,j,Nlead

    complex*16,allocatable :: PhiVecs(:,:),PhiVec(:), VelMat(:,:) , BaseMat(:,:) , tmpMat(:,:) , tmpVec(:)
    complex*16 :: alpha,beta,lambda

    doubleprecision :: current

    alpha = II
    beta  = 0.0

    N     = ito-ifrom+1

    ! There is no need to diagonalize the non degenerate states
    if(N == 1) return

    Nlead = this%no_sites

    allocate(VelMat(N,N))
    allocate(BaseMat(N,N))
    allocate(tmpMat(Nlead,Nlead))
    allocate(tmpVec(Nlead))

    allocate(PhiVec(Nlead))

    tmpMat = 0
    lambda = this%lambdas(mdir,ifrom)
    do i = 1 , Nlead
    do j = 1 , Nlead
        tmpMat(i,j) = + ( this%valsTau(i,j)*lambda &
                  - conjg(this%valsTau(j,i)*lambda))
    enddo
    enddo

    do i = 0 , N-1
    do j = 0 , N-1
       PhiVec = this%modes(mdir,ifrom+j,:)
       call ZGEMV ( 'N', Nlead, Nlead, ALPHA, tmpMat, Nlead,PhiVec, 1, BETA, tmpVec, 1 )
       VelMat(i+1,j+1) = sum(conjg(this%modes(mdir,ifrom+i,:))*tmpVec)
    enddo
    enddo

    call alg_ZHEEV(N,VelMat,BaseMat,tmpVec(1:N))

    do i = 1 , Nlead
        beta = 0
    do j = 1 , N
        tmpVec(j) = sum(this%modes(mdir,ifrom:ito,i)*BaseMat(:,j))
    enddo
        this%modes(mdir,ifrom:ito,i) = tmpVec(1:N)
    enddo
    ! normalize currents
    do i = ifrom , ito
        current = (mode_current_sparse(this%modes(mdir,i,:),lambda,sparse_tau_vals,sparse_tau_rcid,sparse_tau_novals))
        this%modes(mdir,i,:)  = this%modes(mdir,i,:)/sqrt(abs(current))
        this%currents(mdir,i) =-1.0
        if(current > 0) this%currents(mdir,i) = 1.0
    enddo


    deallocate(VelMat)
    deallocate(BaseMat)
    deallocate(tmpMat)
    deallocate(tmpVec)
    deallocate(PhiVec)

end subroutine diagonalize_currents



! ------------------------------------------------------------------------
! Calculate average spins for incoming modes. This function can be
! used for modes with two orbitals per atom, other wise it will not work,
! since Paulli matrices are involved here.
! dir - M_IN=1 or M_OUT=2 - determine the input or output propagating modes
! spins(no_incoming_modes,{1-x,2-y,3-z}) returns the values of spin
!       polarization for {x,y,z} directions.
! ------------------------------------------------------------------------
subroutine calc_average_spins(this,dir,spins)
    class(qlead) :: this
    integer :: dir
    doubleprecision,allocatable :: spins(:,:)
    ! local variables
    complex*16,allocatable :: Chi_up(:),Chi_down(:)
    complex*16 :: YA,YB
    integer :: no_modes,i,m,spin,iter1,iter2

    print*,"SYS::LEAD::Calculation of the average spin polarization"

    if(dir == M_IN ) no_modes = this%no_in_modes
    if(dir == M_OUT) no_modes = this%no_out_modes



    if(allocated(spins)) deallocate(spins)
    allocate(spins   (no_modes,3)) ! (modes,{x,y,z})
    allocate(Chi_up  (this%no_sites))
    allocate(Chi_down(this%no_sites))

    do m = 1 , no_modes
        Chi_up    = 0
        Chi_down  = 0
        iter1     = 0
        iter2     = 0
        do i = 1, this%no_sites
            spin = this%l2g(i,2)
            ! Filling the spinors
            if(spin == 1) then ! up
                iter1 = iter1 + 1
                Chi_up  (iter1) = this%modes(dir,m,i)
            else if(spin == 2) then ! down
                iter2 = iter2 + 1
                Chi_down(iter2) = this%modes(dir,m,i)
            else if(spin > 2) then ! not supported
                print*,"SYS::LEAD::Calc_average_spin works only with the systems with the s=1/2."
                print*,"           So it cannot be used for systems with s>1/2. Finded spin states:",spin
                stop -1
            endif
        enddo ! end of loop over sites in the lead


            ! Direction Z
            YA = sum(abs(Chi_up)**2)
            YB = sum(abs(Chi_down)**2)
            spins(m,3) = (YA-YB)/(YA+YB) ! kierunek z

            ! Direction X
            YA = sum(abs(Chi_up+Chi_down)**2) ! up
            YB = sum(abs(Chi_up-Chi_down)**2) ! down
            spins(m,1) = (YA-YB)/(YA+YB) ! kierunek x

            ! Direction Y
            YA = sum(abs(Chi_up+II*Chi_down)**2)
            YB = sum(abs(Chi_up-II*Chi_down)**2)
            spins(m,2) = (YA-YB)/(YA+YB) ! kierunek y

    enddo



    deallocate(Chi_up,Chi_down)
end subroutine calc_average_spins

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

subroutine alg_ZHEEV(N,Mat,Vecs,Lambdas)
    integer :: N
    complex*16 :: Mat(:,:),Vecs(:,:),Lambdas(:)

    ! -------------------------------------------------
    !                    LAPACK
    ! -------------------------------------------------
    !     .. Local Scalars ..
    INTEGER          INFO, LWORK, LWMAX

    !     .. Local Arrays ..

    DOUBLE PRECISION,allocatable :: RWORK( : )
    COMPLEX*16,allocatable        :: WORK( : )


    ! Initialize lapack and allocate arrays
    LWMAX = N*50

    allocate(RWORK ( 3*N  ))
    allocate(WORK  ( LWMAX ))

    !
    !     Query the optimal workspace.
    !
    LWORK  = -1
    call zheev("N", "L", N, Mat, N, Lambdas, work, lwork, rwork, info)


    if( INFO /= 0 ) then
        print*,"  alg_ZHEEV: Error during solving with info:",INFO
        stop
    endif
    LWORK  = MIN( LWMAX, INT( WORK( 1 ) ) )

    deallocate( WORK)
    allocate(WORK  ( LWORK ))

    Vecs = Mat
    call zheev("V", "L", N, Vecs, N, Lambdas, work, lwork, rwork, info)

    deallocate(RWORK)
    deallocate(WORK)

end subroutine alg_ZHEEV

doubleprecision function cond_number(N,A) result(rval)
  integer :: N
  complex*16,dimension(:,:):: A
  complex*16,allocatable,dimension(:)       :: WORK
  doubleprecision,allocatable,dimension(:)  :: RWORK
  integer,allocatable,dimension (:)    :: IPIV
  complex*16,dimension(:,:),allocatable:: tmpA
  integer info,error
  doubleprecision :: anorm,norm
  external zlange
  doubleprecision zlange
  B_SINGULAR_MATRIX = .false.

  allocate(WORK(2*N),IPIV(N),RWORK(2*N),tmpA(N,N),stat=error)
  if (error.ne.0)then
    print *,"SYS::LEAD::ZGETRF::error:not enough memory"
    stop
  end if

  tmpA = A

  anorm = ZLANGE( '1', N, N, tmpA, N, WORK )

  call ZGETRF(N,N,tmpA,N,IPIV,info)
  call zgecon( '1', N, tmpA, N, anorm, norm, work, rwork, info )
  rval = norm


  deallocate(IPIV,WORK,RWORK,tmpA,stat=error)
  if (error.ne.0)then
    print *,"SYS::LEAD::ZGETRF::error:fail to release"
    stop
  end if
end function cond_number

doubleprecision function cond_SVD(N,A) result(rval)
    integer :: N
    complex*16,dimension(:,:):: A
    complex*16 , allocatable :: tmpA(:,:) , VT(:,:),tU(:,:),tVT(:,:)
    doubleprecision , allocatable :: S(:) , tS(:)
    integer :: i,j,M
    allocate(tmpA(N,N))
    tmpA = A
    call ZSVD(N,tmpA,tU,tS,tVT)
    M = 0
    do i = 1 , N
    if(abs(tS(i)) > 1.0d-20 ) M = M+1
    enddo
    rval = abs(tS(1)/tS(N))
    deallocate(tmpA,tU,tS,tVT)
end function cond_SVD


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
      if(allocated(U)) deallocate(U)
      if(allocated(S)) deallocate(S)
      if(allocated(VT)) deallocate(VT)
      allocate(U( LDU, N ), VT( LDVT, N ), S( N ),WORK(LWMAX) , RWORK(5*N))


      CALL ZGESVD( 'All', 'All', N, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK , RWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

      CALL ZGESVD( 'All', 'All', N, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK , RWORK, INFO )

      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm computing SVD failed to converge.'
         STOP
      END IF


end subroutine ZSVD

subroutine inverse_svd(N,A)
    complex*16 :: A(:,:)
    integer :: N,i,j
    complex*16 ,allocatable , dimension(:,:) :: U,VT
    doubleprecision,allocatable,dimension(:) :: S


    call ZSVD(N,A,U,S,VT)

    do i = 1, N
    do j = 1, N
        A(i,j) = sum( conjg(Vt(:,i))*(S(:)**(-1))*conjg(U(j,:)) )
    enddo
    enddo

    deallocate(U,S,VT)
end subroutine inverse_svd


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

    eps        = 1.0D-16
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

endmodule modlead
