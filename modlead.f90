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
use modalgs

implicit none
private

ENUM, BIND(C)
  ENUMERATOR :: QLEAD_MS_TO_STANDARD   = 1 ! find eigen modes by converting from GEVP to SEVP, if Tau is invertible
  ENUMERATOR :: QLEAD_MS_GENERALIZED   = 2 ! Use generalized eigenvalue problem to find lead modes
END ENUM

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
    complex*16,allocatable      :: QLdcmp_Qmat(:,:),QLdcmp_Lmat(:,:) ! Finding modes transmission with QL decomposition


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
!    procedure,pass(this) :: extract_modes_wfm
endtype qlead

public :: qlead
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
    if(QSYS_DEBUG_LEVEL>0) print"(A,3f12.4,A)"," SYS::LEAD::Lead translation vector XYZ=(",cellBA_vec,")"

    ! -----------------------------------------------------------
    ! Creating the Hamiltonian of the Lead and Coupling matrix
    ! -----------------------------------------------------------
    do i = 1 , no_sites
        atom_id = this%l2g(i,1) ! get atom index , this%l2g(i,2) tells what is the spinID of that atom
        s1      = this%l2g(i,2) ! get atom spinID index
        ! iterate over all bonds in that atom (A)
        do b = 1 , all_atoms(atom_id)%no_bonds
            ns1 = size(all_atoms(atom_id)%bonds(b)%bondMatrix,1)
            ns2 = size(all_atoms(atom_id)%bonds(b)%bondMatrix,2)

            ! get ids of associated atom
            bond_atom_id = all_atoms(atom_id)%bonds(b)%toAtomID
            bond_id      = all_atoms(bond_atom_id)%globalIDs(1)
            ! if tmp_g2l(bond_id) == 0 it means that this atom with id = bond_atom_id
            ! is not in the lead, so it has to be located in the next unit cell.
            if( tmp_g2l(bond_id) == 0 ) then

                ! now we search for an atom in the main uint cell with the same position
                ! what assiociated atom in the next lead.
                cellB_pos = all_atoms(bond_atom_id)%atom_pos
                do j = 1 , no_sites
                    cellA_pos = all_atoms(this%l2g(j,1))%atom_pos + cellBA_vec
                    if( sqrt(sum( (cellB_pos - cellA_pos )**2)) < minimum_distance*1.0D-5 ) then
                        exit
                    endif
                enddo
                if( j > no_sites ) then
                    ! check if atom with id bond_atom_id is in the next unit cell not in the
                    ! previous one.
                    if(QSYS_DEBUG_LEVEL<2) then
                        if( lshape%is_inside(cellB_pos+cellBA_vec) == .true. ) cycle
                    endif
                    print"(A)",              " SYS::LEAD::ERROR::The translation"
                    print"(A,i9,A,3e12.4,A)","                   from atom:",atom_id," at position r=(",all_atoms(atom_id)%atom_pos,")"
                    print"(A,i9,A,3e12.4,A)","                   to atom  :",bond_atom_id," at position r=(",all_atoms(bond_atom_id)%atom_pos,")"
                    print"(A)",              "                   is impossible by applying lead translation vector. Maybe some of the "
                    print"(A)",              "                   atoms have not been inlcuded in the lead. Check lead area and T vector."
                    cycle
                endif
                ! now "j" contains of localID of the state in the lead which is equivalent
                ! to same state but in the next unit cell.
                do s2 = 1 , ns2
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
                enddo
            else ! in case of coupling which occures within lead
                do s2 = 1 , ns2
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
                enddo
            endif ! end of else if "atom is in the lead or outside"
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

    if(allocated(this%QLdcmp_Qmat))  deallocate(this%QLdcmp_Qmat)
    if(allocated(this%QLdcmp_Lmat))  deallocate(this%QLdcmp_Lmat)

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
    DOUBLE PRECISION,allocatable,dimension(:)     :: RWORK( : )

    integer :: k,p,i,j,no_in,no_out,no_e_in,no_e_out,no_modes,no_inf_modes

    INTEGER                                      :: INFO
    COMPLEX*16 , dimension(:) ,     allocatable  :: ALPHA , BETA , QLdcmp_TauVec

    logical, dimension(:), allocatable    :: iselect

    doubleprecision :: tmpc,current,time,dval1,dval2
    COMPLEX*16 ::lambda,one,zero
    logical :: bShurDecompositionForStandEP

    integer :: STABILIZATION_METHOD = QLEAD_MS_TO_STANDARD

    ! Profiling - used to debug program SPEED, with QSYS_DEBUG_LEVEL == 2
    doubleprecision :: T_total
    doubleprecision :: T_solving
    doubleprecision :: T_ordering
    doubleprecision :: T_SE_constr
    doubleprecision :: T_BLOCH_constr

    T_total = get_clock();

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

    if(allocated(this%QLdcmp_Qmat))  deallocate(this%QLdcmp_Qmat)
    if(allocated(this%QLdcmp_Lmat))  deallocate(this%QLdcmp_Lmat)

    allocate(this%modes     (2,N,N))
    allocate(this%lambdas   (2,N))
    allocate(this%SigmaMat  (N,N))
    allocate(this%LambdaMat (N,N))
    allocate(this%currents  (2,N))
    allocate(this%UTildeDagger(N,N))
    allocate(this%QLdcmp_Qmat(N,N))
    allocate(this%QLdcmp_Lmat(N,N))
    allocate(QLdcmp_TauVec(N))


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



    allocate(ALPHA  (2*N))
    allocate(BETA   (2*N))
    allocate(RWORK  (2*N))

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

    T_solving = get_clock()
    if(QSYS_DEBUG_LEVEL>0) then
        print*,"SYS::LEAD::Solving eigenvalue problem for leads."
    endif
    ! ---------------------------------------------------------------------------
    ! If the condition number of coupling matrix is small enough
    ! we can try to convert Generalized eigen value problem to
    ! standard one.
    ! ---------------------------------------------------------------------------
    Z11 = this%valsTau  - this%valsS1*Ef
    CONDA = alg_cond(Z11)
    if(QSYS_FORCE_ZGGEV_TO_FIND_MODES) CONDA = 1.0D16 ! Go to GGEV mehtod if forced
    if( CONDA * qsys_double_error < QSYS_ERROR_EPS ) then

        if(QSYS_DEBUG_LEVEL > 0) then
            print*,"SYS::LEAD::Tau is not ill-conditioned thus solving GEEV."
        endif

        STABILIZATION_METHOD = QLEAD_MS_TO_STANDARD
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


        one  = -1
        zero =  0
        ! Calculate new matrix block element
        ! TauDag^-1 * Tau

        call ZGEMM( 'N', 'N', N, N, N, one ,Z11 , N , &
                            Mdiag,  N, zero,Z21 , N )
        MA(N+1:2*N,1:N) = Z21
        ! filling diag with diagonal matrix
        Mdiag = 0
        do i = 1 , N
            Mdiag(i,i) =  1
        enddo
        MA(1:N,N+1:2*N) =  Mdiag

        ! Calculate new matrix block element
        ! TauDag^-1 * (E - H0)
        one = +1
        Z21 = this%valsS0*Ef - this%valsH0

        call ZGEMM( 'N', 'N', N, N, N, one ,Z11  , N, &
                            Z21, N   , zero,Mdiag, N )
        MA(N+1:2*N,N+1:2*N) = Mdiag

        ! Solve eigen problem normally or use Schur decomposition
        ! Use  Schur decomposition to solve eigen problem
        if(QSYS_FORCE_SCHUR_DECOMPOSITION) then
            if(QSYS_DEBUG_LEVEL > 0) then
                print*,"SYS::LEAD::Using Schur decomposition to find modes."
            endif
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
        else ! Use the simplest approach
            if(QSYS_DEBUG_LEVEL > 0) then
                print*,"SYS::LEAD::Using standard GEEV method to find modes."
            endif
            ! Solve standard Eigen value problem
            call geev(MA, ALPHA, vr=Z ,info=INFO)
            ! Checking solution
            if( INFO /= 0 ) then
                print*,"SYS::LEAD::Cannot solve eigenvalue problem for eigenmodes: ZGEEV info:",INFO
                stop
            endif
        endif ! end if QSYS_FORCE_SCHUR_DECOMPOSITION
        BETA  = 1.0
    ! ------------------------------------------------------------------------------
    else ! Solve problem with Generalized eigen value solver.
    ! ------------------------------------------------------------------------------
        if(QSYS_DEBUG_LEVEL > 0) then
            print*,"SYS::LEAD::Solving GGEV problem: Tau is ill-conditioned?"
        endif
    ! -------------------------------------------------------
    ! Creation of the Generalized eigenvalue problem Eq. (52)
    ! Use this method when coupling matrix is singular.
    ! -------------------------------------------------------
        STABILIZATION_METHOD = QLEAD_MS_GENERALIZED

        Z11 = this%valsS0*Ef - this%valsH0
        Z21 = this%valsTau   - this%valsS1*Ef
        ! Try to convert Generalized eigen value problem to standard one
        ! using the SVD, if this method will lead to possible unstable
        ! solutions global variable B_SINGULAR_MATRIX will be equal false.
        B_SINGULAR_MATRIX = .true. ! if schur factorization if forced skip try_svd_modes_decomposition
        if(.not. QSYS_FORCE_SCHUR_DECOMPOSITION &
        .and. .not. QSYS_FORCE_ZGGEV_TO_FIND_MODES ) call try_svd_modes_decomposition(Z11,Z21,ALPHA,BETA,Z)
        ! Solve directly generalized eigen value problem without tricks
        if(B_SINGULAR_MATRIX) then
            ! this%valsTau contains TauDagger coupling matrix as in Zwierzycki paper
            do i = 1 , N
            do j = 1 , N
                Mdiag(i,j) =  conjg(this%valsTau(j,i)) - Ef*conjg(this%valsS1(j,i)) ! Dag of Tau
            enddo
            enddo

            ! Filling matrices
            ! MA=  (   0       1   )    MB = ( 1  0 )
            !      (  -t^*   E-H   )         ( 0  t )

            MA  = 0
            MB  = 0
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

            ! ---------------------------------------------------------------
            ! Solve problem using Shur decompostion
            ! ---------------------------------------------------------------
            if(QSYS_FORCE_SCHUR_DECOMPOSITION) then
                allocate(Qshur      (2*N,2*N))
                allocate(Zshur      (2*N,2*N))
                allocate(Sshur      (2*N,2*N))
                allocate(Pshur      (2*N,2*N))

                Sshur = MA
                Pshur = MB
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
            else ! end of else QSYS_FORCE_SCHUR_DECOMPOSITION: using GGEV instead
                call GGEV(MA, MB, alpha, beta , vr = Z,info=info)
                if( INFO /= 0 ) then
                    print*,"SYS::LEAD::Cannot solve generalized eigenvalue problem for eigenmodes: ZGGEV info:",INFO
                    stop
                endif

            endif ! end of use QSYS_FORCE_SCHUR_DECOMPOSITION

        endif ! end of if B_SINGULAR_MATRIX after try_svd_modes_decomposition call

    endif ! end of if try to convert to standard eigenvalue problem
    ! ------------------------------------------------------------------------------
     T_solving = get_clock() - T_solving

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
    T_ordering = get_clock()
    call extract_modes_wfm(this,N,c,d,ALPHA,BETA,Z,Z11,Z21,blochF,iselect,no_inf_modes,no_in,no_out,no_e_in,no_e_out)
    T_ordering = get_clock() - T_ordering

    T_BLOCH_constr = get_clock()
    blochF(M_OUT,:,:)  = 0
    select case(STABILIZATION_METHOD)
        ! -----------------------------------------------------------
        case (QLEAD_MS_TO_STANDARD,QLEAD_MS_GENERALIZED)
        ! -----------------------------------------------------------
        ! Construction of the F^-1+ matrix Eq. (57)
            if(QSYS_FORCE_SCHUR_DECOMPOSITION) then
                Z  = Zshur
                ! Order Schur matrices with iselect matrix
                if(STABILIZATION_METHOD == QLEAD_MS_TO_STANDARD) then
                    call trsen(Sshur, select=iselect,q=Z ,info=info)
                else ! if STABILIZATION_METHOD == QLEAD_MS_GENERALIZED
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
                if(B_SINGULAR_MATRIX) then
                    print*,"QSYS::LEAD::Shur Z21 matrix is non-invertible, try to use QTBM to solve your problem"
                    print*,"            by setting:QSYS_SCATTERING_METHOD = QSYS_SCATTERING_QTBM in main.f90"
                    stop
                endif
                ! Calculate F as F = Z11 * Z21^(-1)
                call gemm(Z11, Z21, Mdiag)
                blochF(M_OUT,:,:) = Mdiag
                this%bUseShurDecomposition =.true.

            ! ------------------------------------------------------
            else ! use Schur decomposition
            ! ------------------------------------------------------
            ! Check if matrix Z21 = \lambda*U(-) is well conditioned
            ! otherwise use SVD to invert matrix to improve statilibty

            if(alg_cond(Z21)*qsys_double_error < QSYS_ERROR_EPS) then
                ! Another version of Z21 inversion :D
!                this%bUseShurDecomposition =.true.
!                Mdiag(:,:) = Z21
!                call inverse_matrix(N,Mdiag)
!                blochF(M_OUT,:,:)  = 0
!                one  = 1.0
!                zero = 0.0
!                ! Use blas for matrix multiplication
!                call ZGEMM('N','N',N,N,N,one,Z11,N,Mdiag,N,zero,Z21,N)
!                blochF(M_OUT,:,:) = Z21

                Mdiag(:,:) = this%modes(M_OUT,:,:)
                call inverse_matrix(N,Mdiag)
!                this%UTildeDagger  = Mdiag
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

            else ! use SVD to calculate Bloch matrix since Z21 is ill-conditioned
                call calc_bloch_matrix(Z11,Z21,Mdiag) ! using SVD
                blochF(M_OUT,:,:) = Mdiag
                this%bUseShurDecomposition =.true.
            endif ! if matrix \lambda*U(-) is well conditioned

            endif ! if QSYS_FORCE_SCHUR_DECOMPOSITION

!        ! -----------------------------------------------------------
!        case (QLEAD_MS_GENERALIZED)
!        ! -----------------------------------------------------------
        case default
            print*,"QSYS::QLEAD::There is no such method of modes stabilization: ",STABILIZATION_METHOD
            stop
    end select

    T_BLOCH_constr = get_clock() - T_BLOCH_constr

    if(QSYS_DEBUG_LEVEL > 0) then
        print*,"SYS::LEAD::Calculation of Self Energies."
    endif
    T_SE_constr = get_clock()
    ! Calculating of SigmaMatrix - SELF ENERGY
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

    no_modes = this%no_out_modes + this%no_out_em - no_inf_modes
    deallocate(RWORK)
    allocate(RWORK(this%no_sites))

    ! Calculate the SVD decomposition of modes matrix U(modes)
    ! U(modes) = U_svd * S * V^+ we will use U_svd matrix to
    ! order the eigen values of this matrix

    Mdiag(:,:) = transpose(this%modes(M_OUT,:,:))
    call gesvd(Mdiag,RWORK ,u=Z11,job="All" ,info=info)
    if(info /= 0 ) then
        print*,"LEAD::SVD decomposition error with info=",info
    endif
    this%QLdcmp_Lmat = Z11

    ! Performing QL factorization of the matrix U in order to get more stable
    ! transmision calculation of transmission matrix
    ! Transpose mode matrix to have modes ordered in columns instead of rows
    ! The perform QL factorizaion
    ! We perform QL factorization on transformed U(modes) matrix
    ! U' = U_svd^+ * U(modes)
    Mdiag(:,:) = transpose(this%modes(M_OUT,:,:))
    call gemm(Z11,Mdiag,Z21,transa='C')
    Mdiag            = Z21
    this%QLdcmp_Qmat = Z21

    call geqlf(Mdiag , tau=QLdcmp_TauVec ,info=INFO)
    if(info /= 0 ) then
        print*,"QSYS::LEAD::QR geqlf decomposition error for matrix U(M_OUT,:,:) with info=",INFO
    endif
    ! Get Q matrix.
    Z11 = Mdiag
    call ungql(Z11, QLdcmp_TauVec ,info=INFO)
    if(info /= 0 ) then
        print*,"QSYS::LEAD::QR ungql decomposition error for matrix U(M_OUT,:,:) with info=",INFO
    endif
    ! Get L matrix
    Z21 = this%QLdcmp_Qmat
    call unmql(Mdiag, QLdcmp_TauVec, Z21, side="L" ,trans='C' ,info=INFO)
    if(info /= 0 ) then
        print*,"QSYS::LEAD::QR unmql decomposition error for matrix U(M_OUT,:,:) with info=",INFO
    endif
    ! Back to original space by multiplying obtained matrix Q by U_svd
    ! Q = U_svd * Q'
    call gemm(this%QLdcmp_Lmat,Z11,this%QLdcmp_Qmat)
    this%QLdcmp_Lmat = Z21

    if( QSYS_DEBUG_LEVEL > 1 ) then
        print*,"QSYS::LEAD::Performed QL factorization of the modes matrix "
        print*,"               cond(U)=",alg_cond(this%modes(M_OUT,:,:))
        print*,"               cond(Q)=",alg_cond(this%QLdcmp_Qmat)
        print*,"               cond(L)=",alg_cond(this%QLdcmp_Lmat)
        print*,"        cond(L(modes))=",alg_cond(Z21(1:this%no_out_modes,1:this%no_out_modes))
    endif

    T_SE_constr = get_clock() - T_SE_constr

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
    T_total = get_clock() - T_total

    if(QSYS_DEBUG_LEVEL > 1) then
        print*,"SYS::LEAD::Profiling stats:"
        print*,"           Matrix size                 =",this%no_sites
        print*,"           Eig. problem solve time  [s]=",T_solving
        print*,"           Modes ordering time      [s]=",T_ordering
        print*,"           BLOCH mat. constr. time  [s]=",T_BLOCH_constr
        print*,"           SE mat. constr. time     [s]=",T_SE_constr
        print*,"           TOTAL  time              [s]=",T_total
    endif



    if(allocated(iselect))        deallocate(iselect)
    if(allocated(ALPHA))          deallocate(ALPHA)
    if(allocated(BETA))           deallocate(BETA)
    if(allocated(RWORK))          deallocate(RWORK)

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
    if(allocated(Qshur))          deallocate(Qshur)
    if(allocated(Zshur))          deallocate(Zshur)
    if(allocated(Sshur))          deallocate(Sshur)
    if(allocated(Pshur))          deallocate(Pshur)

    if(QSYS_DEBUG_LEVEL > 0) then
    print*,"SYS::LEAD::Modes calculated in time:", get_clock() - time , "[s]"
    endif
end subroutine calculate_modes


subroutine extract_modes_wfm(this,N,c,d,ALPHA,BETA,Z,Z11,Z21,blochF,iselect,no_inf_modes,no_in,no_out,no_e_in,no_e_out)
    class(qlead)              :: this
    integer :: N
    complex*16, allocatable , dimension(:,:)    :: Z
    complex*16, allocatable , dimension(:,:)    :: d,c,Z11,Z21
    complex*16, allocatable , dimension(:,:,:)  :: blochF
    integer :: k,i,j,no_in,no_out,no_e_in,no_e_out,no_inf_modes
    COMPLEX*16 , dimension(:) ,     allocatable  :: ALPHA , BETA
    logical, dimension(:), allocatable    :: iselect
    doubleprecision                       :: tmpc,current
    COMPLEX*16                            :: lambda

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
!                print*,i,"c=",current,"l=",abs(lambda)
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
    if(allocated(this%Tnm)) deallocate(this%Tnm)
    if(allocated(this%Rnm)) deallocate(this%Rnm)
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
        lambda= (ALPHA(i)/BETA(i))
        if(abs(Beta(i))>1e-16) then
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

            else if( abs(lambda) < 1.0 .and. abs(lambda) > 1.D-16 ) then
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


end subroutine extract_modes_wfm

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
    complex*16 ,allocatable , dimension(:,:) :: leadPhi,tmpVec
    complex*16 ,allocatable , dimension(:) :: tmpTn
    complex*16 :: tmpT
    logical :: bInputLead
    integer :: i,j,la,ls,lg,no_modes

    bInputLead = .false.
    if(present(inputmode)) bInputLead = .true.

    if(this%lead_type == QSYS_LEAD_TYPE_PSEUDO_TRANSPARENT) then
        this%Tnm = 0
        return
    endif

    allocate(leadPhi(this%no_sites,1))

    do i = 1 , this%no_sites
        la = this%l2g(i,1)
        ls = this%l2g(i,2)
        lg = all_atoms(la)%globalIDs(ls)
        if(bInputLead) then
            leadPhi(i,1) = phi(lg) - this%modes(M_IN,inputmode,i)
        else
            leadPhi(i,1) = phi(lg)
        endif
    enddo



    if(QSYS_SCATTERING_METHOD == QSYS_SCATTERING_WFM) then
        ! Use standard equation when Schur vector are not used
!        if(this%bUseShurDecomposition) then
!            no_modes = size(this%UTildeDagger,1)
!            allocate(tmpVec(no_modes))
!
!            do i = 1 , no_modes
!                tmpVec(i) = sum(conjg(this%modes(M_OUT,i,:))*leadPhi)
!            enddo
!            do i = 1 , this%no_out_modes
!                tmpT = sum(this%UTildeDagger(i,:) * tmpVec(:))
!                this%Tnm(i,1) =  abs(tmpT)**2
!            enddo
!            deallocate(tmpVec)
!        else
!            print*,"cond=",this%SVD_S/this%SVD_S(1)
!            print*,"calculation for lead:",bInputLead
!            do i = 1 , this%no_out_modes
!                tmpT = sum(this%UTildeDagger(:,i) * leadPhi(:,1))
!                this%Tnm(i,1) =  abs(tmpT)**2
!                print*,"t",i,"=",abs(tmpT)**2
!            enddo

            no_modes = this%no_out_modes

            allocate(tmpVec (this%no_sites,1))
            allocate(tmpTn  (no_modes))

            tmpVec(:,1) = leadPhi(:,1)
            call gemm(this%QLdcmp_Qmat,leadPhi,tmpVec,transa  = 'C' )

            tmpTn = 0
            do i = 1 , no_modes
                tmpT = 0
                do j = 1 , i-1
                    tmpT = tmpT + tmpTn(j)*this%QLdcmp_Lmat(i,j)
                enddo
                tmpTn(i) =  (tmpVec(i,1) - tmpT)/this%QLdcmp_Lmat(i,i)
                this%Tnm(i,1) =  abs(tmpTn(i))**2
            enddo

            deallocate(tmpVec)
            deallocate(tmpTn)

!            call calculate_Tnm_SVD(this,leadPhi)
!        endif
    else if(QSYS_SCATTERING_METHOD == QSYS_SCATTERING_QTBM) then

        no_modes = size(this%UTildeDagger,1)
        allocate(tmpVec(no_modes,1))
        do i = 1 , no_modes
            tmpVec(i,1) = sum(conjg(this%modes(M_OUT,i,:))*leadPhi(:,1))
        enddo
        do i = 1 , this%no_out_modes
            tmpT = sum(this%UTildeDagger(i,:) * tmpVec(:,1))
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

    complex*16,allocatable :: PhiVec(:), VelMat(:,:) , BaseMat(:,:) , tmpMat(:,:) , tmpVec(:)
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

! -----------------------------------------------------------------------
! Wrapper for stabilized calculation of Bloch matrix, defined as
! F = A * B^(-1)
! B matrix is decomposed with SVD, singular values are removed
! during inversion.
! -----------------------------------------------------------------------
subroutine calc_bloch_matrix(A,B,F)
    complex*16 :: A(:,:),B(:,:),F(:,:)

    complex*16,allocatable      :: svd_U(:,:),svd_Vt(:,:)
    complex*16,allocatable      :: svd_Ur(:,:),svd_Vtr(:,:)
    doubleprecision,allocatable :: svd_S(:)
    integer :: i ,j, n , info , n_svd
    n = size(A,1)

    allocate(svd_U(n,n))
    allocate(svd_Vt(n,n))
    allocate(svd_S(n))

    ! Decompose matrix B = U S V^H
    call gesvd(B,svd_S ,u=svd_U ,vt=svd_Vt ,job="All" ,info=info)
    if(info /= 0 ) then
        print*,"QSYS:::LEAD::SVD decomposition error with info=",info, " in :calc_bloch_matrix"
        stop
    endif
    n_svd = 0
    ! Remove singular values from B matrix, thus reduce the size
    ! of matrices
    do i = 1 , n
!        print*,i,svd_S(i)/svd_S(1)
        if( svd_S(i)/svd_S(1) > QSYS_ERROR_EPS ) then
            n_svd = n_svd + 1
        endif
    enddo

    F = 0
    B = 0
    allocate(svd_Ur(n,n_svd))
    allocate(svd_Vtr(n_svd,n))

    svd_Vtr = svd_Vt(1:n_svd,1:n)

    ! Multiply B' = A * V
    call gemm(A, svd_Vtr, svd_Ur ,transa="N",transb="C" )
!    call gemm(A, svd_Vt(1:n_svd,1:n), B(1:n,1:n_svd) ,transa="N",transb="C" )
    ! Multiply B'' = B' * S (S is diagonal matrix)

    deallocate(svd_Vtr)
    allocate(svd_Vtr(n,n_svd))
    svd_Vtr = svd_Ur
    do i = 1 , n
    do j = 1 , n_svd
        svd_Vtr(i,j) = svd_Vtr(i,j) * svd_S(j)**(-1)
    enddo
    enddo

    svd_Ur = svd_U(1:n,1:n_svd)
    ! Multiply F = B'' * U^H
    call gemm(svd_Vtr, svd_Ur, F,transa="N",transb="C" )
!    call gemm(B(1:n,1:n_svd), svd_U(1:n,1:n_svd), F,transa="N",transb="C" )

    deallocate(svd_U )
    deallocate(svd_S )
    deallocate(svd_Vt)
end subroutine calc_bloch_matrix

! ---------------------------------------------------------- !
! Solves generalized eigenvalue problem of form:
! (     0     1   ) ( c ) =         ( 1    0   ) ( c )
! (   -Tau    H0  ) ( d ) = \lambda ( 0  Tau^+ ) ( d )
! where d = \lambda * c
! by converting it to  standard one using SVD decomposition
! of Tau^+ = U * S * V^+
! ---------------------------------------------------------- !
subroutine try_svd_modes_decomposition(H0,TauDag,ALPHA,BETA,Z)
    complex*16 :: H0(:,:),TauDag(:,:),ALPHA(:),BETA(:),Z(:,:)

    complex*16,allocatable      :: svd_U(:,:),svd_Vt(:,:)
    doubleprecision,allocatable :: svd_S(:)
    complex*16,allocatable      :: svd_hcU(:,:),svd_hcVt(:,:)
    complex*16,allocatable      :: MA(:,:),MB(:,:),Diag(:,:),Tau(:,:),DiagS(:,:),tmpMat(:,:)
    complex*16,allocatable      :: Kaa(:,:),Kac(:,:),Kca(:,:),tmp_Kca(:,:),Kcc(:,:),Kcc_INV(:,:),Kmm(:,:),KMA(:,:)
    complex*16,allocatable      :: geev_a(:),geev_b(:)

    doubleprecision ::   cond_Value
    integer :: n2 , n1 , n_svd , n2_svd
    integer :: i,j,info

    if(QSYS_DEBUG_LEVEL>0) then
        print*,"SYS::LEAD::Finding modes using SVD decomposition of Generalized eigenvalue problem."
    endif
    n1 = size(TauDag,1) ! Size of coupling matrix
    n2 = n1*2           ! Size of Generalized eigenvalue problem
    allocate(MA  (n2,n2))
    allocate(MB  (n2,n2))
    allocate(Diag(n1,n1))
    allocate(DiagS(n1,n1))
    allocate(Tau(n1,n1))
    allocate(tmpMat(n1,n1))
    allocate(svd_hcU(n1,n1))
    allocate(svd_hcVt(n1,n1))
    allocate(svd_U(n1,n1))
    allocate(svd_Vt(n1,n1))
    allocate(svd_S(n1))

    ! Perform SVD of coupling matrix Tau = svd_U * svd_S * svd_Vt
    Tau = TauDag
    call gesvd(Tau,svd_S ,u=svd_U ,vt=svd_Vt ,job="All" ,info=info)
    if(info /= 0 ) then
        print*,"LEAD::SVD decomposition error with info=",info
    endif

    Diag = 0
    DiagS= 0
    Tau  = 0
    MA   = 0
    MB   = 0
    ! Calculate \lambda =\infty subspace size: number of zero rows in Tau matrix
    ! Prepare other matrices:
    n_svd= 0
    do i = 1 , n1
        Diag(i,i)  = 1
        DiagS(i,i) = svd_S(i)

        if( svd_S(i)/svd_S(1) > QSYS_DELTA_SVD ) then
            n_svd = n_svd + 1
        endif
    do j = 1 , n1
        svd_hcU(i,j)  = conjg(svd_U(j,i))
        svd_hcVt(i,j) = conjg(svd_Vt(j,i))
        Tau(i,j)      = conjg(TauDag(j,i))
    enddo
    enddo

    ! Create the generalized eigenvalue problem of form:
    ! (     0         V    ) ( c ) =         ( 1  0 ) ( c )
    ! (-U^+Tsu^+   U^+H0 V ) (~d ) = \lambda ( 0  S ) (~d )
    ! where: ~d = V^+ d, and  d = \lambda c
    ! and S is a diagonal matrix from SVD decomposition.
    call gemm(svd_hcU, Tau, tmpMat )
    MA(n1+1:2*n1,1:n1)      = -tmpMat
    MA(1:n1,n1+1:2*n1)      = svd_hcVt

    call gemm(svd_hcU,        H0, tmpMat )
    call gemm(tmpMat , svd_hcVt , svd_hcU )
    MA(n1+1:2*n1,n1+1:2*n1) = svd_hcU

    MB(1:n1,1:n1)           = Diag
    MB(n1+1:2*n1,n1+1:2*n1) = DiagS

    n_svd = n1 + n_svd
    allocate(Kmm(n_svd,n_svd))
    allocate(geev_a(n_svd))
    allocate(geev_b(n_svd))
    allocate(Kaa(n_svd,n_svd))
    allocate(KMA(n_svd,n_svd))
    allocate(Kac(n_svd,n2-n_svd))
    allocate(Kca(n2-n_svd,n_svd))
    allocate(tmp_Kca(n2-n_svd,n_svd))
    allocate(Kcc(n2-n_svd,n2-n_svd))
    allocate(Kcc_INV(n2-n_svd,n2-n_svd))

    ! Divide matrix into submatrices:
    ! ( Kaa Kac ) ( pa ) =         ( Kmm  0 ) ( pa )
    ! ( Kca Kcc ) ( pb ) = \lambda (  0   0 ) ( pb )
    ! Where pa vector of size (n_svd)
    !   and pb vector of size (n2-n_svd)
    !   matrix Kaa and Kmm has shape (n_svd,n_svd)
    Kmm = MB(1:n_svd,1:n_svd)
    Kaa = MA(1:n_svd,1:n_svd)
    Kac = MA(1:n_svd,1+n_svd:n2)
    Kca = MA(1+n_svd:n2,1:n_svd)
    Kcc = MA(1+n_svd:n2,1+n_svd:n2)

    ! Such system can be transformed to standard eigenvalue
    ! problem if cond(Kcc) < small enough, i.e. is well conditioned
    cond_Value = alg_cond(Kcc)

    B_SINGULAR_MATRIX = .false.
    if( cond_Value*qsys_double_error > 1.0D-6 ) then
        print*,"LEAD::SVD decomposition cannot be use. cond(Kcc)=",cond_Value
        print*,"TRY TO USE FORCE SCHUR DECOMPOSITION MODE"
        deallocate(geev_a  )
        deallocate(geev_b  )
        deallocate(MA  )
        deallocate(MB  )
        deallocate(Diag)
        deallocate(tmpMat)
        deallocate(Tau )
        deallocate(svd_hcU )
        deallocate(svd_hcVt)
        deallocate(svd_U )
        deallocate(svd_S )
        deallocate(svd_Vt)
        deallocate(Kmm)
        deallocate(Kaa)
        deallocate(Kac)
        deallocate(tmp_Kca)
        deallocate(Kca)
        deallocate(Kcc)
        deallocate(Kcc_INV)
        B_SINGULAR_MATRIX = .true. ! This will be used to chose proper method for modes finding i.e:GGEV or Schur Decomp.
        return
    endif
    Kcc     = MA(1+n_svd:n2,1+n_svd:n2)

    ! Invert Kcc matrix
    Kcc_INV = Kcc
    call inverse_matrix(n2-n_svd,Kcc_INV)

    ! Redefine: Kca := Kcc^-1 * Kca
    call gemm(Kcc_INV, Kca, tmp_Kca )
    Kca = tmp_Kca
    ! Calculate: MA := Kaa - Kac * Kcc^-1 * Kca
    call gemm(Kac, Kca, KMA )
    KMA = Kaa - KMA
    ! Multiply by inverse ob B matrix which is diagonal matrix
    ! Generalized eigenvalue problem is transformed to standard one.
    do i = 1 , n_svd
    do j = 1 , n_svd
            KMA(i,j) = KMA(i,j)/KMM(i,i)
    enddo
    enddo

    ! Here we again perform SVD in order to remove singular values
    ! but this time with |lambda| = 0
    deallocate(svd_S )
    deallocate(svd_U )
    deallocate(svd_Vt )
    deallocate(DiagS )
    allocate(svd_S(n_svd) )
    allocate(svd_U(n_svd,n_svd) )
    allocate(svd_Vt(n_svd,n_svd) )
    allocate(DiagS(n_svd,n_svd))
    ! From SVD of matrix A from A x = \lambda x we will get only the
    ! upper part of factorized matrix.
    call gesvd(KMA,svd_S,u=svd_U ,vt=svd_Vt ,job="All" ,info=info)

    ! Finding singular values
    n2_svd = 0
    DiagS  = 0
    do i = 1, n_svd
        if( svd_S(i)/svd_S(1) > QSYS_DELTA_SVD ) then
            n2_svd = n2_svd + 1
            DiagS(i,i) = svd_S(i)
        endif
    enddo
    ! Create new matrix A' = S * V^+ * U
    call gemm(DiagS,svd_Vt,KMA)
    call gemm(KMA,svd_U,svd_Vt)


    if(QSYS_DEBUG_LEVEL > 1) then
        print*,"QSYS::LEAD::SVD cond(Kcc)          =",cond_Value
        print*,"            SVD compression factor =",dble(n2_svd)/n1/2.0
        print*,"            SVD no. singular values=",n2 - n2_svd
        print*,"            SVD cond(Kmm)          =",alg_cond(Kmm)
    endif

    ! And solve eigenproblem for matrix A'

    deallocate(geev_b)
    deallocate(KMA)
    deallocate(Kcc)
    allocate(KMA(1:n2_svd,1:n2_svd))
    allocate(Kcc(1:n2_svd,1:n2_svd))
    allocate(geev_b(1:n2_svd))

    KMA = svd_Vt(1:n2_svd,1:n2_svd)
    call geev(KMA, geev_b , vr=Kcc,info=INFO)


    DiagS                    = 0
    geev_a                   = 0
    geev_a(1:n2_svd)         = geev_b
    ALPHA(1:n_svd)           = geev_a
    DiagS(1:n2_svd,1:n2_svd) = Kcc

    ! Transform back to original space by multiplying result with
    ! X = U_svd * X'
    call gemm(svd_U,DiagS,svd_Vt)
    Kmm(1:n_svd,1:n_svd) = svd_Vt




    Z = 0

    BETA = 1
    BETA(n_svd+1:n2) = 0.0 ! set modes with lambda=infty

    ! Calculate rest of the eigenvectors from
    ! pb(i) = - Kcc^-1 * Kca * pa(i) := - Kca * Kmm
    ! Here we directly locate pb in Kmm matrix Kmm[values,mode]

    call gemm(Kca, Kmm, tmp_Kca )

    Z(1:n_svd,1:n_svd)    =  Kmm
    Z(n_svd+1:n2,1:n_svd) = -tmp_Kca

    ! Back to original basis: d = V * ~d
    ! note that V = svd_hcVt
    deallocate(MA)
    deallocate(MB)
    allocate(MA(n1,n2))
    allocate(MB(n1,n2))

    MA = Z(1+n1:n2,:)
    call gemm(svd_hcVt, MA, MB )
    Z(1+n1:n2,:) = MB

    deallocate(geev_a  )
    deallocate(geev_b  )
    deallocate(MA  )
    deallocate(MB  )
    deallocate(Diag)
    deallocate(Tau )
    deallocate(tmpMat )
    deallocate(svd_hcU )
    deallocate(svd_hcVt)
    deallocate(svd_U )
    deallocate(svd_S )
    deallocate(svd_Vt)
    deallocate(Kmm)
    deallocate(Kaa)
    deallocate(Kac)
    deallocate(tmp_Kca)
    deallocate(Kca)
    deallocate(Kcc)
    deallocate(Kcc_INV)
endsubroutine try_svd_modes_decomposition



endmodule modlead
