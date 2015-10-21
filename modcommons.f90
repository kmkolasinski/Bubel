module modcommons
implicit none

integer,parameter,public  :: QSYS_NO_BONDS_INC_VALUE      = 10     ! For sparse structures
integer,parameter,public  :: QSYS_NO_ATOMS_INC_VALUE      = 10000  !
logical,public            :: QSYS_USE_ZGGEV_TO_FIND_MODES = .false. ! When finding modes use generalied eigenvalue problem, can be more stable
logical,public            :: QSYS_DISABLE_HERMICITY_CHECK = .false.
integer,public            :: QSYS_DEBUG_LEVEL = 0   ! 0 - less messages, 1-more, 2-even more

logical,public            :: B_SINGULAR_MATRIX = .false.
private
! -----------------------------------------------
! Stucture which holds connection between
! Atom A and B.
! -----------------------------------------------
type qbond
    integer     :: toAtomID        ! ID of atom B
    complex*16,allocatable,dimension(:,:)  :: bondMatrix
    complex*16,allocatable,dimension(:,:)  :: overlapMatrix ! overlap matrix in case of LCAO orbitals
    contains
    procedure, public, pass(this) :: destroy_bond
    procedure, public, pass(this) :: copy_bond
endtype qbond

! -----------------------------------------------
! Stucture which holds the information about
! atom
! -----------------------------------------------
type qatom

    doubleprecision :: atom_pos(3)  ! position (x,y,z) in some units
    integer         :: no_in_states ! number of internal states (e.g. spin degree of freedom)
    logical         :: bActive      ! if the site is taken into calculations
    integer         :: flag , flag_aux0         ! arbitrary number, can be used by user e.g. to distinguish two different atoms
    integer,allocatable,dimension(:)     :: globalIDs ! in case of no_bonds > 1 this array contains global ID of atom in spin state
    type(qbond),allocatable,dimension(:) :: bonds     ! contains information about hoping between different atoms
                                                     ! Atom A may have connection with itself
    integer         :: no_bonds                      ! nuber of conetions with different atoms

    contains
    procedure, public, pass(site) :: init
    procedure, public, pass(site) :: destroy
    procedure, public, pass(site) :: add_bond

end type qatom



! ----------------------------------------------------------------
! Structure responsible for  Nearest neigthbour search parameter
! ----------------------------------------------------------------
ENUM , BIND(C)
  ENUMERATOR :: QSYS_NNB_FILTER_BOX       = 1
  ENUMERATOR :: QSYS_NNB_FILTER_CHECK_ALL = 2
  ENUMERATOR :: QSYS_NNB_FILTER_DISTANCE  = 3
END ENUM


type nnb_params
    doubleprecision :: box(3) ! estimated search distance in XYZ directions
    doubleprecision :: distance        ! if NNB_FILTER = QSYS_NNB_FILTER_DISTANCE compare distance
                                       ! between atoms, not coordinates
    integer         :: NNB_FILTER = QSYS_NNB_FILTER_BOX
endtype nnb_params

public :: QSYS_NNB_FILTER_CHECK_ALL
public :: QSYS_NNB_FILTER_DISTANCE
public :: QSYS_NNB_FILTER_BOX
public :: qatom , nnb_params

contains

! ------------------------------------------------------------------------
! Initialize qatom structure. This function does not
! have to be called. Each parameter of atom structure
! can be accessed separetely.
! atom_pos     - position of atom in space. Units are uniportant.
! no_in_states - [optional] number of spin states. Default value is 1.
!                In case of {-1/2,+1/2} electron spin set it to 2.
! bActive      - [optional] each atom can be disactivated before final construction
!                of the lattice. If bActive == false then dis atom will
!                not be taken during the hamiltonian construction.
! flag         - [optional] can be used by user to perform some specific action
! ------------------------------------------------------------------------
subroutine init(site,atom_pos,no_in_states,bActive,flag,flag_aux0)
    class(qatom)     :: site
    doubleprecision :: atom_pos(3)
    integer, optional :: no_in_states , flag ,flag_aux0
    logical, optional :: bActive



    site%no_in_states = 1
    site%bActive      = .true.
    site%atom_pos     = atom_pos
    site%no_bonds     = 0
    site%flag         = 0
    site%flag_aux0    = 0
    if(present(no_in_states)) site%no_in_states  = no_in_states
    if(present(bActive))      site%bActive       = bActive
    if(present(flag))         site%flag          = flag
    if(present(flag_aux0))    site%flag_aux0     = flag_aux0

end subroutine init

! ------------------------------------------------------------------------
! Free allocated memory
! ------------------------------------------------------------------------
subroutine destroy(site)
    class(qatom)  :: site
    integer :: b
    if(allocated(site%globalIDs)) deallocate(site%globalIDs)

    do b = 1,site%no_bonds
        call site%bonds(b)%destroy_bond()
    enddo
    if(allocated(site%bonds))     deallocate(site%bonds)
    site%bActive      = .false.
    site%no_bonds     = 0
    site%flag         = 0
    site%flag_aux0    = 0
end subroutine destroy

! ------------------------------------------------------------------------
! Add new qbonding between two atoms (hoping between A and B).
! fromInnerID - id of spin state of current atom
! toAtomID    - id of atom B
! toInnerID   - id of spin state of atom B
! bondValue   - hoping paremeter
! ------------------------------------------------------------------------
subroutine add_bond(site,toAtomID,bondMatrix,overlapMatrix)
    class(qatom)          :: site
    integer               :: toAtomID
    complex*16,dimension(:,:)              :: bondMatrix
    complex*16,dimension(:,:) , optional  :: overlapMatrix
    integer :: b,nb,ns1,ns2

    ! temporal array
    type(qbond),allocatable,dimension(:) :: tmp_bonds


    ! increase number of bond in atoms
    site%no_bonds = site%no_bonds+1
    nb = size(site%bonds)
    ! adding new bond requires resizing of the bonds array if necessary
    if(.not. allocated(site%bonds)) then
        allocate(site%bonds(QSYS_NO_BONDS_INC_VALUE))
    else if( site%no_bonds > nb ) then


        allocate(tmp_bonds(nb))
        do b = 1 , nb
            call tmp_bonds(b)%copy_bond(site%bonds(b))
            call site%bonds(b)%destroy_bond()
        enddo
        if(allocated(site%bonds)) deallocate(site%bonds)

        allocate(site%bonds(size(tmp_bonds)+QSYS_NO_BONDS_INC_VALUE))

        do b = 1 , nb
            call site%bonds(b)%copy_bond(tmp_bonds(b))
            call tmp_bonds(b)%destroy_bond()
        enddo
        deallocate(tmp_bonds)
    endif


    ! set new bond
    ns1 = size(bondMatrix,1)
    ns2 = size(bondMatrix,2)
    allocate(site%bonds(site%no_bonds)%bondMatrix(ns1,ns2))
    allocate(site%bonds(site%no_bonds)%overlapMatrix(ns1,ns2))

    site%bonds(site%no_bonds)%bondMatrix    = bondMatrix
    site%bonds(site%no_bonds)%toAtomID      = toAtomID

    site%bonds(site%no_bonds)%overlapMatrix = 0.0
    if(present(overlapMatrix)) site%bonds(site%no_bonds)%overlapMatrix = overlapMatrix;

endsubroutine add_bond

subroutine destroy_bond(this)
    class(qbond)  :: this
    this%toAtomID = -1
    if(allocated(this%overlapMatrix)) deallocate(this%overlapMatrix)
    if(allocated(this%bondMatrix))    deallocate(this%bondMatrix)

end subroutine destroy_bond

subroutine copy_bond(this,source)
    class(qbond)  :: this
    type(qbond) :: source
    integer :: ns1,ns2

    ns1 = size(source%bondMatrix,1)
    ns2 = size(source%bondMatrix,2)
    call this%destroy_bond()

    allocate(this%bondMatrix   (ns1,ns2))
    allocate(this%overlapMatrix(ns1,ns2))

    this%toAtomID      = source%toAtomID
    this%bondMatrix    = source%bondMatrix
    this%overlapMatrix = source%overlapMatrix

end subroutine copy_bond

end module modcommons
