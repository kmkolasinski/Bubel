module modcommons
implicit none



private

integer ,parameter :: QSYS_NO_BONDS_INC_VALUE = 10 ! For sparse structures

! -----------------------------------------------
! Stucture which holds connection between
! Atom A and B.
! -----------------------------------------------
type qbond
    integer     :: fromInnerID    ! hoping from source atom with spin at ID fromInnerID
    integer     :: toAtomID       ! ID of atom B
    integer     :: toInnerID      ! ID of spin of atom B
    complex*16  :: bondValue      ! hoping parameter
    complex*16  :: overlapValue   ! overlap value in case of LCAO orbitals
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
    integer,allocatable,dimension(:)    :: globalIDs ! in case of no_bonds > 1 this array contains global ID of atom in spin state
    type(qbond),allocatable,dimension(:) :: bonds     ! contains information about hoping between different atoms
                                                     ! Atom A may have connection with itself
    integer         :: no_bonds                      ! nuber of conetions with different atoms

    contains
    procedure, public, pass(site) :: init
    procedure, public, pass(site) :: destroy
    procedure, public, pass(site) :: add_bond!(site,atomID,innerID,bondValue)

end type qatom

ENUM , BIND(C)
  ENUMERATOR :: QSYS_NNB_FILTER_BOX       = 1
  ENUMERATOR :: QSYS_NNB_FILTER_CHECK_ALL = 2
  ENUMERATOR :: QSYS_NNB_FILTER_DISTANCE  = 3
END ENUM

! ----------------------------------------------------------------
! Structure responsible for  Nearest neigthbour search parameter
! ----------------------------------------------------------------
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
    class(qatom)     :: site
    if(allocated(site%globalIDs)) deallocate(site%globalIDs)
    if(allocated(site%bonds))     deallocate(site%bonds)
    site%bActive      = .false.
    site%no_bonds     = 0
    site%flag         = 0
end subroutine destroy

! ------------------------------------------------------------------------
! Add new qbonding between two atoms (hoping between A and B).
! fromInnerID - id of spin state of current atom
! toAtomID    - id of atom B
! toInnerID   - id of spin state of atom B
! bondValue   - hoping paremeter
! ------------------------------------------------------------------------
subroutine add_bond(site,fromInnerID,toAtomID,toInnerID,bondValue,overlapValue)
    class(qatom)     :: site
    integer         :: toAtomID,fromInnerID,toInnerID
    complex*16      :: bondValue
    complex*16, optional  :: overlapValue


    ! temporal array
    type(qbond),allocatable,dimension(:) :: tmp_bonds
    ! adding new bond requires resizing of the bonds array
!    if(site%no_bonds > 0) then
!        ! allocated tmp array and copy current array to it
!        allocate(tmp_bonds(site%no_bonds))
!        tmp_bonds = site%bonds
!        ! deallocate bonds and allocate with new size , restore the values
!        if(allocated(site%bonds)) deallocate(site%bonds)
!        allocate(site%bonds(site%no_bonds+1))
!        site%bonds(1:site%no_bonds) = tmp_bonds
!        deallocate(tmp_bonds)
!    else
!        if(allocated(site%bonds)) deallocate(site%bonds)
!        allocate(site%bonds(site%no_bonds+1))
!    endif


    ! increase number of bond in atoms
    site%no_bonds = site%no_bonds+1

    if(.not. allocated(site%bonds)) then
        allocate(site%bonds(QSYS_NO_BONDS_INC_VALUE))
    else if( site%no_bonds > size(site%bonds) ) then
!        print*,"SYS::INFO::Resizing the system from:",size(site%bonds)," to ",size(site%bonds)+QSYS_NO_BONDS_INC_VALUE
        allocate(tmp_bonds(size(site%bonds)))
        tmp_bonds = site%bonds
        if(allocated(site%bonds)) deallocate(site%bonds)

        allocate(site%bonds(size(tmp_bonds)+QSYS_NO_BONDS_INC_VALUE))
        site%bonds(1:site%no_bonds) = tmp_bonds
        deallocate(tmp_bonds)
!    else if(site%no_bonds == 0) then

    endif


    ! set new bond
    site%bonds(site%no_bonds)%bondValue     = bondValue
    site%bonds(site%no_bonds)%fromInnerID   = fromInnerID
    site%bonds(site%no_bonds)%toAtomID      = toAtomID
    site%bonds(site%no_bonds)%toInnerID     = toInnerID
    site%bonds(site%no_bonds)%overlapValue  = 0.0;
    if(present(overlapValue)) site%bonds(site%no_bonds)%overlapValue = overlapValue;

endsubroutine add_bond


end module modcommons
