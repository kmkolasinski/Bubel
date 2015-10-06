! ---------------------------------------------------------------------------------------
!  Created by Krzysztof Kolasinski - 2015
!
!
! ---------------------------------------------------------------------------------------
module modatom
private
! -----------------------------------------------
! Stucture which holds connection between
! Atom A and B.
! -----------------------------------------------
type qbond
    integer     :: fromInnerID ! hoping from source atom with spin at ID fromInnerID
    integer     :: toAtomID    ! ID of atom B
    integer     :: toInnerID   ! ID of spin of atom B
    complex*16  :: bondValue   ! hoping parameter
endtype qbond

! -----------------------------------------------
! Stucture which holds the information about
! atom
! -----------------------------------------------
type qatom

    doubleprecision :: atom_pos(3)  ! position (x,y,z) in some units
    integer         :: no_in_states ! number of internal states (e.g. spin degree of freedom)
    logical         :: bActive      ! if the site is taken into calculations
    integer         :: flag         ! arbitrary number, can be used by user e.g. to distinguish two different atoms
    integer,allocatable,dimension(:)    :: globalIDs ! in case of no_bonds > 1 this array contains global ID of atom in spin state
    type(qbond),allocatable,dimension(:) :: bonds     ! contains information about hoping between different atoms
                                                     ! Atom A may have connection with itself
    integer         :: no_bonds                      ! nuber of conetions with different atoms

    contains
    procedure, public, pass(site) :: init
    procedure, public, pass(site) :: destroy
    procedure, public, pass(site) :: add_bond!(site,atomID,innerID,bondValue)

end type qatom
public :: qatom

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
subroutine init(site,atom_pos,no_in_states,bActive,flag)
    class(qatom)     :: site
    doubleprecision :: atom_pos(3)
    integer, optional :: no_in_states , flag
    logical, optional :: bActive

    site%no_in_states = 1
    site%bActive      = .true.
    site%atom_pos     = atom_pos
    site%no_bonds     = 0
    site%flag         = 0
    if(present(no_in_states)) site%no_in_states  = no_in_states
    if(present(bActive))      site%bActive       = bActive
    if(present(flag))         site%flag          = flag

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
subroutine add_bond(site,fromInnerID,toAtomID,toInnerID,bondValue)
    class(qatom)     :: site
    integer         :: toAtomID,fromInnerID,toInnerID
    complex*16      :: bondValue

    ! temporal array
    type(qbond),allocatable,dimension(:) :: tmp_bonds
    ! adding new bond requires resizing of the bonds array
    if(site%no_bonds > 0) then
        ! allocated tmp array and copy current array to it
        allocate(tmp_bonds(site%no_bonds))
        tmp_bonds = site%bonds
        ! deallocate bonds and allocate with new size , restore the values
        if(allocated(site%bonds)) deallocate(site%bonds)
        allocate(site%bonds(site%no_bonds+1))
        site%bonds(1:site%no_bonds) = tmp_bonds
        deallocate(tmp_bonds)
    else
        if(allocated(site%bonds)) deallocate(site%bonds)
        allocate(site%bonds(site%no_bonds+1))
    endif
    ! increase number of bond in atoms
    site%no_bonds = site%no_bonds+1
    ! set new bond
    site%bonds(site%no_bonds)%bondValue     = bondValue
    site%bonds(site%no_bonds)%fromInnerID   = fromInnerID
    site%bonds(site%no_bonds)%toAtomID      = toAtomID
    site%bonds(site%no_bonds)%toInnerID     = toInnerID
endsubroutine add_bond


end module modatom


! ---------------------------------------------------------------------------------------
! Main module responsible for system construction and eigenvalue problem for that system.
! ---------------------------------------------------------------------------------------
module modsys
use modatom
use modutils
implicit none
private
integer,parameter :: QSYS_NO_ATOMS_INC_VALUE      = 10000
logical           :: QSYS_DISABLE_HERMICITY_CHECK = .false.

! ----------------------------------------------------------------
! Structure responsible for  Nearest neigthbour search parameter
! ----------------------------------------------------------------
type nnb_params
    doubleprecision :: max_distance(3) ! estimated search distance in XYZ directions
    logical         :: check_all_atoms ! if true then all possible connections are tested
                                       ! during the hamiltonian construction
endtype nnb_params


! ----------------------------------------------------------------
! Structure responsible for most important lattice generation
! and eigenvalue problem solution for that lattice.
! ----------------------------------------------------------------
type qsys
    type(qatom) , allocatable, dimension(:) :: atoms ! arrays of atoms

    integer     :: no_atoms         ! total number of atoms (active + notActive)
    integer     :: system_size      ! size of the hamiltonian matrix. Number of unknown in linear system equation

    complex*16     ,allocatable,dimension(:,:) :: eigenvecs ! matrix containing eigenvectors. eigenvecs(global_id,state)
    doubleprecision,allocatable,dimension(:  ) :: eigenvals ! array containg eivenvalues. eigenvals(state)
    integer :: no_eigenvalues                               ! number of stated which has been found by FEAST library
                                                            ! if 0 - no state was found

    ! other variables
    type(qatom)        :: qatom      ! auxiliary variable, can be used by user
    type(nnb_params)  :: qnnbparam  ! auxiliary variable, can be used by user

    contains
    procedure, public, pass(sys) :: init
    procedure, public, pass(sys) :: destroy
    procedure, public, pass(sys) :: add_atom!(sys,site)
    procedure, public, pass(sys) :: make_lattice!(sys,connect_procedure)
    procedure, public, pass(sys) :: save_lattice!(filename,innerA,innerB)
    procedure, public, pass(sys) :: calc_eigenproblem!(sys,pEmin,pEmax,NoStates,no_feast_contours,print_info,pmaks_iter)

endtype qsys
public :: qsys,nnb_params
public :: QSYS_DISABLE_HERMICITY_CHECK
contains

! ------------------------------------------------------------------------
! Initialize system. No arguments are required.
! ------------------------------------------------------------------------
subroutine init(sys)
    class(qsys)     :: sys
    print*,"SYS::Initialization of the system."
    sys%no_atoms         = 0
    sys%qnnbparam%check_all_atoms = .false.
    sys%qnnbparam%max_distance    = 0.0D0

    allocate(sys%atoms(QSYS_NO_ATOMS_INC_VALUE))
end subroutine init

! ------------------------------------------------------------------------
! Free memory allocated by qsystem.
! ------------------------------------------------------------------------
subroutine destroy(sys)
    class(qsys)     :: sys
    integer :: i
    print*,"SYS::Freeing memory..."

    ! free memory allocated by atoms first.
    do i = 1 , sys%no_atoms
        call sys%atoms(i)%destroy()
    enddo
    sys%no_atoms      = 0
    sys%no_eigenvalues= 0
    if(allocated(sys%atoms))     deallocate(sys%atoms)
    if(allocated(sys%eigenvecs)) deallocate(sys%eigenvecs)
    if(allocated(sys%eigenvals)) deallocate(sys%eigenvals)

end subroutine destroy

! ------------------------------------------------------------------------
! Add atom to the system. This function add the atom to the atoms array.
! site - atom structer object. All the properties of the atom can be accessed
!        from the atoms array.
! ------------------------------------------------------------------------
subroutine add_atom(sys,site)
    class(qsys)           :: sys
    type(qatom),intent(in) :: site
    type(qatom),allocatable,dimension(:) :: tmp_atoms


    ! increase the atom index
    sys%no_atoms = sys%no_atoms + 1
    ! when the number of atoms is greated than atoms array size
    ! one needs to resize the array
    if( sys%no_atoms > size(sys%atoms) ) then
        print*,"SYS::INFO::Resizing the system from:",size(sys%atoms)," to ",size(sys%atoms)+QSYS_NO_ATOMS_INC_VALUE
        allocate(tmp_atoms(size(sys%atoms)))
        tmp_atoms = sys%atoms
        deallocate(sys%atoms)
        allocate(sys%atoms(size(tmp_atoms)+QSYS_NO_ATOMS_INC_VALUE))
        sys%atoms(1:size(tmp_atoms)) = tmp_atoms
        deallocate(tmp_atoms)
    endif
    ! add new atom to the system
    sys%atoms(sys%no_atoms) = site

end subroutine add_atom


! ------------------------------------------------------------------------
! When all the atoms are created use this function to created connection
! between the atoms. After execution of this function no new atoms should
! added.
! Additionally this function check if created connection satisfy the
! hermicity of the Hamiltonian matrix. IF no program is stopped.
! If there is no need for Hamiltonina to be hermitian one may disable this
! check with changing the variable QSYS_DISABLE_HERMICITY_CHECK.
! Parameters:
! connect_procedure - user provided procedure which check if atomA in spin
!                     state s1 has coupling with atomB in spin state s2.
!                     If not function returns false. If yes function returns
!                     true and override the coupling_val argument with
!                     propper value of hopin energy.
! nnbparams         - nnbparams structure which contains some information
!                     how to search for the nearest neightbours atoms.
! ------------------------------------------------------------------------
subroutine make_lattice(sys,connect_procedure,nnbparams)
    class(qsys) :: sys
    interface
        logical function connect_procedure(atomA,atomB,s1,s2,coupling_val)
            use modatom
            implicit none
            type(qatom) :: atomA,atomB
            integer    :: s1,s2
            complex*16 :: coupling_val
        end function connect_procedure
    endinterface
    type(nnb_params) :: nnbparams

    ! local variables
    integer         :: i,j,k,l,p,no_active_atoms
    doubleprecision :: time_start
    complex*16      :: cpl_value,cpl_value_inverse
    ! bounding box parameters
    integer,parameter :: ix = 1 , iy = 2 , iz = 3 , cmin = 1 , cmax =  2
    double precision  :: bbox(cmin:cmax,ix:iz)
    ! verlet box searching method
    integer              :: verlet_dims(3) , verlet_max_atoms_in_cell , ivpos(3) , vc
    integer              :: vx,vy,vz,vp,vnnba
    integer, allocatable :: verlet_box(:,:,:,:) , verlet_counter(:,:,:)



    time_start = get_clock()

    ! Calulation of the real number of active atoms.
    no_active_atoms = 0
    do i = 1 , sys%no_atoms
        if(sys%atoms(i)%bActive) no_active_atoms = no_active_atoms + 1
    enddo
    print*,"SYS::Calculating connections between:",no_active_atoms," atoms."

    ! setting up global IDs for atoms and inner states (spins)
    k = 0
    do i = 1 , sys%no_atoms
        if(sys%atoms(i)%bActive) then
            if(allocated(sys%atoms(i)%globalIDs)) deallocate(sys%atoms(i)%globalIDs)
            allocate(sys%atoms(i)%globalIDs(sys%atoms(i)%no_in_states))
            do j = 1 , sys%atoms(i)%no_in_states
                k = k + 1
                sys%atoms(i)%globalIDs(j) = k
            enddo
        endif ! end of active atom
    enddo
    sys%system_size = k ! the number of unknowns
    print*,"SYS::Number of uknown variables     :",k

    ! Making the connection between the atoms
    if(nnbparams%check_all_atoms) then
    ! if all posible connections are tested enter this node
    do i = 1 , sys%no_atoms
    if(sys%atoms(i)%bActive == .true.) then
        do k = 1 , sys%atoms(i)%no_in_states
        do j = 1 , sys%no_atoms
            ! both sites have to be active
            if(sys%atoms(j)%bActive == .true.) then
                do l = 1 , sys%atoms(j)%no_in_states
                    ! if they are nnb one may create a bond between them
                    if(connect_procedure(sys%atoms(i),sys%atoms(j),k,l,cpl_value)) then
                         call sys%atoms(i)%add_bond(k,j,l,cpl_value)
                         print"(A,2i6,A,2i6,A,1f6.2)","creating bond: A(",i,k,") -> B(",j,l,") = ",dble(cpl_value)
!                         print*,"   pos=",sys%atoms(j)%atom_pos(1:2)
                    endif
                enddo
            endif
        enddo
        enddo ! end of k
    endif ! end of i bactive i
    enddo ! end of i
    ! ---------------------------------------------------------------------
    ! Else use verlet method to find nnb's
    ! ---------------------------------------------------------------------
    else
    ! calculate the bounding box of the system
    bbox(cmin,:) = sys%atoms(1)%atom_pos ! minimum position
    bbox(cmax,:) = sys%atoms(1)%atom_pos ! maximum position
    do i = 1 , sys%no_atoms
        if(sys%atoms(i)%bActive) then
        do k = ix , iz
        if( sys%atoms(i)%atom_pos(k) < bbox(cmin,k) ) bbox(cmin,k) = sys%atoms(i)%atom_pos(k)
        if( sys%atoms(i)%atom_pos(k) > bbox(cmax,k) ) bbox(cmax,k) = sys%atoms(i)%atom_pos(k)
        enddo
        endif ! end of if active atom
    enddo

    ! Calculate dimensions of the verlet box
    do i = ix , iz
        verlet_dims(i) = NINT(abs(bbox(cmax,i) - bbox(cmin,i))/abs(nnbparams%max_distance(i)))+1
        if( abs(nnbparams%max_distance(i)) < 1.0D-20 ) verlet_dims(i) = 1
    enddo

    !!Debug:
    !print*,"Verlet:",verlet_dims

    ! allocate array which holds the number of atoms in each verlet cell and clear values
    allocate(verlet_counter(verlet_dims(1),verlet_dims(2),verlet_dims(3)))
    verlet_counter           = 0
    verlet_max_atoms_in_cell = 0 ! calulate maximum number of atoms in one verlet cell

    ! count the number of atoms in each verlet cell
    do i = 1 , sys%no_atoms
    if(sys%atoms(i)%bActive) then ! only active atom is taken into acount
        do k = 1 , 3
        ivpos(k) = 1+FLOOR(verlet_dims(k)*(sys%atoms(i)%atom_pos(k) - bbox(cmin,k))/(bbox(cmax,k)-bbox(cmin,k)+1.0D-6))
        enddo
        verlet_counter(ivpos(1),ivpos(2),ivpos(3)) = verlet_counter(ivpos(1),ivpos(2),ivpos(3)) + 1
        ! search for the maximum number of atoms in cell
        if(verlet_counter(ivpos(1),ivpos(2),ivpos(3)) > verlet_max_atoms_in_cell) verlet_max_atoms_in_cell = verlet_counter(ivpos(1),ivpos(2),ivpos(3))
    endif
    enddo

    !!Debug:
    !print*,"Sum atoms in verlet:",sum(verlet_counter)
    !print*,"Max atoms in cell:",verlet_max_atoms_in_cell

    ! Allocate verlet cells. The place we will put atoms
    allocate(verlet_box(verlet_dims(1),verlet_dims(2),verlet_dims(3),verlet_max_atoms_in_cell))

    ! Perform a test which checks if number of active atoms if same as in the verlet method
    ! If not program has to be stopped.
    if(sum(verlet_counter) /= no_active_atoms) then
        print*,"SYS::ERROR::Verlet system paritioning failed. The checksum does not equal total number of active atom in system."
        print*,"            Total number atoms  :",no_active_atoms
        print*,"            Number from checksum:",sum(verlet_counter)
        print*,"SYS::ERROR::The program has stopped. Use can try to set nnbparam ", &
               "            to check_all_atoms =  true if the considered system is small ",&
               "            or change max_distance parameters."
        stop -1
    endif

    ! putting atoms to verlet boxes
    verlet_counter = 0
    do i = 1 , sys%no_atoms
    if(sys%atoms(i)%bActive) then
        do k = 1 , 3
        ivpos(k) = 1+FLOOR(verlet_dims(k)*(sys%atoms(i)%atom_pos(k) - bbox(cmin,k))/(bbox(cmax,k)-bbox(cmin,k)+1.0D-6))
        enddo
        verlet_counter(ivpos(1),ivpos(2),ivpos(3)) = verlet_counter(ivpos(1),ivpos(2),ivpos(3)) + 1
        vc = verlet_counter(ivpos(1),ivpos(2),ivpos(3))
        ! save the atom ID in the verlet box
        verlet_box(ivpos(1),ivpos(2),ivpos(3),vc) = i
    endif
    enddo

    ! find nnb using verlet arrays
    do i = 1 , sys%no_atoms
        if(sys%atoms(i)%bActive) then
        do p = 1 , 3
        ivpos(p) = 1+FLOOR(verlet_dims(p)*(sys%atoms(i)%atom_pos(p) - bbox(cmin,p))/(bbox(cmax,p)-bbox(cmin,p)+1.0D-6))
        enddo

        do k = 1 , sys%atoms(i)%no_in_states

        ! search in nearest verlet cells
        do vx = max(ivpos(1)-1,1),min(ivpos(1)+1,verlet_dims(1))
        do vy = max(ivpos(2)-1,1),min(ivpos(2)+1,verlet_dims(2))
        do vz = max(ivpos(3)-1,1),min(ivpos(3)+1,verlet_dims(3))
            ! loop over nnb verlet boxes
            do vp = 1 , verlet_counter(vx,vy,vz)
                j        = verlet_box(vx,vy,vz,vp)
                if(sys%atoms(j)%bActive == .true.) then
                ! loop around spin states

                do l = 1 , sys%atoms(j)%no_in_states
                    ! if they are nnb one may create a qbond between them
                    if(connect_procedure(sys%atoms(i),sys%atoms(j),k,l,cpl_value)) then
                         call sys%atoms(i)%add_bond(k,j,l,cpl_value)
                    endif
                enddo


                endif
            enddo
        enddo ! }
        enddo ! } nearest verlet cells
        enddo ! }

        enddo ! end of do i atom

        endif ! end of if active atom


    enddo
        ! free verlet arrays
        deallocate(verlet_box)
        deallocate(verlet_counter)
    endif ! end of if else check all atoms - use verlet cells

    ! Perform a check to test the hermiticity of the hamiltonian matrix.
    ! For a given atom A
    ! Take one bond to atom B and check if the value of the hoping parameter
    ! is the same as for the atom B connected with atom A.
    ! If both atoms are in the same spin state.
    if( .not. QSYS_DISABLE_HERMICITY_CHECK) then
    print*,"SYS::INFO::Checking hermiticity of the matrix."
    do i = 1 , sys%no_atoms ! take the atom A
        if(sys%atoms(i)%bActive) then
            do j = 1 , sys%atoms(i)%no_bonds ! take one bond to atom B
                vp = sys%atoms(i)%bonds(j)%toAtomID
                cpl_value = sys%atoms(i)%bonds(j)%bondValue ! take hoping value from A to B
                cpl_value_inverse = 0
                do k = 1 , sys%atoms(vp)%no_bonds ! search for the same
                    if( sys%atoms(vp)%bonds(k)%toAtomID == i .and. &
                     sys%atoms(i)%bonds(j)%toInnerID   == sys%atoms(vp)%bonds(k)%fromInnerID .and. &
                     sys%atoms(i)%bonds(j)%fromInnerID == sys%atoms(vp)%bonds(k)%toInnerID )  then

                     cpl_value_inverse = sys%atoms(vp)%bonds(k)%bondValue
                     exit
                    endif
                enddo
                if( abs(cpl_value - conjg(cpl_value_inverse)) > 10.0D-10 ) then
                    print*,"==============================================================================="
                    print*,"SYS::ERROR::Created matrix will be not Hermitian. Check the provided ", &
                           "            connection function. The qbonding param between atom A->B", &
                           "            has to complex conjugate for connection between B->A."
                    print"(A,3e16.4)","             An error occured for atom A at position :",sys%atoms(i)%atom_pos
                    print"(A,3e16.4)","             and atom B at position                  :",sys%atoms(vp)%atom_pos
                    print*,"            Atom A id=",i ," spin=",sys%atoms(i)%bonds(j)%fromInnerID
                    print*,"            Atom B id=",vp," spin=",sys%atoms(vp)%bonds(k)%fromInnerID
                    print"(A,3e16.4)","             Hoping parameter for A->B:",cpl_value
                    print"(A,3e16.4)","             Hoping parameter for B->A:",cpl_value_inverse
                    print*,"==============================================================================="
                    stop -1
                endif

            enddo
        endif ! end of active atom
    enddo
    endif ! end of if disable hermiticity check


    print*,"SYS::Connections has been found in ", get_clock()-time_start , " sec."

end subroutine make_lattice



! ------------------------------------------------------------------------
! Save current lattice to file. Use this function after execution of make_lattice
! procedure.
! filename - the name of the output file
! innerA   - [optional] plot connection for atoms in state A. Default is 1
! innerB   - [optional] connected with atoms in state B. This can be used to see the
!            connections between different spin stated of atoms. Default is 1.
! ------------------------------------------------------------------------
subroutine save_lattice(sys,filename,innerA,innerB)
    class(qsys)     :: sys
    character(*)    :: filename
    integer,optional:: innerA , innerB

    integer :: s1,s2
    integer :: i,b,ida,ids1,ids2
    integer,parameter :: ix = 1 , iy = 2 , iz = 3 , cmin = 1 , cmax =  2
    double precision :: lWidth , bbox(2,3)
    print*,"SYS::Saving lattice to file:",filename
    s1 = 1
    s2 = 1
    if(present(innerA)) s1 = innerA
    if(present(innerB)) s2 = innerB

    ! calculate the bounding box of the system
    bbox(cmin,:) = sys%atoms(1)%atom_pos ! minimum position
    bbox(cmax,:) = sys%atoms(1)%atom_pos ! maximum position
    do i = 1 , sys%no_atoms
        if(sys%atoms(i)%bActive) then
        do b = ix , iz
        if( sys%atoms(i)%atom_pos(b) < bbox(cmin,b) ) bbox(cmin,b) = sys%atoms(i)%atom_pos(b)
        if( sys%atoms(i)%atom_pos(b) > bbox(cmax,b) ) bbox(cmax,b) = sys%atoms(i)%atom_pos(b)
        enddo
        endif ! end of if active atom
    enddo

    open(unit=765819,file=filename)
    write(765819,*),"# Bounding box to estimate the scale of the system"
    write(765819,*),bbox(cmin,:)
    write(765819,*),bbox(cmax,:)
    write(765819,*),"# Lines which connect the atoms by qbonding."
    write(765819,*),"# Connections between spin ",s1," and ",s2
    do i = 1 , sys%no_atoms
        if(sys%atoms(i)%bActive) then
        do b = 1 , sys%atoms(i)%no_bonds
        ! write only unique connections
        ida  = sys%atoms(i)%bonds(b)%toAtomID
        ids1 = sys%atoms(i)%bonds(b)%toInnerID
        ids2 = sys%atoms(i)%bonds(b)%fromInnerID
        if( ida >= i ) then
            lWidth = 0.1
            if(ida == i) lWidth = 10.0

            if(s1 == ids1 .and. s2 == ids2) then
                write(765819,"(7e20.6)"),sys%atoms(i)%atom_pos,sys%atoms(ida)%atom_pos,lWidth
            endif
        endif
        enddo
        endif ! end of if active atom
    enddo
    close(765819)


end subroutine save_lattice

! ------------------------------------------------------------------------
! Solve eigen problem for generated lattice. Eigen values and eigen vectors
! are saved to sys%eigenvals and sys%eigenvecs arrays.
! Params:
! (pEmin,pEmax)      - energy range search
! NoStates           - expected number of stated in that range
! no_feast_contours  - [optional] Number of countour integrals perfomed by FEAST. Default is 8.
!                      Allowed values are:  {3,4,5,6,8,10,12,16,20,24,32,40,48}. See
!                      FEAST documentation for more info.
! print_info         - [optional] print DEBUG info. Deafault is 0 (i.e. false)
! pmaks_iter         - [optional] maximum number of FEAST iterations. Default is 10.
! ------------------------------------------------------------------------
subroutine calc_eigenproblem(sys,pEmin,pEmax,NoStates,no_feast_contours,print_info,pmaks_iter)
    class(qsys) :: sys
    doubleprecision   :: pEmin, pEmax
    integer           :: NoStates
    integer,optional  :: no_feast_contours,print_info,pmaks_iter


    integer :: fpm(128)
    integer :: i,j,info,itmp,nw,M0,loop,no_evals,ta,ts
    integer :: no_contours,display_info,maks_iter
    doubleprecision :: epsout
    doubleprecision :: Emin, Emax
    doubleprecision :: time_start

    integer,allocatable                          :: HBROWS(:)
    complex*16,dimension(:,:), allocatable       :: EVectors
    complex*16,dimension(:)  , allocatable       :: MATHVALS
    integer   ,dimension(:,:), allocatable       :: ROWCOLID
    double precision,dimension(:), allocatable   :: Evalues,Rerrors

    integer :: NO_NON_ZERO_VALUES , NO_VARIABLES

    time_start = get_clock()
    ! Przejscie do jednostek donorowych
    Emin = pEmin
    Emax = pEmax



    ! setting the default parameters
    if(.not. present(no_feast_contours)) then
        no_contours = 8
    else
        no_contours = no_feast_contours
    endif
    if(.not. present(print_info)) then
        display_info = 0
    else
        display_info = print_info
    endif
    if(.not. present(pmaks_iter)) then
        maks_iter = 20
    else
        maks_iter = pmaks_iter
    endif



    call feastinit(fpm)
    fpm(1)=display_info      ! do not show any information
    fpm(2)=no_contours       ! number of contours
    fpm(3)=12                ! exponent of the error
    fpm(4)=maks_iter         ! maximum number of iteration
    fpm(5)=0                 ! we start with default vector (if 1 then with provided)
    fpm(6)=0                 ! convergece criterium with value of residuum (0 albo 1)

    ! Calculate the number of non-zero elements in Hamiltonian matrix
    itmp = 0
    do i = 1, sys%no_atoms
        if(sys%atoms(i)%bActive) then
            itmp = itmp + sys%atoms(i)%no_bonds
        endif
    enddo
    NO_NON_ZERO_VALUES = itmp
    ! The number of unknows is taken from the last global index
    NO_VARIABLES       = sys%system_size

    print*,"SYS::Calulating eigenvalue problem using FEAST solver"

    allocate(MATHVALS(NO_NON_ZERO_VALUES))
    allocate(ROWCOLID(NO_NON_ZERO_VALUES,2))

    ! Filling matrix and row-col array
    itmp = 0
    do i = 1 ,  sys%no_atoms
        if(sys%atoms(i)%bActive) then
        do j = 1, sys%atoms(i)%no_bonds
        itmp = itmp + 1
        MATHVALS(itmp)   = sys%atoms(i)%bonds(j)%bondValue
        ROWCOLID(itmp,1) = sys%atoms(i)%globalIDs(sys%atoms(i)%bonds(j)%fromInnerID)

        ta = sys%atoms(i)%bonds(j)%toAtomID
        ts = sys%atoms(i)%bonds(j)%toInnerID
        ROWCOLID(itmp,2) = sys%atoms(ta)%globalIDs(ts)
!        print*,ROWCOLID(itmp,:),MATHVALS(itmp)
        enddo
        endif ! end of if active atom
    enddo

    ! auxiliary variable
    nw   = NO_NON_ZERO_VALUES
    if(display_info==1) then
        print*,"--------------------------------------------------"
        print*,"SYS::EIGENVALUE PROBLEM INPUTS"
        print*,"--------------------------------------------------"
        print*,"SYS::FEAST::Number of variables   ",NO_VARIABLES
        print*,"SYS::FEAST::Energy range from     ",Emin," to ",Emax
        print*,"SYS::FEAST::Number of contours    ",no_contours
        print*,"SYS::FEAST::Max. no. of iters     ",maks_iter
        print*,"SYS::FEAST::Expected no. of states",NoStates
    endif


    ! -----------------------------------------------------------
    allocate(HBROWS(NO_VARIABLES+1))
    call convert_to_HB(NO_NON_ZERO_VALUES,ROWCOLID,MATHVALS,HBROWS)

    ! zgadujemy liczbe stanow
    M0  = NoStates
    allocate(EVectors(NO_VARIABLES,M0))
    allocate(Evalues(M0))
    allocate(Rerrors(M0))



    call zfeast_hcsrev('F',&               ! - 'F' oznacza ze podawana jest pelna macierz
                          NO_VARIABLES,&   ! - rozmiar problemu (ile wezlow z flaga B_NORMAL)
                          MATHVALS(1:nw),& ! - kolejne nie zerowe wartosci w macierzy H
                          HBROWS,&         ! - numeracja wierszy (rodzaj zapisu macierzy rzakidch)
                          ROWCOLID(1:nw,2),& ! - indeksy kolumn odpowiadaja tablicy wartosci CMATA
                          fpm,&            ! - wektor z konfiguracja procedury
                          epsout,&         ! - Residuum wyjsciowe
                          loop, &          ! - Koncowa liczba iteracji
                          Emin,&           ! - Minimalna energia przeszukiwania
                          Emax,&           ! - Maksymalna energia
                          M0,&             ! - Spodziewana liczba modow w zakresie (Emin,Emax)
                          Evalues,&        ! - Wektor z otrzymanymi wartosciami wlasnymi
                          EVectors,&       ! - Macierz z wektorami (kolejne kolumny odpowiadaja kolejnym wartoscia z tablicy Evalues)
                          no_evals,&       ! - Liczba otrzymanych wartosci z przedziale (Emin,Emax)
                          Rerrors,&        ! - Wektor z bledami dla kolejnych wartosci wlasnych
                          info)            ! - Ewentualne informacje o bledach



        if(display_info==1) then
            print*,"SYS::FEAST::Output error       ",  epsout
            print*,"SYS::FEAST::No. interations    ",  loop
            print*,"SYS::FEAST::No. states         ",  no_evals
            print*,"SYS::FEAST::Ouput info value   ",  info
            print*,"SYS::FEAST::Calulation time [s]",  get_clock() - time_start
            print*,"--------------------------------------------------"
        endif



        sys%no_eigenvalues = no_evals

        ! ----------------------------------------------------------------------------------
        ! Obsluga bledow:
        ! ----------------------------------------------------------------------------------
        selectcase(info)
        case( 202 )
            print*," Error : Problem with size of the system n (n≤0) "
            stop
        case( 201 )
            print*," Error : Problem with size of initial subspace m0 (m0≤0 or m0>n) "
            stop
        case( 200 )
            print*," Error : Problem with emin,emax (emin≥emax) "
            stop
        case(100:199)
            print"(A,I4,A)"," Error : Problem with ",info-100,"-th value of the input Extended Eigensolver parameter (fpm(i)). Only the parameters in use are checked. "
            sys%no_eigenvalues = 0
            stop
        case( 4 )
            print*," Warning : Successful return of only the computed subspace after call withfpm(14) = 1 "
            sys%no_eigenvalues = 0

        case( 3 )
            print*," Warning : Size of the subspace m0 is too small (m0<m) "
            sys%no_eigenvalues = 0

        case( 2 )
            print*," Warning : No Convergence (number of iteration loops >fpm(4))"
            sys%no_eigenvalues = 0
        case( 1 )
            print*," Warning : No eigenvalue found in the search interval. See remark below for further details. "
            sys%no_eigenvalues = 0
        case( 0 )
            print*,               "---------------------------------------------"
            print"(A,i12)",       " SYS::FEAST:: No. states  :",sys%no_eigenvalues
            print"(A,f12.3)",     "              Time [s]    :",get_clock() - time_start
            print"(A,e12.4)",     "              Error       :",epsout
            print"(A,i12)",       "              No. iters   :",loop
            print*,               "---------------------------------------------"
        case( -1 )
            print*," Error : Internal error for allocation memory. "
            stop
        case( -2 )
            print*," Error : Internal error of the inner system solver. Possible reasons: not enough memory for inner linear system solver or inconsistent input. "
            stop
        case( -3 )
            print*," Error : Internal error of the reduced eigenvalue solver Possible cause: matrix B may not be positive definite. It can be checked with LAPACK routines, if necessary."
            stop
        case(-199:-100)
            print"(A,I4,A)"," Error : Problem with the ",-info-100,"-th argument of the Extended Eigensolver interface. "
            stop
        endselect

        ! -----------------------------------------------------------------
        ! Kopiowanie wynikow do odpowiednich tablic
        ! -----------------------------------------------------------------


        if(allocated(sys%eigenvals)) deallocate(sys%eigenvals)
        if(allocated(sys%eigenvecs)) deallocate(sys%eigenvecs)


        if(sys%no_eigenvalues > 0 ) then
            allocate(sys%eigenvecs(NO_VARIABLES,sys%no_eigenvalues))

            sys%eigenvecs(:,1:sys%no_eigenvalues) = EVectors(:,1:sys%no_eigenvalues)
            allocate(sys%eigenvals(sys%no_eigenvalues))
            sys%eigenvals(1:sys%no_eigenvalues) = Evalues(1:sys%no_eigenvalues)

        endif


        deallocate(MATHVALS)
        deallocate(ROWCOLID)
        deallocate(HBROWS)
        deallocate(EVectors)
        deallocate(Evalues)
        deallocate(Rerrors)

        print*,"SYS::Eigenvalues calculated. Found:",sys%no_eigenvalues
end subroutine calc_eigenproblem



! ------------------------------------------------------------------------
!                             Helper functions
! ------------------------------------------------------------------------


! Conversion from row-col-value sparse matrix storage to
! HB compressed format.
subroutine convert_to_HB(no_vals,rows_cols,matA,out_rows)
      integer,intent(in)                  :: no_vals
      integer   ,intent(inout),dimension(:,:) :: rows_cols
      complex*16,intent(inout),dimension(:) :: matA
      integer,intent(inout),dimension(:)   :: out_rows
      integer :: iterator, irow ,  from , to
      integer :: i, n

      n        = no_vals
      iterator = 0
      irow     = 0
      do i = 1 , n
          if( rows_cols(i,1) /= irow ) then
            iterator = iterator + 1
            out_rows(iterator) = i
            irow = rows_cols(i,1)
          endif
      enddo
      out_rows(iterator+1) = n + 1


!DEC$ IF DEFINED  (USE_UMF_PACK)
  irow = size(out_rows)-1
    ! sortowanie  kolumn
  do i = 1 , irow-1
  from = out_rows(i)
  to   = out_rows(i+1)-1
      call sort_col_vals(rows_cols(from:to,2),matA(from:to))
  enddo

  ! przesuwanie indeksow do zera
  out_rows       = out_rows -1
  rows_cols(:,2) = rows_cols(:,2) -1
!DEC$ ENDIF

!DEC$ IF DEFINED  (USE_PARDISO)

  irow = size(out_rows)-1
    ! sortowanie  kolumn
  do i = 1 , irow-1
  from = out_rows(i)
  to   = out_rows(i+1)-1
      call sort_col_vals(rows_cols(from:to,2),matA(from:to))
  enddo

!DEC$ ENDIF
end subroutine convert_to_HB

subroutine sort_col_vals(cols,vals)
        integer,intent(inout),dimension(:)    :: cols
        complex*16,intent(inout),dimension(:) :: vals
        integer :: tmp_col
        complex*16 :: tmp_val
        integer :: i  , n
        logical :: test
        n = size(cols)

        test = .true.

        ! sortowanie bombelkowe
        do while(test)
          test = .false.
          do i = 1 , n-1
            if( cols(i) > cols(i+1)  ) then
            tmp_col   = cols(i)
            cols(i)   = cols(i+1)
            cols(i+1) = tmp_col

            tmp_val   = vals(i)
            vals(i)   = vals(i+1)
            vals(i+1) = tmp_val

            test = .true.
            exit
            endif
          enddo
        enddo
end subroutine sort_col_vals
end module modsys
