! ---------------------------------------------------------------------------------------
!  Created by Krzysztof Kolasinski - 2015
!
!
! ---------------------------------------------------------------------------------------


! ---------------------------------------------------------------------------------------
! Main module responsible for system construction and eigenvalue problem for that system.
! ---------------------------------------------------------------------------------------
module modsys
use modcommons
use modutils
use modalgs
implicit none
interface
    logical function func_simple(atomA,atomB,coupling_val)
        use modcommons
        implicit none
        type(qatom) :: atomA,atomB
        complex*16 :: coupling_val
    end function func_simple
    logical function func_matrix(atomA,atomB,coupling_val)
        use modcommons
        implicit none
        type(qatom) :: atomA,atomB
        complex*16  :: coupling_val(:,:)
    end function func_matrix
endinterface
private




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
    logical :: bOverlapMatrixEnabled

    ! other variables
    type(qatom)       :: qatom      ! auxiliary variable, can be used by user
    type(nnb_params)  :: qnnbparam  ! auxiliary variable, can be used by user

    contains
    procedure, public, pass(sys) :: init
    procedure, public, pass(sys) :: destroy
    procedure, public, pass(sys) :: add_atom!(sys,site)
    procedure, public, pass(sys) :: make_lattice!(sys,connect_procedure)
    procedure, public, pass(sys) :: update_lattice!(sys,c_default,c_simple,c_matrix)
    procedure, public, pass(sys) :: update_overlaps!(sys,o_default,o_simple,o_matrix)
    procedure, public, pass(sys) :: save_lattice!(filename,innerA,innerB)
    procedure, public, pass(sys) :: save_data!(filename,array2d,array1d,ofunit)
    procedure, public, pass(sys) :: calc_eigenproblem!(sys,pEmin,pEmax,NoStates,no_feast_contours,print_info,pmaks_iter)
    procedure, public, pass(sys) :: calc_linsys!(sys,dvec,zvec,calc_step,pardiso_mtype)


endtype qsys
public :: qsys
public :: convert_to_HB , sort_col_vals , solve_SSOLEQ
contains

! ------------------------------------------------------------------------
! Initialize system. No arguments are required.
! ------------------------------------------------------------------------
subroutine init(sys)
    class(qsys)     :: sys
    print*,"SYS::Initialization of the system."
    sys%no_atoms         = 0

    sys%qnnbparam%box     = 0.0D0
    sys%qnnbparam%distance= 0.0D0
    sys%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_BOX
    sys%bOverlapMatrixEnabled = .false.
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
    type(qatom)           :: site
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
! c_default   - user provided procedure which check if atomA in spin
!                     state s1 has coupling with atomB in spin state s2.
!                     If not function returns false. If yes function returns
!                     true and override the coupling_val argument with
!                     propper value of hopin energy.
! nnbparams         - nnbparams structure which contains some information
!                     how to search for the nearest neightbours atoms.
! ------------------------------------------------------------------------
subroutine make_lattice(sys,nnbparams,c_simple,c_matrix,o_simple,o_matrix)
    class(qsys) :: sys

    procedure(func_simple)  ,optional :: c_simple   ! }
    procedure(func_matrix)  ,optional :: c_matrix   ! } coupling functions


    procedure(func_simple)  ,optional :: o_simple   ! } overlap  functions
    procedure(func_matrix)  ,optional :: o_matrix   ! }

    type(nnb_params) :: nnbparams

    ! local variables
    integer         :: i,j,k,l,p,no_active_atoms,ns1,ns2,s1,s2
    doubleprecision :: time_start
    complex*16      :: cpl_value,cpl_1x1(1,1),cpl_delta
    complex*16,allocatable :: cpl_matrix(:,:)
    ! bounding box parameters
    integer,parameter :: ix = 1 , iy = 2 , iz = 3 , cmin = 1 , cmax =  2
    double precision  :: bbox(cmin:cmax,ix:iz)
    ! verlet box searching method
    integer              :: verlet_dims(3) , verlet_max_atoms_in_cell , ivpos(3) , vc
    integer              :: vx,vy,vz,vp
    integer, allocatable :: verlet_box(:,:,:,:) , verlet_counter(:,:,:)

    if(.not. present(c_simple)  .and. &
       .not. present(c_matrix)) then
        print*,"SYS::ERROR::No connection function has been provided to make_lattice()."
        print*,"            Cannont creat logical connection between atoms."
        stop -1
    endif
    if(nnbparams%NNB_FILTER == 0)then
        print*,"SYS::ERROR:: nnbparams%NNB_FILTER == 0 no filter specified!"
        stop -1
    endif
    sys%bOverlapMatrixEnabled = .false.
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
            do j = 1 , sys%atoms(i)%no_bonds
                call sys%atoms(i)%bonds(j)%destroy_bond()
            enddo
            if(allocated(sys%atoms(i)%bonds))     deallocate(sys%atoms(i)%bonds)
            sys%atoms(i)%no_bonds = 0
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
    if(nnbparams%NNB_FILTER == QSYS_NNB_FILTER_CHECK_ALL) then
    ! if all posible connections are tested enter this node
    do i = 1 , sys%no_atoms
    if(sys%atoms(i)%bActive == .true.) then

        ! -----------------------------------------------
        ! C_SIMPLE
        ! -----------------------------------------------
        if(present(c_simple)) then
        do j = 1 , sys%no_atoms
            if(j<i .and. QSYS_FORCE_HERMITIAN_MATRIX) cycle
            ! both sites have to be active
            if(sys%atoms(j)%bActive == .true.) then
                ! if they are nnb one may create a bond between them
                if(c_simple(sys%atoms(i),sys%atoms(j),cpl_value)) then
                    cpl_1x1 = cpl_value
                    call sys%atoms(i)%add_bond(j,cpl_1x1)

                    if(i/=j .and. QSYS_FORCE_HERMITIAN_MATRIX) then
                        call sys%atoms(j)%add_bond(i,conjg(cpl_1x1))
                        k = sys%atoms(i)%no_bonds
                        p = sys%atoms(j)%no_bonds
                        ! generate ping pong reference
                        sys%atoms(i)%bonds(k)%fromBondID = p
                        sys%atoms(j)%bonds(p)%fromBondID = k
                    endif

                else if(i==j) then ! Force diagonal term for transport
                    cpl_value = cmplx(0.0D0,0.0D0)
                    cpl_1x1   = cpl_value
                    call sys%atoms(i)%add_bond(j,cpl_1x1)
                endif
            endif
        enddo
        endif ! end of if c_simple
        ! -----------------------------------------------
        ! C_MATRIX
        ! -----------------------------------------------
        if(present(c_matrix)) then


        do j = 1 , sys%no_atoms
            if(j<i .and. QSYS_FORCE_HERMITIAN_MATRIX) cycle

            if(sys%atoms(j)%bActive == .true.) then
            if(.not. allocated(cpl_matrix)) allocate(cpl_matrix(sys%atoms(i)%no_in_states,sys%atoms(j)%no_in_states))
            if(size(cpl_matrix,1) /= sys%atoms(i)%no_in_states .or. &
               size(cpl_matrix,2) /= sys%atoms(j)%no_in_states) then
               deallocate(cpl_matrix)
               allocate(cpl_matrix(sys%atoms(i)%no_in_states,sys%atoms(j)%no_in_states))
            endif
            ! if they are nnb one may create a bond between them
            if(c_matrix(sys%atoms(i),sys%atoms(j),cpl_matrix)) then
                ! both sites have to be active
                call sys%atoms(i)%add_bond(j,cpl_matrix)

                if(i/=j .and. QSYS_FORCE_HERMITIAN_MATRIX) then
                    call sys%atoms(j)%add_bond(i,transpose(conjg(cpl_matrix)))
                    k = sys%atoms(i)%no_bonds
                    p = sys%atoms(j)%no_bonds
                    ! generate ping pong reference
                    sys%atoms(i)%bonds(k)%fromBondID = p
                    sys%atoms(j)%bonds(p)%fromBondID = k
                endif

            else if(i==j) then ! Force diagonal term for transport
                cpl_matrix = 0
                call sys%atoms(i)%add_bond(j,cpl_matrix)
            endif ! end of if nnb
            endif ! end of if j active
        enddo ! end of j

        endif ! end of if c_matrix

    endif ! end of i bactive i
    enddo ! end of i
    ! ---------------------------------------------------------------------
    ! Else use verlet method to find nnb's
    ! ---------------------------------------------------------------------
    else

    if(nnbparams%NNB_FILTER == QSYS_NNB_FILTER_DISTANCE) nnbparams%box = 2*nnbparams%distance

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
        verlet_dims(i) = NINT(abs(bbox(cmax,i) - bbox(cmin,i))/abs(nnbparams%box(i)))+1
        if( abs(nnbparams%box(i)) < 1.0D-20 ) verlet_dims(i) = 1
    enddo


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
               "            or change box parameters."
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

        ! -----------------------------------------------
        ! C_SIMPLE
        ! -----------------------------------------------
        if(present(c_simple)) then
        ! search in nearest verlet cells
        do vx = max(ivpos(1)-1,1),min(ivpos(1)+1,verlet_dims(1))
        do vy = max(ivpos(2)-1,1),min(ivpos(2)+1,verlet_dims(2))
        do vz = max(ivpos(3)-1,1),min(ivpos(3)+1,verlet_dims(3))
            ! loop over nnb verlet boxes
            do vp = 1 , verlet_counter(vx,vy,vz)
                j        = verlet_box(vx,vy,vz,vp)
                if(j<i .and. QSYS_FORCE_HERMITIAN_MATRIX) cycle


                if(sys%atoms(j)%bActive == .true.) then
                ! loop around spin states
                if(nnbparams%NNB_FILTER == QSYS_NNB_FILTER_BOX) then
                    ! if they are nnb one may create a qbond between them
                    if(c_simple(sys%atoms(i),sys%atoms(j),cpl_value)) then
                        cpl_1x1 = cpl_value
                        call sys%atoms(i)%add_bond(j,cpl_1x1)

                        if(i/=j .and. QSYS_FORCE_HERMITIAN_MATRIX) then
                            call sys%atoms(j)%add_bond(i,conjg(cpl_1x1))
                            k = sys%atoms(i)%no_bonds
                            p = sys%atoms(j)%no_bonds
                            ! generate ping pong reference
                            sys%atoms(i)%bonds(k)%fromBondID = p
                            sys%atoms(j)%bonds(p)%fromBondID = k
                        endif
                    else if(i==j ) then ! Force diagonal term for transport
                         cpl_value = cmplx(0.0D0,0.0D0)
                         cpl_1x1   = cpl_value
                         call sys%atoms(i)%add_bond(j,cpl_1x1)
                    endif
                else if(nnbparams%NNB_FILTER == QSYS_NNB_FILTER_DISTANCE) then
                    ! check distance before asking
                    if( sqrt(sum((sys%atoms(i)%atom_pos-sys%atoms(j)%atom_pos)**2)) < nnbparams%distance) then
                    ! if they are nnb one may create a qbond between them
                    if(c_simple(sys%atoms(i),sys%atoms(j),cpl_value)) then
                        cpl_1x1 = cpl_value
                        call sys%atoms(i)%add_bond(j,cpl_1x1)

                        if(i/=j .and. QSYS_FORCE_HERMITIAN_MATRIX) then
                            call sys%atoms(j)%add_bond(i,conjg(cpl_1x1))
                            k = sys%atoms(i)%no_bonds
                            p = sys%atoms(j)%no_bonds
                            ! generate ping pong reference
                            sys%atoms(i)%bonds(k)%fromBondID = p
                            sys%atoms(j)%bonds(p)%fromBondID = k
                        endif
                    else if(i==j) then ! Force diagonal term for transport
                         cpl_value = cmplx(0.0D0,0.0D0)
                         cpl_1x1   = cpl_value
                         call sys%atoms(i)%add_bond(j,cpl_1x1)
                    endif
                    endif ! end of if check distance
                endif
                endif! end of if active atom
            enddo
        enddo ! }
        enddo ! } nearest verlet cells
        enddo ! }

        endif ! end of if c_simple

        ! -----------------------------------------------
        ! C_MATRIX
        ! -----------------------------------------------
        if(present(c_matrix)) then



        ! search in nearest verlet cells
        do vx = max(ivpos(1)-1,1),min(ivpos(1)+1,verlet_dims(1))
        do vy = max(ivpos(2)-1,1),min(ivpos(2)+1,verlet_dims(2))
        do vz = max(ivpos(3)-1,1),min(ivpos(3)+1,verlet_dims(3))
            ! loop over nnb verlet boxes
            do vp = 1 , verlet_counter(vx,vy,vz)
                j        = verlet_box(vx,vy,vz,vp)
                if(j<i .and. QSYS_FORCE_HERMITIAN_MATRIX) cycle

                if(sys%atoms(j)%bActive == .true.) then
                ! loop around spin states

                if(.not. allocated(cpl_matrix)) allocate(cpl_matrix(sys%atoms(i)%no_in_states,sys%atoms(j)%no_in_states))
                if(size(cpl_matrix,1) /= sys%atoms(i)%no_in_states .or. &
                   size(cpl_matrix,2) /= sys%atoms(j)%no_in_states) then
                   deallocate(cpl_matrix)
                   allocate(cpl_matrix(sys%atoms(i)%no_in_states,sys%atoms(j)%no_in_states))
                endif

                if(nnbparams%NNB_FILTER == QSYS_NNB_FILTER_BOX) then
                    if(c_matrix(sys%atoms(i),sys%atoms(j),cpl_matrix)) then
                        call sys%atoms(i)%add_bond(j,cpl_matrix)
                        ! create hermitian system automatically
                        if(i/=j .and. QSYS_FORCE_HERMITIAN_MATRIX) then
                            call sys%atoms(j)%add_bond(i,transpose(conjg(cpl_matrix)))
                            k = sys%atoms(i)%no_bonds
                            p = sys%atoms(j)%no_bonds
                            ! generate ping pong reference
                            sys%atoms(i)%bonds(k)%fromBondID = p
                            sys%atoms(j)%bonds(p)%fromBondID = k
                        endif
                    else if( i == j) then
                        cpl_matrix = 0
                        call sys%atoms(i)%add_bond(j,cpl_matrix)
                    endif ! end if is NNB

                else if(nnbparams%NNB_FILTER == QSYS_NNB_FILTER_DISTANCE) then
                    ! check distance before asking
                    if( sqrt(sum((sys%atoms(i)%atom_pos-sys%atoms(j)%atom_pos)**2)) < nnbparams%distance) then
                    if(c_matrix(sys%atoms(i),sys%atoms(j),cpl_matrix)) then
                        call sys%atoms(i)%add_bond(j,cpl_matrix)
                        ! create hermitian system automatically

                       if(i/=j .and. QSYS_FORCE_HERMITIAN_MATRIX) then
                            call sys%atoms(j)%add_bond(i,transpose(conjg(cpl_matrix)))
                            k = sys%atoms(i)%no_bonds
                            p = sys%atoms(j)%no_bonds
                            ! generate ping pong reference
                            sys%atoms(i)%bonds(k)%fromBondID = p
                            sys%atoms(j)%bonds(p)%fromBondID = k
                        endif

                    else if( i == j) then
                        cpl_matrix = 0
                        call sys%atoms(i)%add_bond(j,cpl_matrix)
                    endif ! end if is NNB
                    endif ! end of if check distance
                endif ! end if DISTANCE filter
                endif! end of if active atom
            enddo
        enddo ! }
        enddo ! } nearest verlet cells
        enddo ! }



        endif ! end of if c_matrix
        ! -----------------------------------------------
        ! end of C_MATRIX
        ! -----------------------------------------------

        endif ! end of if active atom
    enddo ! end of i
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
!                cpl_value = sys%atoms(i)%bonds(j)%bondValue ! take hoping value from A to B
!                cpl_value_inverse = 0
                do k = 1 , sys%atoms(vp)%no_bonds ! search for the same
                    if( sys%atoms(vp)%bonds(k)%toAtomID == i )  then
!                     cpl_value_inverse = sys%atoms(vp)%bonds(k)%bondValue
                     exit
                    endif
                enddo
                ns1 = size(sys%atoms(i)%bonds(j)%bondMatrix,1)
                ns2 = size(sys%atoms(i)%bonds(j)%bondMatrix,2)
                cpl_delta = 0.0
                do s1 = 1 , ns1
                do s2 = 1 , ns2
                    cpl_delta = cpl_delta + abs( sys%atoms(i)%bonds(j)%bondMatrix(s1,s2) - conjg(sys%atoms(vp)%bonds(k)%bondMatrix(s2,s1)) )
                enddo
                enddo
                if( abs(cpl_delta) > 10.0D-6 ) then
                    print*,"==============================================================================="
                    print*,"SYS::ERROR::Created matrix will be not Hermitian. Check the provided ", &
                           "            connection function. The qbonding param between atom A->B", &
                           "            has to complex conjugate for connection between B->A."
                    print"(A,3e16.4)","             An error occured for atom A at position :",sys%atoms(i)%atom_pos
                    print"(A,3e16.4)","             and atom B at position                  :",sys%atoms(vp)%atom_pos

                    print"(A)",       "             Hoping matrix for A->B:"
                    do s1 = 1, ns1
                    print"(A,3000e16.6)","             ",sys%atoms(i)%bonds(j)%bondMatrix(s1,:)
                    enddo
                    print"(A)",       "             Hoping matrix for B->A:"
                    do s2 = 1, ns2
                    print"(A,3000e16.6)","             ",sys%atoms(vp)%bonds(k)%bondMatrix(:,s2)
                    enddo
                    print*,"==============================================================================="
                    stop -1
                endif

            enddo
        endif ! end of active atom
    enddo
    endif ! end of if disable hermiticity check

    if(present(c_matrix))deallocate(cpl_matrix)
    call sys%update_overlaps(o_simple,o_matrix)

    print*,"SYS::Connections has been found in ", get_clock()-time_start , " sec."

end subroutine make_lattice


subroutine update_lattice(sys,c_simple,c_matrix,o_simple,o_matrix)
    class(qsys) :: sys

    procedure(func_simple)  ,optional :: c_simple
    procedure(func_matrix)  ,optional :: c_matrix

    procedure(func_simple)  ,optional :: o_simple
    procedure(func_matrix)  ,optional :: o_matrix

    ! locals
    doubleprecision :: time_start
    integer :: i,j,k,vp,b,ta,fa,ns1,ns2,s1,s2
    logical :: bondTest
    complex*16      :: cpl_value,cpl_delta,cpl_1x1(1,1)
    complex*16,allocatable :: cpl_matrix(:,:)

    time_start = get_clock()

    if(.not. present(c_simple)  .and. &
       .not. present(c_matrix)) then
        print*,"SYS::ERROR::No connection function has been provided to make_lattice()."
        print*,"            Cannont creat logical connection between atoms."
        stop -1
    endif

    do i = 1 , sys%no_atoms
        if(sys%atoms(i)%bActive == .true.)then

        do b = 1 , sys%atoms(i)%no_bonds
            fa = i
            ta = sys%atoms(i)%bonds(b)%toAtomID
            if(ta<fa .and. QSYS_FORCE_HERMITIAN_MATRIX) cycle
            cpl_value = 0.0
            if(present(c_simple)) then
                bondTest = c_simple(sys%atoms(fa),sys%atoms(ta),cpl_value)
                cpl_1x1  = cpl_value

                sys%atoms(fa)%bonds(b)%bondMatrix = cpl_1x1

                if(ta/=fa .and. QSYS_FORCE_HERMITIAN_MATRIX) then
                    k = sys%atoms(fa)%bonds(b)%fromBondID
                    sys%atoms(ta)%bonds(k)%bondMatrix = (conjg(cpl_1x1))
                endif
            else if(present(c_matrix)) then
                if(.not. allocated(cpl_matrix)) allocate(cpl_matrix(sys%atoms(fa)%no_in_states,sys%atoms(ta)%no_in_states))
                if(size(cpl_matrix,1) /= sys%atoms(fa)%no_in_states .or. &
                   size(cpl_matrix,2) /= sys%atoms(ta)%no_in_states) then
                   deallocate(cpl_matrix)
                   allocate(cpl_matrix(sys%atoms(fa)%no_in_states,sys%atoms(ta)%no_in_states))
                endif
                bondTest = c_matrix(sys%atoms(fa),sys%atoms(ta),cpl_matrix)
                sys%atoms(fa)%bonds(b)%bondMatrix = cpl_matrix

                if(ta/=fa .and. QSYS_FORCE_HERMITIAN_MATRIX) then
                    k = sys%atoms(fa)%bonds(b)%fromBondID
                    sys%atoms(ta)%bonds(k)%bondMatrix = transpose(conjg(cpl_matrix))
                endif
            endif
        enddo
        endif ! end of if active atome i
    enddo


    ! Perform a check to test the hermiticity of the hamiltonian matrix.
    ! For a given atom A
    ! Take one bond to atom B and check if the value of the hoping parameter
    ! is the same as for the atom B connected with atom A.
    ! If both atoms are in the same spin state.
    if( .not. QSYS_DISABLE_HERMICITY_CHECK) then
    if(QSYS_DEBUG_LEVEL > 0) then
    print*,"SYS::INFO::Checking hermiticity of the matrix..."
    endif
    do i = 1 , sys%no_atoms ! take the atom A
        if(sys%atoms(i)%bActive) then
            do j = 1 , sys%atoms(i)%no_bonds ! take one bond to atom B
                vp = sys%atoms(i)%bonds(j)%toAtomID
!                cpl_value = sys%atoms(i)%bonds(j)%bondValue ! take hoping value from A to B
!                cpl_value_inverse = 0
                do k = 1 , sys%atoms(vp)%no_bonds ! search for the same
                    if( sys%atoms(vp)%bonds(k)%toAtomID == i )  then
!                     cpl_value_inverse = sys%atoms(vp)%bonds(k)%bondValue
                     exit
                    endif
                enddo

                ns1 = size(sys%atoms(i)%bonds(j)%bondMatrix,1)
                ns2 = size(sys%atoms(i)%bonds(j)%bondMatrix,2)
!                print*,"a->b:",i, sys%atoms(i)%bonds(j)%bondMatrix
!                print*,"b->a:",vp, sys%atoms(vp)%bonds(k)%bondMatrix
                cpl_delta = 0.0
                do s1 = 1 , ns1
                do s2 = 1 , ns2
                    cpl_delta = cpl_delta + abs( sys%atoms(i)%bonds(j)%bondMatrix(s1,s2) - conjg(sys%atoms(vp)%bonds(k)%bondMatrix(s2,s1)) )
                enddo
                enddo
                if( abs(cpl_delta) > 10.0D-6 ) then
                    print*,"==============================================================================="
                    print*,"SYS::ERROR::Created matrix will be not Hermitian. Check the provided ", &
                           "            connection function. The qbonding param between atom A->B", &
                           "            has to complex conjugate for connection between B->A."
                    print"(A,3e16.4)","             An error occured for atom A at position :",sys%atoms(i)%atom_pos
                    print"(A,3e16.4)","             and atom B at position                  :",sys%atoms(vp)%atom_pos

                    print"(A)",       "             Hoping matrix for A->B:"
                    do s1 = 1, ns1
                    print"(A,3000e16.6)","             ",sys%atoms(i)%bonds(j)%bondMatrix(s1,:)
                    enddo
                    print"(A)",       "             Hoping matrix for B->A:"
                    do s2 = 1, ns2
                    print"(A,3000e16.6)","             ",sys%atoms(vp)%bonds(k)%bondMatrix(:,s2)
                    enddo
                    print*,"==============================================================================="
                    stop -1
                endif

            enddo
        endif ! end of active atom
    enddo
    endif ! end of if disable hermiticity check

    call sys%update_overlaps(o_simple,o_matrix)
    if(QSYS_DEBUG_LEVEL > 0) then
    print*,"SYS::Lattice update done in", get_clock()-time_start , " sec."
    endif
    if(allocated(cpl_matrix)) deallocate(cpl_matrix)
endsubroutine update_lattice



subroutine update_overlaps(sys,o_simple,o_matrix)
    class(qsys) :: sys

    procedure(func_simple)  ,optional :: o_simple
    procedure(func_matrix)  ,optional :: o_matrix
    ! locals
    doubleprecision :: time_start
    integer :: i,j,k,vp,b,ta,fa,s1,s2,ns1,ns2
    logical :: bondTest
    complex*16      :: cpl_value,cpl_delta,cpl_1x1(1,1)
    complex*16,allocatable :: cpl_matrix(:,:)

    time_start = get_clock()


    do i = 1 , sys%no_atoms
        if(sys%atoms(i)%bActive == .true.)then
        ns1 = sys%atoms(i)%no_in_states
        do b = 1 , sys%atoms(i)%no_bonds
            fa = i
            ta = sys%atoms(i)%bonds(b)%toAtomID
            if(ta<fa .and. QSYS_FORCE_HERMITIAN_MATRIX) cycle

            sys%atoms(fa)%bonds(b)%overlapMatrix = 0.0
            if(present(o_simple)) then
                if(o_simple(sys%atoms(fa),sys%atoms(ta),cpl_value)) then
                    cpl_1x1 = cpl_value
                    sys%atoms(fa)%bonds(b)%overlapMatrix = cpl_1x1
                    ! update lower part of the matrix
                    if(ta/=fa .and. QSYS_FORCE_HERMITIAN_MATRIX) then
                        k = sys%atoms(fa)%bonds(b)%fromBondID
                        sys%atoms(ta)%bonds(k)%overlapMatrix = (conjg(cpl_1x1))
                    endif
                endif
            else if(present(o_matrix)) then
                if(.not. allocated(cpl_matrix)) allocate(cpl_matrix(sys%atoms(fa)%no_in_states,sys%atoms(ta)%no_in_states))
                if(size(cpl_matrix,1) /= sys%atoms(fa)%no_in_states .or. &
                   size(cpl_matrix,2) /= sys%atoms(ta)%no_in_states) then
                   deallocate(cpl_matrix)
                   allocate(cpl_matrix(sys%atoms(fa)%no_in_states,sys%atoms(ta)%no_in_states))
                endif
                if(o_matrix(sys%atoms(fa),sys%atoms(ta),cpl_matrix)) then
                    sys%atoms(fa)%bonds(b)%overlapMatrix = cpl_matrix
                    if(ta/=fa .and. QSYS_FORCE_HERMITIAN_MATRIX) then
                        k = sys%atoms(fa)%bonds(b)%fromBondID
                        sys%atoms(ta)%bonds(k)%overlapMatrix = transpose(conjg(cpl_matrix))
                    endif
                endif
            else if( fa == ta ) then
                do s1 = 1 , ns1
                sys%atoms(fa)%bonds(b)%overlapMatrix(s1,s1) = 1.0
                enddo
            endif
        enddo
        endif ! end of if active atome i
    enddo


    ! Perform a check to test the hermiticity of the hamiltonian matrix.
    ! For a given atom A
    ! Take one bond to atom B and check if the value of the hoping parameter
    ! is the same as for the atom B connected with atom A.
    ! If both atoms are in the same spin state.
    if(present(o_simple)  .or. present(o_matrix)) then
    sys%bOverlapMatrixEnabled = .true.
    if( .not. QSYS_DISABLE_HERMICITY_CHECK) then
    if(QSYS_DEBUG_LEVEL > 0)then
        print*,"SYS::INFO::Checking hermiticity of the overlap matrix."
    endif
    do i = 1 , sys%no_atoms ! take the atom A
        if(sys%atoms(i)%bActive) then
            do j = 1 , sys%atoms(i)%no_bonds ! take one bond to atom B
                vp = sys%atoms(i)%bonds(j)%toAtomID
!                cpl_value = sys%atoms(i)%bonds(j)%bondValue ! take hoping value from A to B
!                cpl_value_inverse = 0
                do k = 1 , sys%atoms(vp)%no_bonds ! search for the same
                    if( sys%atoms(vp)%bonds(k)%toAtomID == i )  then
!                     cpl_value_inverse = sys%atoms(vp)%bonds(k)%bondValue
                     exit
                    endif
                enddo

                ns1 = size(sys%atoms(i)%bonds(j)%overlapMatrix,1)
                ns2 = size(sys%atoms(i)%bonds(j)%overlapMatrix,2)
                cpl_delta = 0.0
                do s1 = 1 , ns1
                do s2 = 1 , ns2
                    cpl_delta = cpl_delta +  abs( sys%atoms(i)%bonds(j)%overlapMatrix(s1,s2) - conjg(sys%atoms(vp)%bonds(k)%overlapMatrix(s2,s1)) )
                enddo
                enddo
                if( abs(cpl_delta) > 10.0D-6 ) then
                    print*,"==============================================================================="
                    print*,"SYS::ERROR::Created matrix will be not Hermitian. Check the provided ", &
                           "            connection function. The qbonding param between atom A->B", &
                           "            has to complex conjugate for connection between B->A."
                    print"(A,3e16.4)","             An error occured for atom A at position :",sys%atoms(i)%atom_pos
                    print"(A,3e16.4)","             and atom B at position                  :",sys%atoms(vp)%atom_pos

                    print"(A)",       "             Hoping matrix for A->B:"
                    do s1 = 1, ns1
                    print"(A,3000e16.6)","             ",sys%atoms(i)%bonds(j)%overlapMatrix(s1,:)
                    enddo
                    print"(A)",       "             Hoping matrix for B->A:"
                    do s2 = 1, ns2
                    print"(A,3000e16.6)","             ",sys%atoms(vp)%bonds(k)%overlapMatrix(:,s2)
                    enddo
                    print*,"==============================================================================="
                    stop -1
                endif

            enddo
        endif ! end of active atom
    enddo
    endif ! end of if disable hermiticity check
    endif ! present overlap check
    if(QSYS_DEBUG_LEVEL > 0 .and. sys%bOverlapMatrixEnabled ) then
    print*,"SYS::Overlap matrix update done in", get_clock()-time_start , " sec."
    endif
    if(allocated(cpl_matrix)) deallocate(cpl_matrix)
endsubroutine update_overlaps


! ------------------------------------------------------------------------
! Save current lattice to file. Use this function after execution of make_lattice
! procedure.
! filename - the name of the output file
! innerA   - [optional] plot connection for atoms in state A. Default is 1
! innerB   - [optional] connected with atoms in state B. This can be used to see the
!            connections between different spin stated of atoms. Default is 1.
! ------------------------------------------------------------------------
subroutine save_lattice(sys,filename,ofunit)
    class(qsys)     :: sys
    character(*)    :: filename
    integer,optional:: ofunit

    integer :: funit
    integer :: i,b,ida,ids1,ids2,j,itmp,ns1,ns2,s1,s2
    integer,parameter :: ix = 1 , iy = 2 , iz = 3 , cmin = 1 , cmax =  2
    double precision :: lWidth , bbox(2,3)
    complex*16 :: cpl_value
    print*,"SYS::Saving lattice to file:",filename


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

    funit = 765819
    if(present(ofunit)) then
        funit = ofunit
    else
        open(unit=funit,file=filename)
    endif

    write(funit,*),"<lattice>"
    write(funit,*),"<atoms>"
    do i = 1 , sys%no_atoms

        if(sys%atoms(i)%bActive) then
            itmp = 1;
        else
            itmp = 0;
        endif

        write(funit,"(A,3e16.6,4i10,A)"),"<d>",sys%atoms(i)%atom_pos,sys%atoms(i)%flag,itmp,sys%atoms(i)%no_in_states,sys%atoms(i)%no_bonds,"</d>"

    enddo
    write(funit,*),"</atoms>"

    write(funit,*),"<connections>"
    do i = 1 , sys%no_atoms
        if(.not. sys%atoms(i)%bActive) cycle
        do j = 1 , sys%atoms(i)%no_bonds
            ns1 = size(sys%atoms(i)%bonds(j)%bondMatrix,1)
            ns2 = size(sys%atoms(i)%bonds(j)%bondMatrix,2)
            do s1 = 1 , ns1
            do s2 = 1 , ns2
            cpl_value = sys%atoms(i)%bonds(j)%bondMatrix(s1,s2)
!            if( abs(cpl_value) > 1.0D-10 )then
            write(funit,"(A,4i10,2e16.6,A)"),"<d>",i,&
                        sys%atoms(i)%bonds(j)%toAtomID,s1,s2,&
                        DBLE(cpl_value),IMAG(cpl_value),"</d>"
!            endif
            enddo
            enddo
        enddo
    enddo
    write(funit,*),"</connections>"
    write(funit,*),"</lattice>"

    if(.not. present(ofunit)) close(funit)

end subroutine save_lattice

subroutine save_data(sys,filename,array2d,array1d,ofunit)
    class(qsys)     :: sys
    character(*)    :: filename
    integer,optional:: ofunit
    doubleprecision,optional :: array2d(:,:) , array1d(:)
    integer         :: funit,no_cols,i,j,gi,iter
    doubleprecision,allocatable, dimension(:) :: values
    character(128) :: sformat , stmp

    no_cols = 0
    if(present(array1d)) no_cols =  1
    if(present(array2d)) no_cols =  no_cols + size(array2d(:,1))
    if(no_cols == 0) then
        print*,"SYS::No data to save:",filename
        return
    endif
    allocate(values(no_cols))

    funit = 765819
    if(present(ofunit)) then
        funit = ofunit
    else
        open(unit=funit,file=filename)
    endif
    write(funit,*),"<sysdata>"
    write(funit,*),"<info>",no_cols,"</info>"
    write(funit,*),"<name>",filename,"</name>"
    write(stmp,"(i)"),no_cols
    sformat = "(A,2i,"//trim(stmp)//"e20.6,A)"

    do i = 1, sys%no_atoms
        if(.not. sys%atoms(i)%bActive) cycle
        do j = 1 , sys%atoms(i)%no_in_states
        gi = sys%atoms(i)%globalIDs(j)
        iter = 0
        if(present(array1d)) then
            iter = iter + 1
            values(iter) = array1d(gi)
        endif
        if(present(array2d)) then
            iter = iter + 1
            values(iter:no_cols) = array2d(:,gi)
        endif

        write(funit,sformat),"<d>",i,j,values,"</d>"
        enddo
    enddo
    write(funit,*),"</sysdata>"
    if(.not. present(ofunit)) close(funit)
    deallocate(values)
end subroutine save_data

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
    integer :: i,j,info,itmp,nw,M0,loop,no_evals,ta,ts,ns1,ns2,s1,s2
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

    if(sys%bOverlapMatrixEnabled) then
        print*,"==============================================================================="
        print*,"SYS::ERROR::Eigenvalue problem not supported for system with "
        print*,"            Overlap matrix different than identity matrix."
        print*,"==============================================================================="
        stop -1
    endif
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
            do j = 1, sys%atoms(i)%no_bonds
            ns1 = size(sys%atoms(i)%bonds(j)%bondMatrix,1)
            ns2 = size(sys%atoms(i)%bonds(j)%bondMatrix,2)
!            if( abs(sys%atoms(i)%bonds(j)%bondValue) > QSYS_COUPLING_CUTOFF ) &
            itmp = itmp + ns1*ns2
            enddo
        endif
    enddo
    NO_NON_ZERO_VALUES = itmp
    ! The number of unknows is taken from the last global index
    NO_VARIABLES       = sys%system_size
    if(QSYS_DEBUG_LEVEL > 0) then
    print*,"SYS::Calulating eigenvalue problem using FEAST solver"
    endif
    allocate(MATHVALS(NO_NON_ZERO_VALUES))
    allocate(ROWCOLID(NO_NON_ZERO_VALUES,2))

    ! Filling matrix and row-col array
    itmp = 0
    do i = 1 ,  sys%no_atoms

        if(sys%atoms(i)%bActive) then

        ns1   = sys%atoms(i)%no_in_states
        do s1 = 1 , ns1

        do j  = 1, sys%atoms(i)%no_bonds
        ns2 = size(sys%atoms(i)%bonds(j)%bondMatrix,2)


        do s2 = 1 , ns2
            itmp = itmp + 1
            MATHVALS(itmp)   = sys%atoms(i)%bonds(j)%bondMatrix(s1,s2)
            ROWCOLID(itmp,1) = sys%atoms(i)%globalIDs(s1)
            ta = sys%atoms(i)%bonds(j)%toAtomID
            ROWCOLID(itmp,2) = sys%atoms(ta)%globalIDs(s2)
        enddo

        ! remove zero ellements
!        if( abs(MATHVALS(itmp)) < QSYS_COUPLING_CUTOFF )  itmp = itmp - 1

        enddo
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
        if(QSYS_DEBUG_LEVEL > 0) then
        print*,"SYS::Eigenvalues calculated. Found:",sys%no_eigenvalues
        endif
end subroutine calc_eigenproblem



! ------------------------------------------------------------------------
! Solve sparse system of linear equations based on created system.
! Params:
! dvec,zvec [optional] - vector b of size of number of variables: Ax=b
!           dvec - double precision vector
!           zvec - complex vector
! calc_step [def=4,optional] - divide solver on steps. 1 - factorize, 2- solve, 3-free memory, 4-all steps
! pardiso_mtype [optional] - default QSYS_LINSYS_PARDISO_REAL_NON_SYM for real numbers,
!           and QSYS_LINSYS_PARDISO_CMPLX_NON_SYM for complex numbers. This is used when PARDISO solver is used.
! ------------------------------------------------------------------------
subroutine calc_linsys(sys,dvec,zvec,calc_step,pardiso_mtype)
    class(qsys) :: sys
    doubleprecision,optional,dimension(:),intent(inout) :: dvec
    complex*16,optional,dimension(:),intent(inout)      :: zvec

    integer,optional :: calc_step, pardiso_mtype

    integer :: step,mtype

    integer,save :: NO_NON_ZERO_VALUES , NO_VARIABLES
    integer :: i,j,itmp,ns1,ns2,s1,s2,ta

    logical :: bUseComplex,bTestVecSize

    doubleprecision :: time_start,timer_factorization

    complex*16     ,dimension(:)  , allocatable,save  :: ZMATHVALS
    doubleprecision,dimension(:)  , allocatable,save  :: DMATHVALS
    integer   ,dimension(:,:), allocatable,save       :: ROWCOLID
    integer,allocatable,save                          :: HBROWS(:)

    if( .not. present(dvec) .and. .not. present(zvec) ) then
        print*,"SYS::Linear system solver cannot be initialized, vector b is not defined in the calc_linsys call"
        return
    endif

    ! Choose type of numerics
    bUseComplex = .false.
    bTestVecSize= .true.
    if(present( zvec )) bUseComplex = .true.

    step = QSYS_LINSYS_ALL_STEPS
    if(present(calc_step)) step = calc_step

    mtype = QSYS_LINSYS_PARDISO_REAL_NON_SYM
    if(bUseComplex) then
        mtype = QSYS_LINSYS_PARDISO_CMPLX_NON_SYM
    endif
    if(present(pardiso_mtype)) mtype = pardiso_mtype

    if(bUseComplex) then
        if(size(zvec) /= sys%system_size ) bTestVecSize = .false.
        ta = size(zvec)
    else
        if(size(dvec) /= sys%system_size ) bTestVecSize = .false.
        ta = size(dvec)
    endif

    if(.not. bTestVecSize) then
        print*,"SYS::The size of the vector b in linear solver is different than number of variables."
        print*,"     System cannot be solved."
        print*,"     Num. variables:",sys%system_size
        print*,"     Vector b size :",ta
        return
    endif

    time_start = get_clock()


    if(step == QSYS_LINSYS_STEP_FACTORIZE .or. step == QSYS_LINSYS_ALL_STEPS) then
    ! Calculate the number of non-zero elements in matrix
    itmp = 0
    do i = 1, sys%no_atoms
        if(sys%atoms(i)%bActive) then
            do j = 1, sys%atoms(i)%no_bonds
                ns1 = size(sys%atoms(i)%bonds(j)%bondMatrix,1)
                ns2 = size(sys%atoms(i)%bonds(j)%bondMatrix,2)
                itmp = itmp + ns1*ns2
            enddo
        endif
    enddo
    NO_NON_ZERO_VALUES = itmp
    ! The number of unknows is taken from the last global index
    NO_VARIABLES       = sys%system_size

    if(QSYS_DEBUG_LEVEL > 0) then
        print*,"SYS::Solving system of linear equations"
    endif

    allocate(ROWCOLID(NO_NON_ZERO_VALUES,2))
    if(bUseComplex) then
        allocate(ZMATHVALS(NO_NON_ZERO_VALUES))
    else
        allocate(DMATHVALS(NO_NON_ZERO_VALUES))
    endif
    ! Filling matrix and row-col array
    itmp = 0
    do i = 1 ,  sys%no_atoms
        if(sys%atoms(i)%bActive) then
        ns1   = sys%atoms(i)%no_in_states
        do s1 = 1 , ns1
        do j  = 1, sys%atoms(i)%no_bonds
        ns2 = size(sys%atoms(i)%bonds(j)%bondMatrix,2)
        do s2 = 1 , ns2
            itmp = itmp + 1
            if(bUseComplex) then
                ZMATHVALS(itmp) = sys%atoms(i)%bonds(j)%bondMatrix(s1,s2)
            else
                DMATHVALS(itmp) = DBLE(sys%atoms(i)%bonds(j)%bondMatrix(s1,s2))
            endif
            ROWCOLID(itmp,1) = sys%atoms(i)%globalIDs(s1)
            ta = sys%atoms(i)%bonds(j)%toAtomID
            ROWCOLID(itmp,2) = sys%atoms(ta)%globalIDs(s2)
        enddo ! end of do s2
        enddo ! end of do j
        enddo ! end of do s1
        endif ! end of if active atom
    enddo ! end of do i atoms


    if(QSYS_DEBUG_LEVEL>0) then
        print*,"--------------------------------------------------"
        print*,"SYS::LINEAR SYSTEM"
        print*,"--------------------------------------------------"
        print*,"SYS::Number of variables            :",NO_VARIABLES
        print*,"SYS::Total number of non-zero values:",NO_NON_ZERO_VALUES
    endif


    allocate(HBROWS(NO_VARIABLES+1))
    if(bUseComplex) then
        call convert_to_HB(NO_NON_ZERO_VALUES,ROWCOLID,ZMATHVALS,HBROWS)
    else ! other wise use double precision equivalent
        call dalg_convert2HB(NO_NON_ZERO_VALUES,ROWCOLID,DMATHVALS,HBROWS)
    endif

    timer_factorization  = get_clock()
    if(bUseComplex) then
    call solve_SSOLEQ(NO_VARIABLES, &
                      NO_NON_ZERO_VALUES,&
                      ROWCOLID(:,2),HBROWS,&
                      ZMATHVALS(:),zvec,1,mtype)

    else
    call dalg_SSOLEQ(NO_VARIABLES, &
                      NO_NON_ZERO_VALUES,&
                      ROWCOLID(:,2),HBROWS,&
                      DMATHVALS(:),dvec,1,mtype)

    endif
    timer_factorization = get_clock() - timer_factorization
    if(QSYS_DEBUG_LEVEL > 0) then
        print*,"SYS::Linear solver::Factorization time:",timer_factorization
    endif

    endif ! end of if system factorize


    if(step == QSYS_LINSYS_STEP_SOLVE .or. step == QSYS_LINSYS_ALL_STEPS) then
        timer_factorization  = get_clock()
        if(bUseComplex) then
        call solve_SSOLEQ(NO_VARIABLES, &
                          NO_NON_ZERO_VALUES,&
                          ROWCOLID(:,2),HBROWS,&
                          ZMATHVALS(:),zvec,2,mtype)

        else
        call dalg_SSOLEQ(NO_VARIABLES, &
                          NO_NON_ZERO_VALUES,&
                          ROWCOLID(:,2),HBROWS,&
                          DMATHVALS(:),dvec,2,mtype)

        endif
        timer_factorization = get_clock() - timer_factorization
        if(QSYS_DEBUG_LEVEL > 0) then
            print*,"SYS::Linear solver::Solve time:",timer_factorization
        endif
    endif ! end of if solve system


    if(step == QSYS_LINSYS_STEP_FREE_MEMORY .or. step == QSYS_LINSYS_ALL_STEPS) then
        if(bUseComplex) then
        call solve_SSOLEQ(NO_VARIABLES, &
                          NO_NON_ZERO_VALUES,&
                          ROWCOLID(:,2),HBROWS,&
                          ZMATHVALS(:),zvec,3,mtype)
        deallocate(ZMATHVALS)
        else
        call dalg_SSOLEQ(NO_VARIABLES, &
                          NO_NON_ZERO_VALUES,&
                          ROWCOLID(:,2),HBROWS,&
                          DMATHVALS(:),dvec,3,mtype)
        deallocate(DMATHVALS)
        endif
        deallocate(ROWCOLID)
        deallocate(HBROWS)
    endif ! end of if free memory

end subroutine calc_linsys

end module modsys
