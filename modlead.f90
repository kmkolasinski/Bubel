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
implicit none
private
complex*16,parameter :: II = cmplx(0.0D0,1.0D0)
type qlead
    type(qshape) :: lead_shape                         ! contains the area/volume in which are
                                                       ! located atoms which belong to the lead.
    doubleprecision :: lead_vector(3)                  ! direction of propagation of the lead in the periodic structur

    complex*16,dimension(:,:) , allocatable :: valsH0  ! Hamiltonian of the lead - unit cell
    complex*16,dimension(:,:) , allocatable :: valsTau ! Coupling matrix to the next unit cell

    integer :: no_sites                                ! number of unknowns (number of sites in the lead,
                                                       ! including spin states)

    integer,dimension(:,:),allocatable        ::      l2g ! mapping from local id in lead to global id (:,1)  = atom_ID , (:,2) = spin_ID
    integer,dimension(:,:),allocatable        :: next_l2g ! the same as above but mapping to the atoms in the next unit cell

    contains
    procedure,pass(this) :: init_lead!()
    procedure,pass(this) :: bands!(this,filename,kmin,kmax,dk,Emin,Emax)
    procedure,pass(this) :: destroy!()
    procedure,pass(this) :: print_lead!(this,filename,all_atoms)

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
    integer :: i,j,b,b_child
    integer :: no_sites,system_size, atom_id,bond_id,bond_atom_id,bond_id_child,bond_atom_id_child
    integer :: irow,icol,aid,gid,sid,lid
    integer :: no_atoms
    complex :: cval
    integer,allocatable :: next_cell_atoms(:)
    doubleprecision,dimension(3) :: cellA_pos,cellB_pos,tmp_pos,cellBA_vec

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



    no_atoms = 0
    do i = 1 , size(all_atoms) ! Search in all active atoms
        if(all_atoms(i)%bActive) then
            ! if atom is in the unit cell
            tmp_pos = all_atoms(i)%atom_pos - lvec
            if( lshape%is_inside(tmp_pos) == .true. ) then
                ! search for atom with the same position
                do j = 1 , no_sites
                    cellA_pos = all_atoms(this%l2g(j,1))%atom_pos
                    if( sqrt(sum( tmp_pos - cellA_pos )**2) < 1.0E-10 ) then
                        exit
                    endif
                enddo
                if( j > no_sites ) then
                    print"(A)",              " SYS::LEAD::ERROR::The translation"
                    print"(A,i9,A,3e10.2,A)","                   from atom:",i," at position r=(",all_atoms(i)%atom_pos,")"
                    print"(A,i9,A,3e10.2,A)","                   to atom  :",i," at position r=(",all_atoms(i)%atom_pos,")"
                    print"(A)",              "                   is impossible by applying lead translation vector. Maybe some of the "
                    print"(A)",              "                   atoms have not been inlcuded in the lead. Check lead area and T vector."
                    stop -1
                endif
                ! j keeps local ID of the same atom in unit
                this%next_l2g(j,1) = i
                this%next_l2g(j,2) = this%l2g(j,2)

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
                    if( sqrt(sum( cellB_pos - cellA_pos )**2) < 1.0E-10 ) then
                        exit
                    endif
                enddo
                if( j > no_sites ) then
                    print"(A)",              " SYS::LEAD::ERROR::The translation"
                    print"(A,i9,A,3e10.2,A)","                   from atom:",atom_id," at position r=(",all_atoms(atom_id)%atom_pos,")"
                    print"(A,i9,A,3e10.2,A)","                   to atom  :",bond_atom_id," at position r=(",all_atoms(bond_atom_id)%atom_pos,")"
                    print"(A)",              "                   is impossible by applying lead translation vector. Maybe some of the "
                    print"(A)",              "                   atoms have not been inlcuded in the lead. Check lead area and T vector."
                    stop -1
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


            else ! in case of coupling occures within lead

                aid = bond_atom_id ! id of connected atom
                sid = all_atoms(atom_id)%bonds(b)%toInnerID ! and its spin ID
                gid = all_atoms(aid)%globalIDs(sid) ! convert to global state id (atom+spin)
                lid = tmp_g2l(gid) ! remap to local ID in lead
                irow = i
                icol = lid
                cval = all_atoms(atom_id)%bonds(b)%bondValue
                this%valsH0(irow,icol) = cval
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
    if(allocated(this%valsH0))  deallocate(this%valsH0)
    if(allocated(this%valsTau)) deallocate(this%valsTau)
    if(allocated(this%l2g))     deallocate(this%l2g)
    if(allocated(this%next_l2g))deallocate(this%next_l2g)
    this%no_sites = 0
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

subroutine print_lead(this,filename,all_atoms)
    class(qlead) :: this
    character(*) :: filename
    type(qatom),dimension(:) :: all_atoms

    integer,parameter :: funit = 5437629
    integer         :: i,j,id_atom_a,id_spin_a,id_atom_b,id_spin_b
    doubleprecision :: max_abs_matrix_element,normalized_value
    doubleprecision :: CUTOFF_LEVEL,weight
    doubleprecision :: tmp_pos(3)

    CUTOFF_LEVEL = 1.0D-10
    print*,"SYS::LEAD::Writing lead data to file:",filename
    open(unit=funit,file=filename)
    write(funit,"(A)"),"<lead>"
    call this%lead_shape%flush_shape_data_to_file(funit)
    write(funit,"(A)"),"<lead_vector>"
    write(funit,"(3e20.6)"),this%lead_vector
    write(funit,"(A)"),"</lead_vector>"

    ! ----------------------------------------------------------------------------
    ! LEAD DATA
    ! ----------------------------------------------------------------------------
    write(funit,"(A)"),"<lead_data>"
    max_abs_matrix_element = 0
    do i = 1 , this%no_sites
        do j = i , this%no_sites
            if( abs(this%valsH0(i,j)) > max_abs_matrix_element ) max_abs_matrix_element = abs(this%valsH0(i,j))
        enddo
    enddo
    do i = 1 , this%no_sites
        id_atom_a = this%l2g(i,1)

        if( all_atoms(id_atom_a)%bActive) then

        id_spin_a = this%l2g(i,2)
        do j = i , this%no_sites
        id_atom_b = this%l2g(j,1)
        id_spin_b = this%l2g(j,2)
        normalized_value = abs(this%valsH0(i,j))/max_abs_matrix_element
        if( normalized_value > CUTOFF_LEVEL ) then
            write(funit,"(A,7e20.6,2i5,A)"),"   <data>",all_atoms(id_atom_a)%atom_pos,all_atoms(id_atom_b)%atom_pos,normalized_value,id_spin_a,id_spin_b,"</data>"
        endif
        enddo
        endif
    enddo
    write(funit,"(A)"),"</lead_data>"
    ! ----------------------------------------------------------------------------
    ! NEXT UNIT CELL DATA
    ! ----------------------------------------------------------------------------
    write(funit,"(A)"),"<next_cell_lead_data>"
    do i = 1 , this%no_sites
        id_atom_a = this%next_l2g(i,1)

        if(this%next_l2g(i,1) == 0) then
            exit
        endif
        if(all_atoms(id_atom_a)%bActive) then
        id_spin_a = this%next_l2g(i,2)
        do j = i , this%no_sites
        id_atom_b = this%next_l2g(j,1)
        id_spin_b = this%next_l2g(j,2)
        normalized_value = abs(this%valsH0(i,j))/max_abs_matrix_element
        if( normalized_value > CUTOFF_LEVEL ) then
            write(funit,"(A,7e20.6,2i5,A)"),"   <data>",all_atoms(id_atom_a)%atom_pos,all_atoms(id_atom_b)%atom_pos,normalized_value,id_spin_a,id_spin_b,"</data>"
        endif
        enddo
        endif
    enddo
    write(funit,"(A)"),"</next_cell_lead_data>"

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
        if( all_atoms(id_atom_a)%bActive) then
        id_spin_a = this%l2g(i,2)
        do j = 1 , this%no_sites
        id_atom_b = this%next_l2g(j,1)
        id_spin_b = this%next_l2g(j,2)
        normalized_value = abs(this%valsTau(i,j))/max_abs_matrix_element
        if( normalized_value > CUTOFF_LEVEL ) then
            write(funit,"(A,7e20.6,2i5,A)"),"   <data>",all_atoms(id_atom_a)%atom_pos,all_atoms(id_atom_b)%atom_pos,normalized_value,id_spin_a,id_spin_b,"</data>"
        endif
        enddo
        endif
    enddo
    write(funit,"(A)"),"</lead_coupling>"

    ! ----------------------------------------------------------------------------
    ! SET OF ATOMS NEAR TO LEAD
    ! ----------------------------------------------------------------------------
    write(funit,"(A)"),"<nearest_atoms>"

    do i = 1 , size(all_atoms)
        if(all_atoms(i)%bActive) then

        do j = -5 , 5 ! check five nearest unit cells
            tmp_pos = all_atoms(i)%atom_pos - (j-1)*this%lead_vector
            if(this%lead_shape%is_inside(tmp_pos)) then
                weight = 1-(abs(j)-1)/(5.0-1.0)
                if(weight >= 0.5) weight = 1
                write(funit,"(A,4e20.6,A)"),"   <data>",all_atoms(i)%atom_pos,weight,"</data>"
                exit
            endif
        enddo

        endif

    enddo
    write(funit,"(A)"),"</nearest_atoms>"


    write(funit,"(A)"),"</lead>"
    close(funit)

endsubroutine print_lead

endmodule modlead
