! ------------------------------------------------------------ !
! Finding selected elements of the inverse of sparce matrix.
! Implementation based on the Ref. "Knitting algorithm for calculating
! Green functions in quantum systems".
! * K. Kazymyrenko and X. Waintal, Phys. Rev. B 77, 115119
! How to use:
! 1. Call selinv%init_memory(no_nodes,no_vals,vals,rc)
! 2. and find inverses for selected elements: call selinv%invert(elements,evals)
! Bugs:
! - does not work for some matrices
! - for now cannot be used to find first element of the matrix
! ------------------------------------------------------------ !
module modskminv
use modcommons
implicit none
private
integer :: SKMINV_DEBUG_LEVEL = 0 ! show 1, hide 0 debug information

! -----------------------------------------------------------
! Structure which holds the value of the n-th element in the
! matrix and the coupling values to neightbours. This type
! is generated from the rows(),cols(),values() information of
! the matrix.
! -----------------------------------------------------------
type node
    complex*16,allocatable :: Tij(:) ! hoping from current i-th element to neightboors T(i->j)
    complex*16,allocatable :: Tji(:) ! inverse of above T(j->i)
    integer   ,allocatable :: Iij(:) ! global indices of neightbours
    integer   ,allocatable :: Iij_less(:) ! global indices of neightbours with id smaller than current node
    complex*16 :: H0                 ! on site value
    integer    :: no_bonds           ! number of hopings, size of the Tij,Tji, and Iij vectors
    integer    :: no_used            ! number of connections in dynamic interface (used in algorithm)
    integer    :: fl_id              ! dynamically contains mapping from global id to local id in the interface (here called frontline)
    integer    :: no_bonds_less      ! number of neightbours with id smaller than current node id. Size of Iij_less
    contains
    procedure,pass(this) :: reset_node ! reset variables and deallocate arrays
    procedure,pass(this) :: alloc_mem  ! allocate all arrays
end type node

! -----------------------------------------------------------
! Describes the moving interface of not fully connected nodes
! Here we use one structure to define the moving frontline and
! Greens function of selected elements.
! -----------------------------------------------------------
type frontline
    complex*16,allocatable :: Gij(:,:) ! Greens function of the moving interface/frontline
    integer   ,allocatable :: Iij(:)   ! global id-s of nodes in the interface/frontline
    integer    :: no_sites             ! number of sites in the frontline
!    complex*16 :: Gp,Gxa,Gax           ! Gp - holds greens function at specific point (i,j), Gxa,Gax - are temporal variables
!    integer    :: i,j                  ! defines the element to be inversed
    contains
    procedure,pass(this) :: free_frontline!free memory
    procedure,pass(this) :: init_frontline!allocate memory
!    procedure,pass(this) :: init_point!allocate memory for element (i,j)
    procedure,pass(this) :: copy_frontline!(this,scr_ft)
end type frontline

type element
    complex*16,allocatable :: Gij(:,:) ! Greens function of the moving interface/frontline
    integer    :: no_sites             ! number of sites in the frontline
    complex*16 :: Gp,Gxa,Gax           ! Gp - holds greens function at specific point (i,j), Gxa,Gax - are temporal variables
    integer    :: i,j                  ! defines the element to be inversed
    contains
    procedure,pass(this) :: free_element!free memory
    procedure,pass(this) :: init_element!allocate memory for element (i,j)
    procedure,pass(this) :: copy_element!(this,scr_ft)
end type element

! -----------------------------------------------------------
! Main type responsible for the selected inversion of the
! sparse matrix. Contains structure of the matrix represented
! by nodes and seleted elements to be inversed.
! -----------------------------------------------------------
type skminv
    type(node) , allocatable    :: nodes(:)            ! node representation of the inversed matrix
    type(frontline)             :: flNew,flOld         ! dynamic frontline
    type(element),allocatable :: elOld(:),elNew(:)   ! dynamic Greens functions of selected elements
    integer                     :: no_nodes,no_elements! size of the matrix, and number of points to be inversed
    contains
    procedure,pass(this) :: init_memory
    procedure,pass(this) :: invert
    procedure,pass(this) :: free_memory
end type skminv


! -----------------------------------------------------------
! Make public some of the variables
! -----------------------------------------------------------
public :: SKMINV_DEBUG_LEVEL
public :: skminv
contains

! -----------------------------------------------------------
!                           SKMINV
! Free all allocated memory: matrix and selected points.
! -----------------------------------------------------------
subroutine free_memory(this)
    class(skminv) :: this
    integer :: i
    ! deallocate matrix
    do i = 1 , this%no_nodes
        call this%nodes(i)%reset_node()
    enddo
    if(allocated(this%nodes))deallocate(this%nodes)
    ! just in case deallocate selected points, they
    ! should be deallocated in invert function
    do i = 1 , this%no_elements
        call this%elNew(i)%free_element()
        call this%elOld(i)%free_element()
    enddo
    if(allocated(this%elNew))deallocate(this%elNew)
    if(allocated(this%elOld))deallocate(this%elOld)

    ! free moving frontlines
    call this%flNew%free_frontline()
    call this%flOld%free_frontline()
    this%no_nodes    = 0
    this%no_elements = 0
end subroutine free_memory

! -----------------------------------------------------------
!                           SKMINV
! Find the selected elements of the inverse of the input matrix.
! elements(no_elements,2) - array which contains pairs (row,col)
!               to elements which are suppose to be inverted
! values(no_elements) - contains values of selected elements.
!
! TODO:
! call init_memory before this function to allocate memory.
! -----------------------------------------------------------
subroutine invert(this,elements,values)
    class(skminv)           :: this
    integer,dimension(:,:)  :: elements
    complex*16,dimension(:) :: values
    ! internal variables:
    integer    :: no_elems,i
    complex*16 :: GAA,GXA,GAX,tmpC
    integer    :: c_node,c_elem,no_fl,k,i1,i2,s1,s2,p,q
    type(element) :: tmp_element
    ! profiling variables
    doubleprecision :: progress,progress_it
    doubleprecision :: t_step0,t_step1,t_step1e,t_step2,t_step2e,t_tmp,t_step3,t_step4,t_total

    complex(16),allocatable :: oGij1(:,:),nGij1(:,:)
    complex(16),allocatable :: oGij2(:,:),nGij2(:,:)
    complex(16),allocatable :: e_Gax(:),e_Gxa(:),e_Gp(:)



    if(SKMINV_DEBUG_LEVEL > 0) then
        print*,"skminv::starting calculation for:"
        print*,"        Num. elements:",size(values,1)
        print*,"        Matrix. size :",this%no_nodes
    endif
    t_step0  = 0
    t_step1  = 0
    t_step1e = 0
    t_step2  = 0
    t_step2e = 0
    t_step3  = 0
    t_step4  = 0
    t_total  = 0

    t_tmp   = get_clock()
    t_total = t_tmp
    ! allocate elements
    no_elems         = size(values,1)
    this%no_elements = no_elems
    allocate(this%elNew(no_elems))
    allocate(this%elOld(no_elems))

    allocate(oGij1(no_elems,1))
    allocate(oGij2(no_elems,1))
    allocate(nGij1(no_elems,1))
    allocate(nGij2(no_elems,1))
    allocate(e_Gax(no_elems))
    allocate(e_Gxa(no_elems))
    allocate(e_Gp(no_elems ))
    oGij1  = 0
    oGij2  = 0
    nGij1  = 0
    nGij2  = 0
    e_Gp  = 0
    e_Gxa = 0
    e_Gax = 0

    ! initialize elements to invert
    do i = 1,no_elems
        call this%elOld(i)%init_element(1)
        this%elOld(i)%i      = elements(i,1)
        this%elOld(i)%j      = elements(i,2)
        this%elOld(i)%Gp     = 0
!        this%elOld(i)%Iij(1) = 1 ! set starting element
        call this%elNew(i)%init_element(1)
        this%elNew(i)%i  = this%elOld(i)%i
        this%elNew(i)%j  = this%elOld(i)%j
        this%elNew(i)%Gp = 0
!        this%elNew(i)%Iij(1) = 1 ! set starting element
    enddo

    ! initialize first element to be inverted
    call this%flNew%free_frontline()
    call this%flOld%free_frontline()
    call this%flNew%init_frontline(1)
    call this%flOld%init_frontline(1)

    this%flNew%Gij(1,1)  = 1/this%nodes(1)%H0
    this%flNew%Iij(1)    = 1
    this%nodes(1)%fl_id  = 1
    call this%flOld%copy_frontline(this%flNew)

    t_step0 = get_clock() - t_tmp
    ! Perform recursive iteration through all the matrix
    progress_it = 0.1
    do c_node = 2 , this%no_nodes
        progress = c_node/(this%no_nodes+0.0)
        if(progress > progress_it) then
            progress_it = progress_it + 0.05
            if(SKMINV_DEBUG_LEVEL == 2) then
                print"(A,I3,A)"," skminv::solving stage:",nint(progress*100),"%"
            endif
        endif
        t_tmp = get_clock()
        ! Calculate the Gaa from Eq. (0)
        tmpC = 0
        do s1 = 1 , this%nodes(c_node)%no_bonds_less
            i1 = this%nodes(this%nodes(c_node)%Iij_less(s1))%fl_id
            do s2 = 1 , this%nodes(c_node)%no_bonds_less
                i2   = this%nodes(this%nodes(c_node)%Iij_less(s2))%fl_id
                tmpC = tmpC + this%nodes(c_node)%Tij(s1)*this%nodes(c_node)%Tji(s2)*this%flOld%Gij(i2,i1)
            enddo ! end of s2
        enddo ! end of s1

        ! Updade first step of algorithm
        GAA = 1.0/( this%nodes(c_node)%H0 - tmpC )

        ! Actualize the number of coupling on the frontline after addition
        ! of new site.
        do s1 = 1 , this%nodes(c_node)%no_bonds_less
            k = this%nodes(c_node)%Iij_less(s1)
            this%nodes(k)%no_used      = this%nodes(k)%no_used      + 1
            this%nodes(c_node)%no_used = this%nodes(c_node)%no_used + 1
        enddo

        if(SKMINV_DEBUG_LEVEL == 1) then
            ! For the last step this should be equal to last element
            print"(A,i,A,2f12.6,A)","Current node:",c_node," Gaa value=(",GAA,")"
        endif
        t_step1  = t_step1 + get_clock() - t_tmp
        t_tmp    = get_clock()
        ! -----------------------------------------------
        ! Updating the Gaa values of the selected elements
        do c_elem = 1 , no_elems
            GXA   = 0
            GAX   = 0
            ! special treatment is done if current element is equal c_node
            ! if yes Eq. (9) cannot be use and we use (8) instead
            if( this%elOld(c_elem)%j == c_node  )then
                GXA = GAA
            else ! normal treatement
                do s1 = 1 , this%nodes(c_node)%no_bonds_less
                    i2  = this%nodes(this%nodes(c_node)%Iij_less(s1))%fl_id
                    GXA = GXA - this%nodes(c_node)%Tij(s1)*this%elOld(c_elem)%Gij(i2,1)*GAA
                enddo
            endif
            ! same for i-th component
            if( this%elOld(c_elem)%i == c_node  )then
                GAX = GAA
            else
                do s1 = 1 , this%nodes(c_node)%no_bonds_less
                    i2 = this%nodes(this%nodes(c_node)%Iij_less(s1))%fl_id
                    GAX = GAX - this%nodes(c_node)%Tji(s1)*this%elOld(c_elem)%Gij(i2,2)*GAA
                enddo
            endif
            ! Update element value and temporal Gax,Gxa for later use.
            this%elOld(c_elem)%Gax = GAX
            this%elOld(c_elem)%Gxa = GXA
            this%elOld(c_elem)%Gp  = this%elOld(c_elem)%Gp + GXA*GAX/GAA


            e_Gax(c_elem) = GAX
            e_Gxa(c_elem) = GXA
            e_Gp (c_elem) = e_Gp (c_elem) + GXA*GAX/GAA
        enddo ! end of loop over elements


        ! -----------------------------------------------
        ! Updating the number of sites on the frontline
        k = 0
        do p = 1 , this%flOld%no_sites
            ! count only those sites which number of used connection is smaller than number of bounds
            if(this%nodes(this%flOld%Iij(p))%no_used < this%nodes(this%flOld%Iij(p))%no_bonds )then
                k = k + 1
            endif
        enddo
        if(this%nodes(c_node)%no_used < this%nodes(c_node)%no_bonds) k = k + 1
        no_fl = k ! new number of sites in the fronline


        if(no_fl == 0) exit ! this will exit for the last element where frontline is closed
        call this%flNew%init_frontline(no_fl) ! allocate memory for new fronline
        k = 0
        do p = 1, this%flOld%no_sites
            if(this%nodes(this%flOld%Iij(p))%no_used < this%nodes(this%flOld%Iij(p))%no_bonds )then
                k = k + 1
                this%flNew%Iij(k) = this%flOld%Iij(p)
            endif
        enddo
        ! Add new site to the frontline
        if(this%nodes(c_node)%no_used < this%nodes(c_node)%no_bonds) then
            k = k + 1
            this%flNew%Iij(k)         = c_node
        endif

        if(SKMINV_DEBUG_LEVEL == 1) then
            ! For the last step this should be equal to last element
            print"(A,i,A,10000i)","New fron line size:",no_fl," fl. ids:",this%flNew%Iij
        endif

        t_step1e = t_step1e + get_clock() - t_tmp
        t_tmp    = get_clock()
        ! Update of the rest elements of the frontline Greens function:
        ! See Eq. (9) and (10), note that we have minus sign. Which
        ! was working for us for any matrix.
        do k = 1 , no_fl-1
            p   = this%flNew%Iij(k) ! atom do ktorego nalezy nowy gj
            GXA = 0
            GAX = 0
            i1  = this%nodes(p)%fl_id
            do s1 = 1 , this%nodes(c_node)%no_bonds_less
                q = this%nodes(c_node)%Iij_less(s1)
                i2 = this%nodes(q)%fl_id
                GXA = GXA + this%nodes(c_node)%Tij(s1)*this%flOld%Gij(i1,i2)
                GAX = GAX + this%nodes(c_node)%Tji(s1)*this%flOld%Gij(i2,i1)
            enddo
            this%flNew%Gij(k,no_fl) = -GXA*GAA
            this%flNew%Gij(no_fl,k) = -GAX*GAA
        enddo
        ! Update rest of the matrix: See Eq. (11)
        do k = 1 , no_fl-1
            s1 = this%flNew%Iij(k)
            i1 = this%nodes(s1)%fl_id
        do p = 1 , no_fl-1
            s2 = this%flNew%Iij(p)
            i2 = this%nodes(s2)%fl_id
            this%flNew%Gij(k,p) = this%flOld%Gij(i1,i2)  &
                                + this%flNew%Gij(k,no_fl)*this%flNew%Gij(no_fl,p)/GAA
        enddo
        enddo
        ! Fill last element of the matrix
        this%flNew%Gij(no_fl,no_fl) = GAA

        t_step2 = t_step2 + get_clock() - t_tmp
        t_tmp   = get_clock()

        ! Same calculation as above but for the selected elements

        do c_elem = 1 , no_elems
            call tmp_element%init_element(no_fl)
            tmp_element%Gij(no_fl,1) = e_Gxa(c_elem) ! this%elOld(c_elem)%Gxa
            tmp_element%Gij(no_fl,2) = e_Gax(c_elem) !this%elOld(c_elem)%Gax
!            ! here we use  equation: (11)
            GAX = tmp_element%Gij(no_fl,1)/GAA
            GXA = tmp_element%Gij(no_fl,2)/GAA
            do k = 1 , no_fl-1
                s1 = this%flNew%Iij(k)
                i1 = this%nodes(s1)%fl_id
                tmp_element%Gij(k,1) = this%elOld(c_elem)%Gij(i1,1) + GAX*this%flNew%Gij(no_fl,k)
                tmp_element%Gij(k,2) = this%elOld(c_elem)%Gij(i1,2) + this%flNew%Gij(k,no_fl)*GXA
            enddo
!            ! copy new to old
            call this%elOld(c_elem)%init_element(no_fl)
            this%elOld(c_elem)%Gij = tmp_element%Gij
        enddo


        t_step2e = t_step2e + get_clock() - t_tmp
        t_tmp   = get_clock()
        ! recalculate mapping from global to local in the new frontline
        k = 0
        do p = 1, this%flOld%no_sites
            if(this%nodes(this%flOld%Iij(p))%no_used < this%nodes(this%flOld%Iij(p))%no_bonds )then
                k = k + 1
                this%nodes(this%flOld%Iij(p))%fl_id = k
            endif
        enddo
        if(this%nodes(c_node)%no_used < this%nodes(c_node)%no_bonds) then
            k = k + 1
            this%nodes(c_node)%fl_id = k
        endif
        ! copy new to old
        call this%flOld%copy_frontline(this%flNew)
        t_step3 = t_step3 + get_clock() - t_tmp
    enddo ! end of loop c_node

    ! Copy calulated values to output array
    do c_elem = 1 , no_elems
        values(c_elem) = this%elOld(c_elem)%Gp
    enddo

    ! Freeing memory
    do i = 1 , this%no_elements
        call this%elNew(i)%free_element()
        call this%elOld(i)%free_element()
    enddo
    if(allocated(this%elNew))deallocate(this%elNew)
    if(allocated(this%elOld))deallocate(this%elOld)
    this%no_elements = 0
    t_total = get_clock() - t_total
    if(SKMINV_DEBUG_LEVEL == 2) then
        print*,"skminv::Profilling info calculation times for selected stages:"
        print*,"    Initialization time  :",t_step0,"[s]"
        print*,"    Step 1. for interface:",t_step1,"[s]"
        print*,"    Step 1. for elements :",t_step1e,"[s]"
        print*,"    Step 2. for interface:",t_step2,"[s]"
        print*,"    Step 2. for elements :",t_step2e,"[s]"
        print*,"    Renumbering nodes    :",t_step3,"[s]"
        print*,"    Total time           :",t_total,"[s]"
    endif
    if(SKMINV_DEBUG_LEVEL > 0) then
        print*,"skminv::finished"
    endif
!    if(allocated(Gij)) deallocate(Gij)
end subroutine invert

! ------------------------------------------------------------------------
!                           SKMINV
! Initialize memory, and convert matrix from row-col representation
! to arrays of nodes.
! no_nodes - number of nodes, i.e. size of the matrix.
! no_vals  - lenght of the array vals(:) and rc(:,1:2)
! vals(:)  - array which contains non-zero elements of the matrix
! rc(:,row/col)- contains row=1 and col=2 id to non zero values of matrix.
! ------------------------------------------------------------------------
subroutine init_memory(this,no_nodes,no_vals,vals,rc)
    class(skminv) :: this
    integer       :: no_nodes,no_vals
    complex*16    :: vals(:)
    integer       :: rc(:,:)
    ! local variables
    integer    :: i,r,c,k

    ! Free memory before initializating another one
    call this%free_memory()
    this%no_nodes = no_nodes
    allocate(this%nodes(no_nodes))

    ! clear memory of nodes
    do i = 1 , no_nodes
        call this%nodes(i)%reset_node()
    enddo
    ! caclulate number of hopings for each row/node
    do i = 1 , no_vals
        r = rc(i,1)
        c = rc(i,2)
        if(r /= c) this%nodes(r)%no_bonds = this%nodes(r)%no_bonds + 1
    enddo
    ! reallocate memory to calculated number of bonds in previous step
    ! alloc_mem uses this%nodes(i)%no_bonds to allocate memory
    do i = 1 , no_nodes
        call this%nodes(i)%alloc_mem()
    enddo

    ! Fill matrices using provided arrays
    do i = 1 , no_vals
        r = rc(i,1)
        c = rc(i,2)
        if( r == c ) then
            this%nodes(r)%H0 = vals(i)
        else
            this%nodes(r)%no_used = this%nodes(r)%no_used + 1
            k = this%nodes(r)%no_used
            this%nodes(r)%Tij(k) = vals(i)
            this%nodes(r)%Iij(k) = c
        endif
    enddo
    ! Fill the Tji vectors, this increases the number of memory
    ! usage but in the main loop one does not have to search for
    ! the Tij coupling.
    do i = 1 , no_nodes
        this%nodes(i)%no_used = 0
        do k = 1 , this%nodes(i)%no_bonds
            r = this%nodes(i)%Iij(k)
            do c = 1 , this%nodes(r)%no_bonds
                if( this%nodes(r)%Iij(c) == i ) exit
            enddo
            this%nodes(i)%Tji(k) = this%nodes(r)%Tij(c)
        enddo
    enddo
    ! Calculate number of neightboors with id smaller than current node id.
    do i = 1 , no_nodes
        this%nodes(i)%no_bonds_less = 0
        do k = 1 , this%nodes(i)%no_bonds
            if(this%nodes(i)%Iij(k) < i) this%nodes(i)%no_bonds_less = this%nodes(i)%no_bonds_less + 1
        enddo
        deallocate(this%nodes(i)%Iij_less)
        if(this%nodes(i)%no_bonds_less > 0) allocate(this%nodes(i)%Iij_less(this%nodes(i)%no_bonds_less))
    enddo
    ! Fill Iij_less array.
    do i = 1 , no_nodes
        r = 0
        if(this%nodes(i)%no_bonds_less /= 0) then
        do k = 1 , this%nodes(i)%no_bonds
            if(this%nodes(i)%Iij(k) < i) then
                r = r + 1
                this%nodes(i)%Iij_less(r) = this%nodes(i)%Iij(k)
            endif
        enddo
        endif
    enddo
    if(SKMINV_DEBUG_LEVEL == 1) then
        print*,"skminv::printing couplings values:"
        do i = 1 , no_nodes
            print"(A,i,A,1000i)","row=",i," Iij:", this%nodes(i)%Iij
            do k = 1 , this%nodes(i)%no_bonds
                print"(A,i,A,i,A,2f12.5)","from ",i," to:",this%nodes(i)%Iij(k)," t(i,j)=",(this%nodes(i)%Tij(k))
            enddo
            do k = 1 , this%nodes(i)%no_bonds
                print"(A,i,A,i,A,2f12.5)","from ",this%nodes(i)%Iij(k)," to:",i," t(j,i)=",(this%nodes(i)%Tji(k))
            enddo
        enddo
    endif
    call this%flNew%init_frontline(1)
    call this%flOld%init_frontline(1)


end subroutine init_memory

! -------------------------------------------------------------- !
!
! -------------------------------------------------------------- !
subroutine reset_node(this)
    class(node) :: this
    this%no_bonds      = 0
    this%no_used       = 0
    this%no_bonds_less = 0
    this%fl_id         = 0
    this%H0            = 0
    if(allocated(this%Iij))      deallocate(this%Iij)
    if(allocated(this%Tij))      deallocate(this%Tij)
    if(allocated(this%Tji))      deallocate(this%Tji)
    if(allocated(this%Iij_less)) deallocate(this%Iij_less)
endsubroutine reset_node

! -------------------------------------------------------------- !
!
! -------------------------------------------------------------- !
subroutine alloc_mem(this)
    class(node) :: this
    allocate(this%Iij(this%no_bonds))
    allocate(this%Iij_less(this%no_bonds_less))
    allocate(this%Tij(this%no_bonds))
    allocate(this%Tji(this%no_bonds))
    this%Iij      = 0
    if(this%no_bonds_less > 0) this%Iij_less = 0

    this%Tij = 0
    this%Tji = 0
endsubroutine alloc_mem

! -------------------------------------------------------------- !
!
! -------------------------------------------------------------- !
subroutine init_frontline(this,ns)
    class(frontline) :: this
    integer :: ns
    if(ns /= this%no_sites) then
        this%no_sites = ns
        if(allocated(this%Iij)) deallocate(this%Iij)
        if(allocated(this%Gij)) deallocate(this%Gij)
        allocate(this%Iij(this%no_sites))
        allocate(this%Gij(this%no_sites,this%no_sites))
    endif
    this%Iij = 0
    this%Gij = 0
endsubroutine init_frontline


subroutine free_frontline(this)
    class(frontline) :: this
    this%no_sites = 0
    if(allocated(this%Iij)) deallocate(this%Iij)
    if(allocated(this%Gij)) deallocate(this%Gij)
endsubroutine free_frontline
! -------------------------------------------------------------- !
!
! -------------------------------------------------------------- !
subroutine copy_frontline(this,scr_ft)
    class(frontline) :: this
    type(frontline) :: scr_ft
    call this%init_frontline(scr_ft%no_sites)
    this%Gij      = scr_ft%Gij
    this%Iij      = scr_ft%Iij
    this%no_sites = scr_ft%no_sites
endsubroutine copy_frontline

! -------------------------------------------------------------- !
!
! -------------------------------------------------------------- !
subroutine init_element(this,ns)
    class(element) :: this
    integer :: ns
    if(ns /= this%no_sites) then
        this%no_sites = ns
        if(allocated(this%Gij)) deallocate(this%Gij)
        allocate(this%Gij(this%no_sites,2))
    endif
!    this%Gij = 0
endsubroutine init_element

subroutine free_element(this)
    class(element) :: this
    this%no_sites = 0
    this%i        = 0
    this%j        = 0
    this%Gp       = 0
    this%Gax      = 0
    this%Gxa      = 0
    if(allocated(this%Gij)) deallocate(this%Gij)
endsubroutine free_element

subroutine copy_element(this,scr_ft)
    class(element) :: this
    type(element) :: scr_ft

    call this%init_element(scr_ft%no_sites)
    this%Gij      = scr_ft%Gij
    this%no_sites = scr_ft%no_sites
endsubroutine copy_element

end module modskminv

