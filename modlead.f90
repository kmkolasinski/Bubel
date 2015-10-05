! ---------------------------------------------------------------------------------------
!  Created by Krzysztof Kolasinski - 2015
!
!
! ---------------------------------------------------------------------------------------
module modlead
use modshape
use modsys
use modatom
implicit none
private

type qlead
    type(qshape) :: lead_shape
    complex*16,dimension(:,:) , allocatable :: valsH0 ! Hamiltonian of the lead
!    integer,   dimension(:,:) , allocatable :: rowcolH0

    complex*16,dimension(:,:) , allocatable :: valsTau ! Hamiltonian of the lead
  !  integer,   dimension(:,:) , allocatable :: rowcolTau
    integer,dimension(:),allocatable        :: l2g,l2a ! local id to global id
    contains
    procedure,pass(this) :: init_lead!()
    procedure,pass(this) :: destroy!()
endtype qlead

public :: qlead
contains

subroutine init_lead(this,lshape,all_atoms)
    class(qlead) :: this
    type(qshape) :: lshape
    type(qatom),dimension(:) :: all_atoms
    integer ,allocatable ,dimension(:) :: tmp_g2l
    integer :: i,j,b
    integer :: no_sites,system_size

    this%lead_shape = lshape

    system_size = 0
    do i = 1 , size(all_atoms)
        if(all_atoms(i)%bActive) then
            system_size = system_size + all_atoms(i)%no_in_states
        endif ! end of active atom
    enddo

    allocate(tmp_g2l(system_size))

    no_sites = 0
    do i = 1 , size(all_atoms)
        if(all_atoms(i)%bActive) then
            if( lshape%is_inside(all_atoms(i)%atom_pos) == .true. ) then
                do j = 1 , all_atoms(i)%no_in_states
                    no_sites = no_sites + 1
                    tmp_g2l(all_atoms(i)%globalIDs(j)) = no_sites
                enddo
            endif
        endif
    enddo
    print*,"H0 size:",no_sites

    if(allocated(this%valsH0))deallocate(this%valsH0)
    !if(allocated(this%rowcolH0))deallocate(this%rowcolH0)

    if(allocated(this%valsTau))deallocate(this%valsTau)
    if(allocated(this%l2g))deallocate(this%l2g)
    if(allocated(this%l2a))deallocate(this%l2a)
    !if(allocated(this%rowcolTau))deallocate(this%rowcolTau)

    allocate(this%valsH0 (no_sites,no_sites))
    allocate(this%valsTau(no_sites,no_sites))
    allocate(this%l2g(no_sites))
    allocate(this%l2a(no_sites))
    this%valsH0   = 0
    this%valsTau  = 0
    this%l2g = 0

    no_sites = 0
    do i = 1 , size(all_atoms)
        if(all_atoms(i)%bActive) then
            if( lshape%is_inside(all_atoms(i)%atom_pos) == .true. ) then
                do j = 1 , all_atoms(i)%no_in_states
                    no_sites = no_sites + 1
                    this%l2g(no_sites) = all_atoms(i)%globalIDs(j)
                    this%l2a(no_sites) = i
                enddo
            endif
        endif
    enddo


    do i = 1 , no_sites
        do b = 1 , all_atoms(this%l2g(i))%
        if( lshape%is_inside(all_atoms(i)%atom_pos) == .true. ) then
            do j = 1 , all_atoms(i)%no_in_states

            enddo
        endif

    enddo



    deallocate(tmp_g2l)

end subroutine init_lead


subroutine destroy(this)
    class(qlead) :: this
    if(allocated(this%valsH0))deallocate(this%valsH0)
!    if(allocated(this%rowcolH0))deallocate(this%rowcolH0)
    if(allocated(this%valsTau))deallocate(this%valsTau)
    if(allocated(this%l2g))deallocate(this%l2g)
    if(allocated(this%l2a))deallocate(this%l2a)
!    if(allocated(this%rowcolTau))deallocate(this%rowcolTau)
end subroutine destroy


endmodule modlead
