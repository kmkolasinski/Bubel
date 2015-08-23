! ---------------------------------------------------------------------------------------
!
!
!
! ---------------------------------------------------------------------------------------
module modatom

type atom
    doubleprecision :: atom_pos(3) ! position (x,y,z)
    integer         :: atom_id ! unique ID of atom
    integer         :: no_in_states ! number of internal states (e.g. spin degree of freedom)
end type atom

contains



end module modatom

! ---------------------------------------------------------------------------------------
!
!
!
! ---------------------------------------------------------------------------------------
module modunitcell
use modatom
implicit none
    integer,parameter :: NO_ATOMS_INC_VALUE = 50
type unitcell
    doubleprecision :: unit_vectors(3,3) ! unit base vectors (n,(x,y,z)) by default 0
    integer         :: no_atoms ! number of atoms in the unit cell
    type(atom),allocatable,dimension(:) :: atoms ! list of all atom in the unit cell
    contains
    procedure, public, pass(cell) :: init
end type unitcell

contains

subroutine init_cell(cell)
    class(unitcell)     :: cell
    cell%no_atoms     = 0
    cell%unit_vectors = 0
    allocate(cell%atoms(NO_ATOMS_INC_VALUE))
end subroutine init_cell

subroutine free_cell(cell)
    class(unitcell)     :: cell
    cell%no_atoms     = 0
    cell%unit_vectors = 0
    deallocate(cell%atoms)
end subroutine free_cell

end module modunitcell

! ---------------------------------------------------------------------------------------
!
!
!
! ---------------------------------------------------------------------------------------
module modsys
use modatom
implicit none


contains

end module modsys
