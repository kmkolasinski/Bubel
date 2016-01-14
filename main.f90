! ------------------------------------------------------ !
! Quantulaba - simple_sqlat.f90 - Krzysztof Kolasinski 2015
!
! We want to find few first eigenvalues of rectangular
! system defined within effective mass Shroedinger
! equation. Magnetic field will be included into account.
! We assume that after the finite discretization nodes
! of the mesh are separated by dx=5nm distance in
! each direction.
! ------------------------------------------------------ !

program transporter
use modunits
use modscatter
implicit none
character(*),parameter :: output_folder = "./"
type(qscatter)             :: qt
doubleprecision            :: a_dx,x,y,hx,hy
integer ,parameter         :: nx = 100
integer ,parameter         :: ny = 100
doubleprecision,parameter  :: dx = 5.0 ! [nm]
integer , dimension(nx,ny) :: gindex ! converts local index (i,j) to global index
doubleprecision ,dimension(nx*ny) :: rho
integer :: i,j,k


! Use atomic units in effective band model -> see modunit.f90 for more details
QSYS_DEBUG_LEVEL = 1
call modunits_set_GaAs_params()
a_dx = dx * L2LA ! convert it to atomic units
hx = dx * nx / 2
hy = dx * ny / 2

! Initalize system
call qt%init_system()

! ----------------------------------------------------------
! 1. Create mesh - loop over width and height of the lattice
! ----------------------------------------------------------
k      = 0
gindex = 0
do i = 1 , nx
do j = 1 , ny
    x = (i-1) * dx
    y = (j-1) * dx
    call qt%qatom%init((/ x , y , 0.0 * dx /))
    ! Add atom to the system.
    call qt%qsystem%add_atom(qt%qatom)
    k           = k + 1
    gindex(i,j) = k
    rho(k) = exp( -0.01*( (x-hx+50)**2 + (y-hy)**2 ) ) - exp( -0.01*( (x-hx-50)**2 + (y-hy)**2 ) )
enddo
enddo

! ----------------------------------------------------------
! 2. Construct logical connections between sites on the mesh.
! ----------------------------------------------------------
! Set criterium for the nearest neightbours "radius" search algorithm.
! Same as above qt%qnnbparam is a auxiliary variable of type(nnb_params) - more details in modsys.f90
! This structure is responsible for different criteria of nearest neighbour searching
qt%qnnbparam%box = (/2*dx,2*dx,0.0D0/) ! do not search for the sites far than (-dx:+dx) direction
! Setup connections between sites with provided by you function "connect", see below for example.
call qt%qsystem%make_lattice(qt%qnnbparam,c_simple=connect)

call qt%qsystem%calc_linsys(dvec=rho,pardiso_mtype=QSYS_LINSYS_PARDISO_REAL_NON_SYM)

open(unit=11,file="rho.txt")
do i = 1 , nx
do j = 1 , ny
    x = (i-1) * dx
    y = (j-1) * dx
    write(11,"(3f)"),x,y,rho(gindex(i,j))
enddo
    write(11,*),""
enddo
! ----------------------------------------------------------
! X. Clean memory...
! ----------------------------------------------------------
call qt%destroy_system()

contains

! ---------------------------------------------------------------------------
! This function decides if site A (called here atomA) with spin s1 has hoping
! to atom B with spin s2, and what is the value of the coupling.
! If there is no interaction between them returns false, otherwise true.
! ---------------------------------------------------------------------------
logical function connect(atomA,atomB,coupling_val)
    use modcommons
    implicit none
    type(qatom) :: atomA,atomB

    complex*16  :: coupling_val ! you must overwrite this variable
    ! local variables
    integer         :: xdiff,ydiff
    doubleprecision :: dydiff,dxdiff,t0,y

    ! Calculate distance between atoms in units of dx.
    dxdiff = (atomA%atom_pos(1)-atomB%atom_pos(1))/dx
    dydiff = (atomA%atom_pos(2)-atomB%atom_pos(2))/dx


    ! Convert it to integers
    xdiff = NINT(dxdiff)
    ydiff = NINT(dydiff)

    ! default return value
    connect = .false.

    ! hoping parameter
    t0 = 1/(2*m_eff*a_dx**2)

    if( xdiff == 0 .and. ydiff == 0 ) then
        connect      = .true.
        coupling_val = 4*t0
    else if( abs(xdiff) ==  1 .and. ydiff == 0 ) then
        connect = .true.
        coupling_val = -t0
    else if( xdiff ==  0 .and. abs(ydiff) == 1 ) then
        connect = .true.
        coupling_val = -t0
    endif


end function connect

end program transporter
