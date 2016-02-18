! --------------------------------------------------------------- !
! Bubel - simple_galerkin1d.f90 - Krzysztof Kolasinski 2015
!
! We want to construct hamiltonian using expansion in series
! of gaussian distributions in the 1D wire. In order to obtain
! the TB like hamiltonian one has to expand solution in terms
! narrow gaussians. \Psi(x) = \sum c_k \gauss(k,x)
! Putting this into the Shroedinger equation then by projecting
! On \gauss(l,x) we obtain the algebraic equation for c_k coeficients.
! Note that in this case on the right side of the H\Psi = E\Psi
! equation one has to also calulate so-caled overlap intergrals
! Slk = \int gauss(l,x) gauss(k,x) which leads to generalized
! eigen value problem.
! --------------------------------------------------------------- !
program galerkin1d
use modunits
use modshape
use modscatter
implicit none
character(*),parameter     :: output_folder = "simple_galerkin1d_output/"
type(qshape)               :: lead_shape
type(qscatter)             :: qt
integer ,parameter         :: nx = 400 ! number of gaussian points in wire
doubleprecision,parameter  :: dx = 0.1 ! spatial distribution of gaussians [nm]
integer                    :: i,j,k,m,l,nnsize
doubleprecision            :: lead_translation_vec(3),x,y,z,dens,xa,xb
doubleprecision            :: a_Emin,a_Emax,Ef,xhalf

doubleprecision            :: mSkl(-10:10),mTkl(-10:10)
integer :: update = 0


! Initalize system
call qt%init_system()
! ----------------------------------------------------------
! 1. Create wire
! ----------------------------------------------------------
do i = 1 , nx
    call qt%qatom%init((/ (i-1)*dx , 0.0D0, 0.0D0 /))
    call qt%qsystem%add_atom(qt%qatom)
enddo
xhalf = nx*dx/2
! ----------------------------------------------------------
! 2. Construct logical connections between sites on the mesh.
! ----------------------------------------------------------
nnsize = 3 ! restric base to three nearest neigthbours
! Precalculate Kinetic energy term and Gaussian overlaps
do i = -10 , 0 ,1
    xa = 0
    xb = i*dx
    mSkl(i) = Skl(xa,xb,dx/4)
    mTkl(i) = Tkl(xa,xb,dx/4)
enddo
do i = 1 , 10 ,1
    mTkl(i) = mTkl(-i)
    mSkl(i) = mSkl(-i)
enddo
! Restric NNB looking area - here 1D
qt%qnnbparam%box = (/(nnsize+2)*dx,0.0D0,0.0D0/)
! Contruct logical connections between "atoms"
! o_simple=connect_overlaps - corresponds to similar function which calculates the
!                             overlaps between "atoms"
call qt%qsystem%make_lattice(qt%qnnbparam,c_simple=connect,o_simple=connect_overlaps)

! ----------------------------------------------------------
! 3. Add lead to the system - note the unit cell size is now
! larger because of the larger number of nearest neigbours
! taken into account: nnsize*dx
! ----------------------------------------------------------
lead_translation_vec = (/  nnsize*dx , 0.0D0 , 0.0D0 /)
call lead_shape%init_range_3d((/ -0.5*dx , 0.0D0 , 0.0D0 /),lead_translation_vec)
call qt%add_lead(lead_shape,lead_translation_vec)

! To check save band structure file for the first lead
a_Emin =   0.0
a_Emax = 500.0
call qt%leads(1)%bands(output_folder//"bands.dat",-M_PI,+M_PI,M_PI/160.0,a_Emin,a_Emax)

lead_translation_vec = (/ -nnsize*dx , 0.0D0 , 0.0D0 /)
call lead_shape%init_range_3d((/ (+0.5+nx-1)*dx , 0.0D0 , 0.0D0 /),lead_translation_vec)
call qt%add_lead(lead_shape,lead_translation_vec)

! ----------------------------------------------------------
! 4. Save all the informations to one file. Use viewer to
! debug created lattice.
! ----------------------------------------------------------
call qt%save_system(output_folder//"system.xml")
open(unit=111,file=output_folder//"T.dat")
do Ef = 28.5 , 32.0 , 0.02

    call qt%qsystem%update_lattice(c_simple=connect,o_simple=connect_overlaps)
    call qt%calculate_modes(Ef);
    call qt%solve(1,Ef); ! solve scattering problem for the first lead
    ! Save transmission probability to file
    write(111,*),Ef,sum(qt%Tn(:))
enddo
close(111)
! ----------------------------------------------------------
! X. Clean memory...
! ----------------------------------------------------------
call qt%destroy_system()

print*,"Generating plots..."
print*,"Plotting band structure..."
call system("cd "//output_folder//"; ./plot_bands.py")
print*,"Plotting Transmission..."
call system("cd "//output_folder//"; ./plot_T.py")
print*,"Use Quantulaba-Viewer to see the structure: system.xml."


contains

! ----------------------------------------------------------
! Implementation of the coupling between atoms
! ----------------------------------------------------------
logical function connect(atomA,atomB,cval)
    use modcommons
    implicit none
    type(qatom) :: atomA,atomB
    complex*16  :: cval ! you must overwrite this variable
    ! local variables
    integer         :: xdiff
    doubleprecision :: dxdiff,xa,xb

    ! Calculate distance between atoms in units of dx.
    dxdiff = (atomB%atom_pos(1)-atomA%atom_pos(1))/dx
    ! Convert it to integers
    xdiff = NINT(dxdiff)

    xa= atomA%atom_pos(1)
    xb= atomB%atom_pos(1)

    ! default return value
    connect = .false.
    cval    = 0.0
    ! in case of nearest neightbours set calulated hopings
    if(abs(xdiff)<=nnsize-1) then
        cval   = mTkl(xdiff)
        ! in the scattering region add potential
        if(xa > 2 .and. xa < nx*dx-2)  cval = cval + Vkl(min(xa,xb),max(xa,xb),dx/4)
        connect = .true.
    endif
end function connect


logical function connect_overlaps(atomA,atomB,cval)
    use modcommons
    implicit none
    type(qatom) :: atomA,atomB
    complex*16  :: cval ! you must overwrite this variable
    ! local variables
    integer         :: xdiff
    doubleprecision :: dxdiff

    ! Calculate distance between atoms in units of dx.
    dxdiff = (atomB%atom_pos(1)-atomA%atom_pos(1))/dx
    ! Convert it to integers
    xdiff = NINT(dxdiff)
    ! default return value
    connect_overlaps = .false.
    cval   = 0.0

    if(abs(xdiff) <= nnsize-1) then
        cval   = mSkl(xdiff)
        connect_overlaps = .true.
    endif
end function connect_overlaps

! ----------------------------------------------------------
! Gaussian basis function
! ----------------------------------------------------------
doubleprecision function gauss(x,x0,dx)
    doubleprecision  :: x,x0,dx
    gauss = sqrt(1.0/(dx*sqrt(2*M_PI))*exp(-0.5*((x-x0)/dx)**2))
end function gauss
! ----------------------------------------------------------
! Second derivative of gaussian function
! ----------------------------------------------------------
doubleprecision function d2gaussdx(x,x0,dx)
    doubleprecision  :: x,x0,dx
    doubleprecision,parameter :: ddx = 0.0001
    d2gaussdx = gauss(x+ddx,x0,dx) - 2*gauss(x,x0,dx) + gauss(x-ddx,x0,dx)
    d2gaussdx = d2gaussdx/ddx/ddx
end function d2gaussdx
! ----------------------------------------------------------
! Calculation of the overlap elements
! ----------------------------------------------------------
doubleprecision function Skl(xa,xb,dx)
    doubleprecision  :: xa,xb,dx
    doubleprecision xs,xf,deltax,x
    integer :: nsigma
    Skl     = 0
    nsigma  = 5 ! restrict area of integral estimation
    xs      = min(xa,xb)-nsigma*dx
    xf      = max(xa,xb)+nsigma*dx
    deltax  = dx/20
    do x = xs , xf , deltax
        Skl = Skl + gauss(x,xa,dx)*gauss(x,xb,dx)*deltax
    enddo
end function Skl
! ----------------------------------------------------------
! Kinetic energy term
! ----------------------------------------------------------
doubleprecision function Tkl(xa,xb,dx)
    doubleprecision  :: xa,xb,dx
    doubleprecision xs,xf,deltax,x
    integer :: nsigma
    Tkl     = 0
    nsigma  = 5
    xs      = min(xa,xb)-nsigma*dx
    xf      = max(xa,xb)+nsigma*dx
    deltax  = dx/20
    do x = xs , xf , deltax
        Tkl = Tkl - 0.5*gauss(x,xa,dx)*d2gaussdx(x,xb,dx)*deltax
    enddo
end function Tkl
! ----------------------------------------------------------
! Example potential in the system: Two gaussians
! ----------------------------------------------------------
doubleprecision function pot(x)
    doubleprecision :: x,xhalf
    xhalf = nx*dx/2
    pot =  2*(exp(-0.2*((x-xhalf+10)**2)) + exp(-0.2*((x-xhalf-10)**2)))
end function pot
! ----------------------------------------------------------
! Overlaps with the potential
! ----------------------------------------------------------
doubleprecision function Vkl(xa,xb,dx)
    doubleprecision  :: xa,xb,dx
    doubleprecision xs,xf,deltax,x
    integer :: nsigma
    Vkl     = 0
    nsigma  = 5
    xs      = min(xa,xb)-nsigma*dx
    xf      = max(xa,xb)+nsigma*dx
    deltax  = dx/20
    do x = xs , xf , deltax
        Vkl = Vkl + gauss(x,xa,dx)*gauss(x,xb,dx)*pot(x)*deltax
    enddo
end function Vkl
end program galerkin1d
