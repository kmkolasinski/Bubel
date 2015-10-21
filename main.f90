! ------------------------------------------------------ !
! Quantulaba - strange_graphene.f90 - Krzysztof Kolasinski 2015
!
! Graphene unit cell consists of two atoms A and B, those
! atoms can be viewed as a two spin states of some quantum
! object. In that case Graphene lattice can be presented
! as a simple square lattice with propper coupling between
! nearest neightbours.
! ------------------------------------------------------ !

program graphen

use modunits    ! unit conversion tools
use modscatter  ! eigen values and transport
use modlead     ! bandgap structure
use modshape
implicit none
character(*),parameter :: output_folder = "./"
type(qscatter)             :: qt
type(qshape)               :: lead_shape

doubleprecision            :: zeros_vector(200)
doubleprecision            :: a_dx,a_Emin,a_Emax
integer                    :: no_expected_states
integer ,parameter         :: nx = 50
integer ,parameter         :: ny = 40
doubleprecision,parameter  :: dx = 5.0 ! [nm]
doubleprecision,parameter  :: alpha30 =  30.0/180.0*M_PI

integer , dimension(nx,ny) :: gindex ! converts local index (i,j) to global index
integer :: i,j,k,m
doubleprecision            :: lead_translation_vec(3),x,y,z,Ef
character(300) :: line


! Initalize system
call qt%init_system()
! ----------------------------------------------------------
! 1. Create mesh - loop over width and height of the lattice
! ----------------------------------------------------------
k      = 0
gindex = 0
do i = 1 , nx
do j = 1 , ny
    ! Initalize atom structure with position of the atom. Now instead of atoms A and B we have
    ! one atom with two spins.
    call qt%qatom%init((/ (i-1) * dx , (j-1) * dx + (i-1) * dx * 0.0 , 0.0 * dx /),no_in_states=2)
    ! Add atom to the system.
    if(j == 5) qt%qatom%flag = 1
    call qt%qsystem%add_atom(qt%qatom)
    k           = k + 1
    gindex(i,j) = k
enddo
enddo
! ----------------------------------------------------------
! 2. Construct logical connections between sites on the mesh.
! ----------------------------------------------------------
qt%qnnbparam%box = (/2*dx,2*dx,0.0D0/)
call qt%qsystem%make_lattice(qt%qnnbparam,c_matrix=coupling)


! ----------------------------------------------------------
! 3. Find eigenvalues of the system.
! You must provide range of energy (Emin,Emax) and
! number of expected states in that region. Eigen problem
! is solved with FEAST library.
! ----------------------------------------------------------
!no_expected_states = 40
!! we choose small window
!a_Emin =-3.0+0.1
!a_Emax =-2.9+0.1
!call qt%qsystem%calc_eigenproblem(a_Emin,a_Emax,no_expected_states)
!
!! Write founded eigenstates to file if there are states in the range (Emin,Emax)
!if(qt%qsystem%no_eigenvalues > 0) then
!zeros_vector = 0
!open(unit=222,file=output_folder//"eigenvecs.dat")
!do i = 1 , nx
!do j = 1 , ny
!    ! Remapping to graphen lattice
!    x = (qt%qsystem%atoms(gindex(i,j))%atom_pos(1) + qt%qsystem%atoms(gindex(i,j))%atom_pos(2)*sin(alpha30))/dx
!    y = (qt%qsystem%atoms(gindex(i,j))%atom_pos(2) - sin(alpha30)*qt%qsystem%atoms(gindex(i,j))%atom_pos(1))/dx
!    z = 0
!    if(qt%qsystem%atoms(gindex(i,j))%bActive)then
!        ! get unique ID of current site
!        k = qt%qsystem%atoms(gindex(i,j))%globalIDs(2)  ! spin = 2
!        m = qt%qsystem%atoms(gindex(i,j))%globalIDs(1)  ! spin = 1
!        write(222,"(500e20.6)"),x,y,z,abs(qt%qsystem%eigenvecs(k,:))**2+abs(qt%qsystem%eigenvecs(m,:))**2
!    else
!        ! fill with zeros empty spaces
!        write(222,"(500e20.6)"),x,y,z,zeros_vector(1:qt%qsystem%no_eigenvalues)
!    endif
!enddo
!write(222,*),"" ! GNUPLOT necessary line break
!enddo
!close(222)
!endif ! end of if are there eigenstates?



! ----------------------------------------------------------
! 4. Use generated mesh to calculate the band structure
! in the region of homogenous lead.
! ----------------------------------------------------------
! Setup shape and initialize lead
lead_translation_vec = (/  dx , 0.0D0 , 0.0D0 /)
call lead_shape%init_range_3d((/ -0.5*dx , 0.0D0 , 0.0D0 /),lead_translation_vec)
call qt%add_lead(lead_shape,lead_translation_vec)

lead_translation_vec = (/  -dx , 0.0D0 , 0.0D0 /)
call lead_shape%init_range_3d((/ (nx+0.5-1)*dx , 0.0D0 , 0.0D0 /),lead_translation_vec)
call qt%add_lead(lead_shape,lead_translation_vec)



a_Emin =-3.0
a_Emax = 3.0
call qt%leads(1)%bands(output_folder//"bands.dat",-M_PI,+M_PI,M_PI/160.0,a_Emin,a_Emax)


! Save lattice to file to see if constructed system is OK!
call qt%save_system(output_folder//"system.xml")

Ef = 0.1
QSYS_DEBUG_LEVEL = 1 ! tell me more
call qt%calculate_modes(Ef)
call qt%solve(1,Ef)

! Save calculated electron density to file
do i = 1 , size(qt%qauxvec)
    qt%qauxvec(i) = sum(qt%densities(:,i))
enddo
call qt%qsystem%save_data(output_folder//"densities.xml",array2d=qt%densities,array1d=qt%qauxvec)


print*,"Performing energy scan..."
open(unit=111,file=output_folder//"T.dat")
QSYS_DEBUG_LEVEL = 0 ! show me less
do Ef = -3.0 , 3.025 , 0.051
    ! Update hamiltonian elemenents value
    call qt%qsystem%update_lattice(c_matrix=coupling)
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)

    print*,"Energy:",Ef,sum(qt%Tn(:))
    write(111,"(100f20.6)"),Ef,sum(qt%Tn(:))
enddo
close(111)


! ----------------------------------------------------------
! X. Clean memory...
! ----------------------------------------------------------

call qt%destroy_system()

print*,"Generating plots..."
print*,"Plotting band structure..."
call system("cd "//output_folder//"; ./plot_bands.py")
print*,"Plotting eigenvectors..."
call system("cd "//output_folder//"; ./plot_eigenvecs.py")
print*,"Plotting Transmission..."
call system("cd "//output_folder//"; ./plot_T.py")
print*,"Use Viewer program to see the structure and created leads."

contains

! ----------------------------------------------------------
! Define coupling between atoms
! ----------------------------------------------------------
logical function coupling(atomA,atomB,coupling_mat)
    use modcommons
    implicit none
    type(qatom) :: atomA,atomB
    complex*16  :: coupling_mat(:,:) ! you must overwrite this variable
    ! local variables
    integer         :: xdiff,ydiff
    doubleprecision :: dydiff,dxdiff,t0
    ! Calculate distance between atoms in units of dx.
    dxdiff = (atomA%atom_pos(1)-atomB%atom_pos(1))/dx
    dydiff = (atomA%atom_pos(2)-atomB%atom_pos(2))/dx
    ! Convert it to integers
    xdiff = NINT(dxdiff)
    ydiff = NINT(dydiff)

    ! default return value
    coupling       = .false.
    coupling_mat   = 0.0
    t0 = -1
    if( xdiff == 0 .and. ydiff == 0 ) then
        coupling      = .true.
        coupling_mat = t0 * MAT_SX
    else if( abs(xdiff) ==  1 .and. ydiff == 0 ) then
        coupling = .true.
        coupling_mat = 0.5*t0 *(MAT_SX + xdiff*II*MAT_SY)
    else if( xdiff ==  0 .and. abs(ydiff) == 1 ) then
        coupling = .true.
        coupling_mat = 0.5*t0 *(MAT_SX + ydiff*II*MAT_SY)
    endif
end function coupling

end program graphen


!! ------------------------------------------------------ !
!! Quantulaba - simple_sqlat.f90 - Krzysztof Kolasinski 2015
!!
!! We want to find few first eigenvalues of rectangular
!! system defined within effective mass Shroedinger
!! equation. Magnetic field will be included into account.
!! We assume that after the finite discretization nodes
!! of the mesh are separated by dx=5nm distance in
!! each direction.
!! ------------------------------------------------------ !
!
!program transporter
!
! use modunits   ! unit conversion tools
! use modscatter ! eigen values and transport
! use modlead    ! bandgap structure
! use modshape
! implicit none
! character(*),parameter :: output_folder = "./"
! type(qscatter)             :: qt
! type(qshape)               :: lead_shape
! doubleprecision            :: zeros_vector(200)
! doubleprecision            :: a_dx,a_Emin,a_Emax,a_Bz,t0
! integer                    :: no_expected_states
! integer ,parameter         :: nx = 60
! integer ,parameter         :: ny = 50
! doubleprecision,parameter  :: dx = 5.0 ! [nm]
! integer , dimension(nx,ny) :: gindex ! converts local index (i,j) to global index
! integer :: i,j,k,idup,iddown
! doubleprecision            :: lead_translation_vec(3),lead_base_vec(3),Ef,y,z,w,r
! complex*16                 :: cval
!
!
! ! Use atomic units in effective band model -> see modunit.f90 for more details
! call modunits_set_GaAs_params()
! a_dx = dx * L2LA ! convert it to atomic units
! a_Bz = BtoAtomicB(1.0D0) ! 1 Tesla to atomic units
! t0 = 1/(2*m_eff*a_dx**2)
! ! Initalize system
! call qt%init_system()
! QSYS_DEBUG_LEVEL = 1 ! show more info
! ! ----------------------------------------------------------
! ! 1. Create mesh - loop over width and height of the lattice
! ! ----------------------------------------------------------
! k      = 0
! gindex = 0
! do i = 1 , nx
! do j = 1 , ny
!    ! Initalize atom structure with position of the atom.
!    ! For more details see modsys.f90 and qatom structure parameters.
!    ! We assume that the 2D lattice lies at z=0.
!    ! qt%qatom - is a auxiliary variable of type(qatom), you can use your own
!    !                 if you want.
!    w = 2*M_PI*(j-1.0)/(ny)
!    r = dx*ny/2/M_PI
!    y = r*cos(w + 1)
!    z = r*sin(w + 1)
!
!!    z = 0
!!    y = (j-1.0)*dx
!    call qt%qatom%init((/ (i-1) * dx , y , z /))
!    if(j==20) qt%qatom%flag = 1
!!    if(j==ny/4) qt%qatom%bActive = .false.
!    !here e.g. you can diable some part of atoms
!    if( abs(j - ny/2.0) > 10 .and. i < 20 ) then
!!        qt%qatom%bActive = .false. ! do not include those atoms in calculations
!    endif
!
!    ! Add atom to the system.
!    call qt%qsystem%add_atom(qt%qatom)
!    k           = k + 1
!    gindex(i,j) = k
! enddo
! enddo
!
! ! ----------------------------------------------------------
! ! 2. Construct logical connections between sites on the mesh.
! ! ----------------------------------------------------------
!
! qt%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_DISTANCE
! qt%qnnbparam%distance   = dx*1.4
!call qt%qsystem%make_lattice(qt%qnnbparam,c_matrix=connect)
!call qt%save_system(output_folder//"system.xml")
!!cval = -t0
!!do i = 1 , nx
!!    idup   = gindex(i,ny)
!!    iddown = gindex(i,1)
!!    if(qt%qsystem%atoms(iddown)%bActive) then
!!        call qt%qsystem%atoms(iddown)%add_bond(1,idup,1,cval)
!!        call qt%qsystem%atoms(idup)%add_bond(1,iddown,1,cval)
!!    endif
!!enddo
!
! ! ----------------------------------------------------------
! ! 3. Find eigenvalues of the system.
! ! You must provide range of energy (Emin,Emax) and
! ! number of expected states in that region. Eigen problem
! ! is solved with FEAST library.
! ! ----------------------------------------------------------
!! print*,"Finding eigenvalues..."
!! no_expected_states = 80
!! a_Emin = 0.0  / E0 / 1000.0 ! converting from [meV] to atomic units
!! a_Emax = 1.0  / E0 / 1000.0 ! converting from [meV] to atomic units
!! call qt%qsystem%calc_eigenproblem(a_Emin,a_Emax,no_expected_states,print_info=1)
!!
!!
!! ! Write founded eigenstates to file if there are states in the range (Emin,Emax)
!! if(qt%qsystem%no_eigenvalues > 0) then
!! zeros_vector = 0
!! open(unit=222,file=output_folder//"eigenvecs.dat")
!! do i = 1 , nx
!! do j = 1 , ny
!!    if(qt%qsystem%atoms(gindex(i,j))%bActive)then
!!    ! get unique ID of current site
!!    k = qt%qsystem%atoms(gindex(i,j))%globalIDs(1) ! globalIDs(1) - get first spin/orbital state. In case when
!!                                                    ! when you want to add multiple states per site
!!                                                    ! you can use this array to separate sites with different
!!                                                    ! spins.
!!    write(222,"(500e20.6)"),qt%qsystem%atoms(gindex(i,j))%atom_pos(:),abs(qt%qsystem%eigenvecs(k,:))**2
!!    else
!!    ! fill with zeros empty spaces
!!    write(222,"(500e20.6)"),qt%qsystem%atoms(gindex(i,j))%atom_pos(:),zeros_vector(1:qt%qsystem%no_eigenvalues)
!!    endif
!! enddo
!!    write(222,*),"" ! GNUPLOT necessary line break
!! enddo
!! close(222)
!! endif ! end of if are there eigenstates?
!
!
!
!! ----------------------------------------------------------
!! 4. Use generated mesh to calculate the band structure
!! in the region of homogenous lead.
!! ----------------------------------------------------------
!! Setup shape and initialize lead
!
!! Lead needs to know where it is (lead_shape) and using this information it will
!! create propper matrices and connections using list of atoms
!lead_translation_vec = (/  dx , 0.0D0 , 0.0D0 /)
!lead_base_vec = (/-0.5*dx,0.0D0,0.0D0/)
!call lead_shape%init_range_3d(lead_base_vec,lead_translation_vec)
!call qt%add_lead(lead_shape,lead_translation_vec)
!
!! Add second lead at the end of the system
!
!
!lead_translation_vec = (/  -dx , 0.0D0 , 0.0D0 /)
!lead_base_vec = (/(nx+0.5-1)*dx,0.0D0,0.0D0/)
!call lead_shape%init_range_3d(lead_base_vec,lead_translation_vec)
!call qt%add_lead(lead_shape,lead_translation_vec)
!
!a_Emin =-0.1 / A0 / 1000.0 ! converting from [meV] to atomic units
!a_Emax = 0.3 / A0 / 1000.0 ! converting from [meV] to atomic units
!call qt%leads(1)%bands(output_folder//"bands.dat",-M_PI,+M_PI,M_PI/50.0,a_Emin,a_Emax)
!
!call qt%save_system(output_folder//"system.xml")
!
!!! Solve scattering problem for Ef=0.001
!Ef = 1.*t0
!call qt%calculate_modes(Ef)
!call qt%solve(1,Ef)
!! Save calculated electron density to file
!do i = 1 , size(qt%qauxvec)
!    qt%qauxvec(i) = sum(qt%densities(:,i))
!enddo
!call qt%qsystem%save_data(output_folder//"densities.xml",array2d=qt%densities,array1d=qt%qauxvec)
!!! Perform scan in function of Energy
!!open(unit=111,file=output_folder//"T.dat")
!!print*,"Performing energy scan..."
!!QSYS_DEBUG_LEVEL = 1 ! show more info
!!do Ef = 0.0 , 0.001 , 0.000025
!!    ! Update hamiltonian elemenents value
!!    call qt%qsystem%update_lattice(c_default=connect)
!!    call qt%calculate_modes(Ef)
!!    call qt%solve(1,Ef)
!!
!!    print*,Ef
!!    write(111,"(100f20.6)"),Ef,sum(qt%Tn(:))
!!enddo
!!close(111)
!
!
!! ----------------------------------------------------------
!! X. Clean memory...
!! ----------------------------------------------------------
!call lead_shape%destroy_shape()
!call qt%destroy_system()
!print*,"Generating plots..."
!print*,"Plotting band structure..."
!call system("cd "//output_folder//"; ./plot_bands.py")
!!print*,"Plotting eigenvectors..."
!!call system("cd "//output_folder//"; ./plot_eigenvecs.py")
!!print*,"Plotting Transmission..."
!!call system("cd "//output_folder//"; ./plot_T.py")
!print*,"Use Viewer program to see the structure and crated leads."
!contains
!
!! ---------------------------------------------------------------------------
!! This function decides if site A (called here atomA) with spin s1 has hoping
!! to atom B with spin s2, and what is the value of the coupling.
!! If there is no interaction between them returns false, otherwise true.
!! ---------------------------------------------------------------------------
!logical function connect(atomA,atomB,coupling_val)
!    use modcommons
!    implicit none
!    type(qatom) :: atomA,atomB
!    complex*16,dimension(:,:)  :: coupling_val ! you must overwrite this variable
!    ! local variables
!    integer         :: xdiff,ydiff
!    doubleprecision :: dydiff,dxdiff,dzdiff,y,dr
!
!    ! Calculate distance between atoms in units of dx.
!    dxdiff = (atomA%atom_pos(1)-atomB%atom_pos(1))/dx
!    dydiff = (atomA%atom_pos(2)-atomB%atom_pos(2))/dx
!    dzdiff = (atomA%atom_pos(3)-atomB%atom_pos(3))/dx
!
!    dr = sqrt(dxdiff**2 + dydiff**2 + dzdiff**2)
!    ! Convert it to integers
!    xdiff = NINT(dxdiff)
!    ydiff = NINT(dydiff)
!
!    ! default return value
!    connect = .false.
!
!
!
!    if( dr < 0.001*dx ) then
!        connect      = .true.
!        coupling_val = 4*t0 + atomA%flag*0.000001*t0
!
!!    else if( abs(xdiff) ==  1 .and.  abs(ydiff) == 0) then
!    else if( abs(xdiff) ==  1  ) then
!        connect = .true.
!        ! magnetic Field enters the hamiltonian by phase: EXP(i*DX*Bx*y)
!        y = (atomA%atom_pos(2) - ny/2*dx ) * L2LA ! convert to atomic units
!        coupling_val = -t0 * EXP(II*xdiff*a_dx*a_Bz*y)
!    else
!!    else if(abs(ydiff) == 1 .and. abs(xdiff) ==  0) then
!        connect = .true.
!        coupling_val = -t0
!    endif
!
!
!end function connect
!
!end program transporter
