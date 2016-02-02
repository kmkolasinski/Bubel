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

 use modunits   ! unit conversion tools
 use modscatter ! eigen values and transport
 use modlead    ! bandgap structure
 use modshape
 use modutils
 implicit none
 character(*),parameter :: output_folder = "plots/"
 type(qscatter)             :: qt
 type(qshape)               :: lead_shape
 doubleprecision            :: zeros_vector(200)
 doubleprecision            :: a_dx,a_Emin,a_Emax,a_Bz
 integer                    :: no_expected_states
 integer ,parameter         :: nx = 50
 integer ,parameter         :: ny = 50
 doubleprecision,parameter  :: dx = 5.0 ! [nm]
 integer , dimension(nx,ny) :: gindex ! converts local index (i,j) to global index
 integer :: i,j,k,p
 doubleprecision            :: lead_translation_vec(3),Ef



 ! Use atomic units in effective band model -> see modunit.f90 for more details
 call modunits_set_GaAs_params()
 a_dx = dx * L2LA ! convert it to atomic units
 a_Bz = BtoAtomicB(0.0D0) ! 1 Tesla to atomic units

 ! Initalize system
 call qt%init_system()

 ! ----------------------------------------------------------
 ! 1. Create mesh - loop over width and height of the lattice
 ! ----------------------------------------------------------
 k      = 0
 gindex = 0
 do i = 1 , nx
 do j = 1 , ny
    ! Initalize atom structure with position of the atom.
    ! For more details see modsys.f90 and qatom structure parameters.
    ! We assume that the 2D lattice lies at z=0.
    ! qt%qatom - is a auxiliary variable of type(qatom), you can use your own
    !                 if you want.

    call qt%qatom%init((/ (i-1) * dx , (j-1) * dx + (i-1) * dx * 0.0 , 0.0 * dx /))

    !here e.g. you can diable some part of atoms
!    if( sqrt(abs(i - nx/2.0-1)**2 + abs(j - ny/2.0-1)**2)*dx < 20 ) then
!        qt%qatom%bActive = .false. ! do not include those atoms in calculations
!    endif

    ! Add atom to the system.
    call qt%qsystem%add_atom(qt%qatom)
    k           = k + 1
    gindex(i,j) = k
 enddo
 enddo

! ----------------------------------------------------------
! 2. Construct logical connections between sites on the mesh.
! ----------------------------------------------------------
! Set criterium for the nearest neightbours "radius" search algorithm.
! Same as above qt%qnnbparam is a auxiliary variable of type(nnb_params) - more details in modsys.f90
! This structure is responsible for different criteria of nearest neighbour searching
!qt%qnnbparam%box = (/1.2*dx,1.2*dx,0.0D0/) ! do not search for the sites far than (-dx:+dx) direction
! Setup connections between sites with provided by you function "connect", see below for example.
qt%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_DISTANCE
qt%qnnbparam%distance = 4*dx
QSYS_DEBUG_LEVEL = 1
p = 2
call qtools_fd_template(v_dnFdxn,2,p)
call qtools_fd_template(v_dnFdyn,2,p)


call qt%qsystem%make_lattice(qt%qnnbparam,c_simple=connect)



! ----------------------------------------------------------
! 4. Use generated mesh to calculate the band structure
! in the region of homogenous lead.
! ----------------------------------------------------------
! Setup shape and initialize lead
call lead_shape%init_rect(SHAPE_RECTANGLE_XY,-0.5*dx,0.5*dx,-0.5*dx,(ny+1)*dx)
! Lead needs to know where it is (lead_shape) and using this information it will
! create propper matrices and connections using list of atoms
lead_translation_vec = (/  p*dx , 0.0D0 , 0.0D0 /)
call lead_shape%init_range_3d((/-0.5*dx,0.0D0,0.0D0/),lead_translation_vec)

call qt%add_lead(lead_shape,lead_translation_vec)

! Add second lead at the end of the system
call lead_shape%init_rect(SHAPE_RECTANGLE_XY,(-0.5+nx-1)*dx,(0.5+nx-1)*dx,-0.5*dx,(ny+1)*dx)

lead_translation_vec = (/  -p*dx , 0.0D0 , 0.0D0 /)
call lead_shape%init_range_3d((/(nx-0.5)*dx,0.0D0,0.0D0/),lead_translation_vec)
call qt%add_lead(lead_shape,lead_translation_vec)

a_Emin =  0.0 / E0 / 1000.0 ! converting from [meV] to atomic units
a_Emax = 300.0 / E0 / 1000.0 ! converting from [meV] to atomic units
call qt%leads(1)%bands(output_folder//"bands.dat",-M_PI,+M_PI,M_PI/50.0,a_Emin,a_Emax)

call qt%save_system(output_folder//"system.xml")

! Solve scattering problem for Ef=0.001
Ef = 5/E0/1000.0
QSYS_DEBUG_LEVEL = 1 ! show more info
call qt%calculate_modes(Ef)
call qt%solve(1,Ef)
! Save calculated electron density to file
do i = 1 , size(qt%qauxvec)
    qt%qauxvec(i) = sum(qt%densities(:,i))
enddo
call qt%qsystem%save_data(output_folder//"densities.xml",array2d=qt%densities,array1d=qt%qauxvec)
! Perform scan in function of Energy



!open(unit=111,file=output_folder//"T.dat")
!print*,"Performing energy scan..."
!QSYS_DEBUG_LEVEL = 1 ! show more info
!do Ef = 0.0 , 0.001 , 0.000025
!    ! Update hamiltonian elemenents value
!    call qt%qsystem%update_lattice(c_simple=connect)
!    call qt%calculate_modes(Ef)
!    call qt%solve(1,Ef)
!
!    print*,Ef
!    write(111,"(100f20.6)"),Ef,sum(qt%Tn(:))
!enddo
!close(111)


! ----------------------------------------------------------
! X. Clean memory...
! ----------------------------------------------------------
call lead_shape%destroy_shape()
call qt%destroy_system()
!print*,"Generating plots..."
!print*,"Plotting band structure..."
call system("cd "//output_folder//"; python plot_bands.py")
!print*,"Plotting eigenvectors..."
!call system("cd "//output_folder//"; ./plot_eigenvecs.py")
!print*,"Plotting Transmission..."
!call system("cd "//output_folder//"; ./plot_T.py")
!print*,"Use Viewer program to see the structure and created leads."
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
    integer         :: idx,idy
    doubleprecision :: dydiff,dxdiff,t0,y

    ! Calculate distance between atoms in units of dx.
    dxdiff = (atomA%atom_pos(1)-atomB%atom_pos(1))/dx
    dydiff = (atomA%atom_pos(2)-atomB%atom_pos(2))/dx
    ! Convert it to integers
    idx = NINT(dxdiff)
    idy = NINT(dydiff)
    ! default return value
    connect      = .true.
    coupling_val = 0
    ! hoping parameter
    t0 = 1/(2*m_eff*a_dx**2)


    coupling_val = -t0*(dnFdxn(2,idx,idy,0) + dnFdyn(2,idx,idy,0))

    if(abs(coupling_val) < qsys_double_error) connect = .false.

end function connect
end program transporter
