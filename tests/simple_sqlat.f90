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

 use modunits ! unit conversion tools
 use modsys   ! eigen values
 use modlead  ! bandgap structure
 use modshape
 implicit none
 character(*),parameter :: output_folder = "simple_sqlat_output/"
 type(qsys)                 :: qsystem
 type(qshape)               :: lead_shape
 type(qlead)                :: lead
 doubleprecision            :: zeros_vector(200)
 doubleprecision            :: a_dx,a_Emin,a_Emax,a_Bz
 integer                    :: no_expected_states
 integer ,parameter         :: nx = 50
 integer ,parameter         :: ny = 25
 doubleprecision,parameter  :: dx = 5.0 ! [nm]
 integer , dimension(nx,ny) :: gindex ! converts local index (i,j) to global index
 integer :: i,j,k
 doubleprecision            :: lead_translation_vec(3)


 ! Use atomic units in effective band model -> see modunit.f90 for more details
 call modunits_set_GaAs_params()
 a_dx = dx * L2LA ! convert it to atomic units
 a_Bz = BtoAtomicB(1.0D0) ! 1 Tesla to atomic units

 ! Initalize system
 call qsystem%init()

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
    ! qsystem%qatom - is a auxiliary variable of type(qatom), you can use your own
    !                 if you want.

    call qsystem%qatom%init((/ (i-1) * dx , (j-1) * dx + (i-1) * dx * 0.0 , 0.0 * dx /))

    !here e.g. you can diable some part of atoms
!    if( sqrt(abs(i - nx/2.0-1)**2 + abs(j - ny/2.0-1)**2)*dx < 20 ) then
!        qsystem%qatom%bActive = .false. ! do not include those atoms in calculations
!    endif

    ! Add atom to the system.
    call qsystem%add_atom(qsystem%qatom)
    k           = k + 1
    gindex(i,j) = k
 enddo
 enddo

 ! ----------------------------------------------------------
 ! 2. Construct logical connections between sites on the mesh.
 ! ----------------------------------------------------------
 ! Set criterium for the nearest neightbours "radius" search algorithm.
 ! Same as above qsystem%qnnbparam is a auxiliary variable of type(nnb_params) - more details in modsys.f90
 ! This structure is responsible for different criteria of nearest neighbour searching
  qsystem%qnnbparam%box = (/2*dx,2*dx,0.0D0/) ! do not search for the sites far than (-dx:+dx) direction
 ! Setup connections between sites with provided by you function "connect", see below for example.
 call qsystem%make_lattice(qsystem%qnnbparam,c_default=connect)

 ! Save lattice to file to see if constructed system is OK!
 ! Use plot_lattice.py to see the results.
 call qsystem%save_lattice(output_folder//"lattice.xml")

 ! ----------------------------------------------------------
 ! 3. Find eigenvalues of the system.
 ! You must provide range of energy (Emin,Emax) and
 ! number of expected states in that region. Eigen problem
 ! is solved with FEAST library.
 ! ----------------------------------------------------------
 no_expected_states = 80
 a_Emin = 0.0  / E0 / 1000.0 ! converting from [meV] to atomic units
 a_Emax = 5.0  / E0 / 1000.0 ! converting from [meV] to atomic units
 call qsystem%calc_eigenproblem(a_Emin,a_Emax,no_expected_states)


 ! Write founded eigenstates to file if there are states in the range (Emin,Emax)
 if(qsystem%no_eigenvalues > 0) then
 zeros_vector = 0
 open(unit=222,file=output_folder//"eigenvecs.dat")
 do i = 1 , nx
 do j = 1 , ny
    if(qsystem%atoms(gindex(i,j))%bActive)then
    ! get unique ID of current site
    k = qsystem%atoms(gindex(i,j))%globalIDs(1) ! globalIDs(1) - get first spin state. In case when
                                                ! when you want to add multiple states per site
                                                ! you can use this array to separate sites with different
                                                ! spins.
    write(222,"(500e20.6)"),qsystem%atoms(gindex(i,j))%atom_pos(:),abs(qsystem%eigenvecs(k,:))**2
    else
    ! fill with zeros empty spaces
    write(222,"(500e20.6)"),qsystem%atoms(gindex(i,j))%atom_pos(:),zeros_vector(1:qsystem%no_eigenvalues)
    endif
 enddo
    write(222,*),"" ! GNUPLOT necessary line break
 enddo
 close(222)
 endif ! end of if are there eigenstates?


 ! ----------------------------------------------------------
 ! 4. Use generated mesh to calculate the band structure
 ! in the region of homogenous lead.
 ! ----------------------------------------------------------
 ! Setup shape and initialize lead
 call lead_shape%init_rect(SHAPE_RECTANGLE_XY,-0.5*dx,0.5*dx,-0.5*dx,(ny+1)*dx)
 ! Lead needs to know where it is (lead_shape) and using this information it will
 ! create propper matrices and connections using list of atoms
 lead_translation_vec = (/  dx , 0.0D0 , 0.0D0 /)
 call lead%init_lead(lead_shape,lead_translation_vec,qsystem%atoms)
 a_Emin = 0.0 / A0 / 1000.0 ! converting from [meV] to atomic units
 call lead%print_lead(output_folder//"lead.xml",qsystem%atoms)
 a_Emax = 0.3 / A0 / 1000.0 ! converting from [meV] to atomic units
 call lead%bands(output_folder//"bands.dat",-M_PI,+M_PI,M_PI/50.0,a_Emin,a_Emax)
 ! ----------------------------------------------------------
 ! X. Clean memory...
 ! ----------------------------------------------------------
 call lead%destroy()
 call qsystem%destroy()
 print*,"Generating plots..."
 print*,"Plotting band structure..."
 call system("cd "//output_folder//"; ./plot_bands.py")
 print*,"Plotting eigenvectors..."
 call system("cd "//output_folder//"; ./plot_eigenvecs.py")
 print*,"Use Viewer program to see the structure and crated lead."

 contains

! ---------------------------------------------------------------------------
! This function decides if site A (called here atomA) with spin s1 has hoping
! to atom B with spin s2, and what is the value of the coupling.
! If there is no interaction between them returns false, otherwise true.
! ---------------------------------------------------------------------------
logical function connect(atomA,atomB,s1,s2,coupling_val)
    use modatom
    implicit none
    type(qatom) :: atomA,atomB
    integer     :: s1,s2
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

    ! We assume that there is no spin taken into account so in our
    ! case is always s1 = 1 and s2 = 1, thus s1 and s2 are not used
    ! here.
    t0 = 1/(2*m_eff*a_dx**2)

    if( xdiff == 0 .and. ydiff == 0 ) then
        connect      = .true.
        coupling_val = 4*t0
    else if( abs(xdiff) ==  1 .and. ydiff == 0 ) then
        connect = .true.
        ! magnetic Field enters the hamiltonian by phase: EXP(i*DX*Bx*y)
        y = (atomA%atom_pos(2) - ny/2*dx ) * L2LA ! convert to atomic units
        coupling_val = -t0 * EXP(II*xdiff*a_dx*a_Bz*y)
    else if( xdiff ==  0 .and. abs(ydiff) == 1 ) then
        connect = .true.
        coupling_val = -t0
    endif


end function connect

end program transporter
