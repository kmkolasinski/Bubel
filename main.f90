! ------------------------------------------------------ !
! Bubel - simple_sqlat.f90 - Krzysztof Kolasinski 2015
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
 use modskminv
 implicit none
 character(*),parameter :: output_folder = "./"
 type(qscatter)             :: qt
 type(qshape)               :: lead_shape
 type(skminv)               :: selinv
 doubleprecision            :: zeros_vector(200)
 doubleprecision            :: a_dx,a_Emin,a_Emax,a_Bz
 integer                    :: no_expected_states
 integer ,parameter         :: nx = 50
 integer ,parameter         :: ny = 20
 doubleprecision,parameter  :: dx = 5.0 ! [nm]
 integer , dimension(nx,ny) :: gindex ! converts local index (i,j) to global index
 integer :: i,j,k,melem,p
 doubleprecision            :: lead_translation_vec(3),Ef
 integer,allocatable,dimension(:,:) :: elements
 complex*16,allocatable,dimension(:) :: evals




 ! Use atomic units in effective band model -> see modunit.f90 for more details
 call modunits_set_GaAs_params()
 a_dx = dx * L2LA ! convert it to atomic units
 a_Bz = BtoAtomicB(1.5D0) ! 1 Tesla to atomic units

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
    if( sqrt(abs(i - nx/2.0-1)**2 + abs(j - ny/2.0-1)**2)*dx < 5 ) then
        qt%qatom%bActive = .false. ! do not include those atoms in calculations
    endif
    if(qt%qatom%bActive) then
        ! Add atom to the system.
        call qt%qsystem%add_atom(qt%qatom)
        k           = k + 1
        gindex(i,j) = k
    endif
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


! ----------------------------------------------------------
! 4. Use generated mesh to calculate the band structure
! in the region of homogenous lead.
! ----------------------------------------------------------
! Setup shape and initialize lead
call lead_shape%init_rect(SHAPE_RECTANGLE_XY,-0.5*dx,0.5*dx,-0.5*dx,(ny+1)*dx)
! Lead needs to know where it is (lead_shape) and using this information it will
! create propper matrices and connections using list of atoms
lead_translation_vec = (/  dx , 0.0D0 , 0.0D0 /)
call qt%add_lead(lead_shape,lead_translation_vec)

! Add second lead at the end of the system
call lead_shape%init_rect(SHAPE_RECTANGLE_XY,(-0.5+nx-1)*dx,(0.5+nx-1)*dx,-0.5*dx,(ny+1)*dx)
lead_translation_vec = (/  -dx , 0.0D0 , 0.0D0 /)
call qt%add_lead(lead_shape,lead_translation_vec)

a_Emin =  0.0 / E0 / 1000.0 ! converting from [meV] to atomic units
a_Emax = 150.0 / E0 / 1000.0 ! converting from [meV] to atomic units
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
call qt%solve(2,Ef)
do i = 1 , size(qt%qauxvec)
    qt%qauxvec(i) = (qt%qauxvec(i) +  sum(qt%densities(:,i)))/2
enddo
call qt%qsystem%save_data(output_folder//"densities.xml",array2d=qt%densities,array1d=qt%qauxvec)
! Perform scan in function of Energy
!SKMINV_DEBUG_LEVEL = 1
call selinv%init_memory(qt%NO_VARIABLES,qt%NO_NON_ZERO_VALUES,qt%MATHVALS,qt%ROWCOLID)
!QSYS_FORCE_SCHUR_DECOMPOSITION = .true.



do i = 1 , nx
do j = 1 , ny
    p = gindex(i,j)
    if(p > 0)then
        write(112,*)i,j,qt%qauxvec(p)
    else
        write(112,*)i,j,0.0
    endif
enddo
    write(112,*),""
enddo

melem = 500!qt%NO_VARIABLES
allocate(evals(melem))
allocate(elements(melem,2))
do i = 1, melem
    elements(i,:) = (/ i+1+200,i+1+200/)
enddo
SKMINV_DEBUG_LEVEL = 2
call selinv%invert(elements,evals)

qt%densities(1,:) = qt%qauxvec
do i = 1, melem
    qt%qauxvec(elements(i,1))     = -dble(evals(i))
    qt%densities(1,elements(i,1)) = -imag(evals(i))
enddo

do i = 1 , nx
do j = 1 , ny
    p = gindex(i,j)
    if( p > 0)then
        write(111,*),i,j,qt%qauxvec(p),qt%densities(1,p)
    else
        write(111,*),i,j,0.0,0.0
    endif
enddo
    write(111,*),""
enddo


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
!call system("cd "//output_folder//"; ./plot_bands.py")
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
        ! magnetic Field enters the hamiltonian by phase: EXP(i*DX*Bx*y)
        y = (atomA%atom_pos(2) - ny/2*dx ) * L2LA ! convert to atomic units
        coupling_val = -t0 * EXP(II*xdiff*a_dx*a_Bz*y)
    else if( xdiff ==  0 .and. abs(ydiff) == 1 ) then
        connect = .true.
        coupling_val = -t0
    endif


end function connect

end program transporter

!program p_Knitting
!use modunits
!use modalgs
!use modskminv
!use ifport
!implicit none
!
!
!integer,parameter :: nx = 1000 , no_elems = nx-60
!complex*16        :: vals(nx*nx),matC(nx,nx),evals(no_elems)
!integer           :: rc(nx*nx,2),elements(no_elems,2)
!integer :: i,j,iter
!doubleprecision :: t1, t2, sum_delta
!type(skminv) :: selinv
!
!matC = 0
!do i = 1 , nx
!    do j = i , i + 4
!        if( j <= nx ) then
!            matC(i,j) = 3*i
!            if(i /= j) matC(i,j) = +1 !+ rand() * II! + rand()*0.8
!            matC(j,i) = matC(i,j)
!            if(i /= j) matC(i,j) = +1 ! - II! + rand()*0.8
!            if(j-i == 3) then
!                matC(i,j) = 0
!                matC(j,i) = matC(i,j)
!            endif
!
!!            if(j-i == 5) then
!!                matC(i,j) = 0
!!                matC(j,i) = matC(i,j)
!!            endif
!        endif
!    enddo
!enddo
!!matC(1,1) =-5
!!matC(nx/2,nx/2+1) = 0.001
!iter = 0
!do i = 1, nx
!do j = 1, nx
!    if(abs(matC(i,j)) > 0.00001 ) then
!        iter = iter + 1
!        vals(iter) = matC(i,j)
!        rc(iter,1) = i
!        rc(iter,2) = j
!    endif
!enddo
!enddo
!
!call write_matrix(matC,"A.txt","(400f10.5)")
!
!t1 = get_clock()
!print*, "cond(A)=",alg_cond(matC)
!call inverse_matrix(nx,matC)
!t1 = get_clock() - t1
!
!t2 = get_clock()
!print*,"sel inv:"
!SKMINV_DEBUG_LEVEL = 2
!call selinv%init_memory(nx,iter,vals,rc)
!do i = 1 , no_elems
!    elements(i,:) = (/i+1+NINT(rand()*50),i+1+NINT(rand()*50)/)
!enddo
!call selinv%invert(elements,evals)
!
!t2 = get_clock() - t2
!sum_delta = 0
!!print*,"mat last:",matC(nx,nx)
!do i = 1 , size(evals)
!    sum_delta = sum_delta + abs(evals(i) - matC(elements(i,1),elements(i,2)))
!!    print"(i4,A,2f10.6,A,2f10.6,A,f10.6)",i," m=",matC(elements(i,1),elements(i,2))," v=",evals(i)," diff=",abs(evals(i)-matC(elements(i,1),elements(i,2)))
!enddo
!call selinv%free_memory()
!print*,"Time    inv:",t1
!print*,"Time selinv:",t2
!print*,"sun_delta  :",sum_delta
!
!
!contains
!
!
!subroutine write_matrix(matrix,filename, mformat)
!    complex*16 ,intent(in) :: matrix(:,:)
!    character(*) :: mformat,filename
!    integer :: i,j,sx,sy
!
!    sx = size(matrix,1)
!    sy = size(matrix,2)
!
!    open(unit=32123,file=filename)
!    do i = 1 , sx
!        write(32123, mformat),dble(matrix(i,:))
!    enddo
!    close(32123)
!
!end subroutine write_matrix
!
!
!subroutine print_matrix(matrix, mformat)
!    complex*16 ,intent(in) :: matrix(:,:)
!    character(*) :: mformat
!    integer :: i,j,sx,sy
!
!    sx = size(matrix,1)
!    sy = size(matrix,2)
!    do i = 1 , sx
!        print mformat,dble(matrix(i,:))
!    enddo
!
!end subroutine print_matrix
!
!end program p_Knitting
