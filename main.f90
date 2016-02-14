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
 use modalgs
 implicit none
 character(*),parameter :: output_folder = "plots/"
 type(qscatter)             :: qt
 type(qshape)               :: lead_shape
 doubleprecision            :: zeros_vector(200)
 doubleprecision            :: a_dx,a_Emin,a_Emax,a_Bz
 integer                    :: no_expected_states
 integer ,parameter         :: nx = 100
 integer ,parameter         :: ny = 1
 doubleprecision,parameter  :: dx = 5.0 ! [nm]
 integer , dimension(nx,ny) :: gindex ! converts local index (i,j) to global index
 integer :: i,j,k,xtip
 doubleprecision            :: lead_translation_vec(3),Ef,Vtip,Vbar
 complex*16 :: H0(nx,nx),G0(nx,nx)


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
Ef = 4.0/E0/1000.0
Vtip = Ef/2
Vbar = 1
QSYS_DEBUG_LEVEL = 1 ! show more info

call qt%calculate_modes(Ef)

call qt%solve(1,Ef)
! Save calculated electron density to file
do i = 1 , size(qt%qauxvec)
    qt%qauxvec(i) = sum(qt%densities(:,i))
enddo
call qt%qsystem%save_data(output_folder//"densities.xml",array2d=qt%densities,array1d=qt%qauxvec)
! Perform scan in function of Energy


!open(unit=111,file=output_folder//"Gx.dat")
!print*,"Performing energy scan..."
!QSYS_DEBUG_LEVEL = 1 ! show more info
!do xtip = nx/2-40 , nx/2+40
!    ! Update hamiltonian elemenents value
!    call qt%qsystem%update_lattice(c_simple=connect)
!    call qt%solve(1,Ef)
!    print*,xtip
!    write(111,"(100f20.6)"),xtip*dx,sum(qt%Tn(:)),sum(qt%Rn(:))
!enddo
!close(111)
!
!call calc_sgm()


Vtip = 0
open(unit=111,file=output_folder//"T.dat")
open(unit=112,file=output_folder//"Tn.dat")
print*,"Performing energy scan..."
QSYS_DEBUG_LEVEL = 0 ! show more info
do Ef = 0.0000001 , 5.0/E0/1000 , 0.01/E0/1000
    ! Update hamiltonian elemenents value
    Vbar = 1
    call qt%qsystem%update_lattice(c_simple=connect)
    call qt%calculate_modes(Ef)
    call qt%solve(1,Ef)
    write(111,"(100e20.6)"),Ef*1000.0*E0,sum(qt%Tn(:))

    call calc_double_barier()

enddo
close(111)


! ----------------------------------------------------------
! X. Clean memory...
! ----------------------------------------------------------
call lead_shape%destroy_shape()
call qt%destroy_system()
!print*,"Generating plots..."
!print*,"Plotting band structure..."
!call system("cd "//output_folder//"; python plot_bands.py")
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
    doubleprecision :: dydiff,dxdiff,t0,y,vpot
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
        vpot = 0
        if( NINT(atomA%atom_pos(1)/dx)+1 == nx/2 - 20 ) vpot = t0/2*Vbar
        if( NINT(atomA%atom_pos(1)/dx)+1 == nx/2 + 20 ) vpot = t0/2*Vbar
        if( NINT(atomA%atom_pos(1)/dx)+1 == xtip )      vpot = vpot + Vtip
        coupling_val = 2*t0 + vpot
    else if( abs(xdiff) ==  1 .and. ydiff == 0 ) then
        connect      = .true.
        coupling_val = -t0
    endif
end function connect


subroutine calc_sgm()

    doubleprecision :: total_Tn , total_Rn
    integer         :: i , leadID
    complex*16      :: vphi0(nx),dvphi(nx),Gnp,Gpp,Vp,Yp,dYn
    leadID = 1
    print*,"SGM simple method:"
    ! ----------------------------------------------------
    !               Calculate transmission from
    !                      wave function
    ! ----------------------------------------------------
    Vtip = 0
    call qt%qsystem%update_lattice(c_simple=connect)
    call qt%construct_matrix(Ef)
    H0   = 0
    do i = 1 , qt%NO_NON_ZERO_VALUES
        j = qt%ROWCOLID(i,1)
        k = qt%ROWCOLID(i,2)
        H0(j,k) = qt%MATHVALS(i)
    enddo
    G0 = -H0
    call inverse_matrix(nx,G0)
    call qt%solve(1,Ef)
    vphi0 = qt%wavefunc(1,:)
    dvphi = 0
!    print*,"cond(H0)=",alg_cond(H0)
!    call write_matrix(H0,output_folder//"H0.dat","(500f10.4)")
!    call write_matrix(G0,output_folder//"G0.dat","(500f10.4)")
    Vtip = Ef/2
open(unit=111,file=output_folder//"Gxp.dat")
do xtip = nx/2-40 , nx/2+40
    ! Update hamiltonian elemenents value
    Gnp = G0(nx,xtip)
    Gpp = G0(xtip,xtip)
    Vp  = Vtip
    Yp  = vphi0(xtip)
    dYn = -Gnp*Vp*Yp + Gnp*Vp*(1+Gpp*Vp)**(-1)*Gpp*Vp*Yp
    dvphi     = vphi0
    dvphi(nx) = dvphi(nx) + dYn

!    print*,xtip
!    write(111,"(100f20.6)"),xtip*dx,sum(qt%Tn(:))


    total_Tn = 0
    total_Rn = 0
    do i = 1 , qt%no_leads
        qt%leads(i)%totalT = 0
        ! Skip tranmission from input lead
        if( i /= leadID) then
            call qt%leads(i)%calculate_Tnm(qt%qsystem%atoms,dvphi)
            total_Tn = total_Tn  + qt%leads(i)%modeT
        else
            call qt%leads(i)%calculate_Tnm(qt%qsystem%atoms,dvphi,1)
            total_Rn = total_Rn  + qt%leads(i)%modeT
        endif
        qt%leads(i)%totalT = qt%leads(i)%totalT + qt%leads(i)%modeT
    enddo
    write(111,"(100f20.6)"),(xtip)*dx,total_Tn,total_Rn
enddo
end subroutine calc_sgm


subroutine calc_double_barier()

    doubleprecision :: total_Tn , total_Rn , t0
    integer         :: i , leadID , pnts(2) , p2,p1
    complex*16      :: vphi0(nx),dvphi(nx),Gnp(1,2),Gpp(2,2),Vp(2,2),Yp(2,1)
    complex*16      :: Ypn(2,1),dYn,P(2,2),Pn(2,2),OP(2,2),L(1,2),S1,S2
    leadID = 1

    ! ----------------------------------------------------
    !               Calculate transmission from
    !                      wave function
    ! ----------------------------------------------------
    Vtip = 0
    Vbar = 0
    t0 = 1/(2*m_eff*a_dx**2)
    call qt%qsystem%update_lattice(c_simple=connect)
    call qt%solve(1,Ef)
    H0   = 0
    do i = 1 , qt%NO_NON_ZERO_VALUES
        j = qt%ROWCOLID(i,1)
        k = qt%ROWCOLID(i,2)
        H0(j,k) = qt%MATHVALS(i)
    enddo
    G0    = -H0
    call inverse_matrix(nx,G0)
!    call write_matrix(G0,output_folder//"G0.dat","(500f10.4)")
    vphi0 = qt%wavefunc(1,:)
    dvphi = 0
    Vtip  = 0

    pnts = (/ nx/2 - 20 , nx/2 + 20 /)
    Vp = 0
    do p1 = 1 , 2
        Gnp(1,p1)  = G0(nx,pnts(p1))
        Vp(p1,p1)  = t0/2
        Yp(p1,1)   = vphi0(pnts(p1))

    do p2 = 1 , 2
        Gpp(p1,p2) = G0(pnts(p1),pnts(p2))
!        if(p1 /= p2) Gpp(p1,p2) = Gpp(p1,p2)

    enddo
    enddo

!    Vp(1,1) = Vp(1,1)/100

    ! Update hamiltonian elemenents value

    call gemm(Gpp,Vp,P)
    call gemm(Gnp,Vp,L)
    OP = P
    do p1 = 1 , 2
        OP(p1,p1) = OP(p1,p1) + 1
    enddo

    call inverse_matrix(2,OP)
    call gemm(OP,P,Pn)
    OP = Pn
    do p1 = 1 , 2
        OP(p1,p1) = OP(p1,p1) - 1
    enddo

    call gemm(OP,Yp,Ypn)

    dYn = sum(L(1,:)*Ypn(:,1))

    S1 = - Gnp(1,1)*Vp(1,1)/(1 + Vp(1,1)*Gpp(1,1))*Yp(1,1)
    S2 = - Gnp(2,2)*Vp(2,2)/(1 + Vp(2,2)*Gpp(2,2))*Yp(2,2)


    dvphi     = vphi0
    dvphi(nx) = dvphi(nx) + dYn


    total_Tn = 0
    total_Rn = 0
    do i = 1 , qt%no_leads
        qt%leads(i)%totalT = 0
        ! Skip tranmission from input lead
        if( i /= leadID) then
            call qt%leads(i)%calculate_Tnm(qt%qsystem%atoms,dvphi)
            total_Tn = total_Tn  + qt%leads(i)%modeT
        else
            call qt%leads(i)%calculate_Tnm(qt%qsystem%atoms,dvphi,1)
            total_Rn = total_Rn  + qt%leads(i)%modeT
        endif
        qt%leads(i)%totalT = qt%leads(i)%totalT + qt%leads(i)%modeT
    enddo
    write(112,"(100e20.6)"),Ef*1000*E0,total_Tn,total_Rn

end subroutine calc_double_barier

subroutine print_matrix(matrix, mformat)
    complex*16 ,intent(in) :: matrix(:,:)
    character(*) :: mformat
    integer :: i,j,sx,sy

    sx = size(matrix,1)
    sy = size(matrix,2)
    print*,"no rows:",sx
    print*,"no cols:",sy
    print*,"format :",mformat
    do i = 1 , sx
        print mformat,abs(matrix(i,:))
    enddo

end subroutine print_matrix

subroutine write_matrix(matrix,filename, mformat)
    complex*16 ,intent(in) :: matrix(:,:)
    character(*) :: mformat,filename
    integer :: i,j,sx,sy

    sx = size(matrix,1)
    sy = size(matrix,2)
    print*,"no rows:",sx
    print*,"no cols:",sy
    print*,"format :",mformat
    open(unit=32123,file=filename)
    do i = 1 , sx
        write(32123, mformat),dble(matrix(i,:))
    enddo
    close(32123)

end subroutine write_matrix

end program transporter
