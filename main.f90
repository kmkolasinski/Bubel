
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
 use modalgs
 implicit none
 character(*),parameter :: output_folder = "plots/"
 type(qscatter)             :: qt
 type(qshape)               :: lead_shape
 doubleprecision            :: zeros_vector(200)
 doubleprecision            :: a_dx,a_Emin,a_Emax,a_Bz
 integer                    :: no_expected_states
 integer ,parameter         :: nx = 200
 integer ,parameter         :: ny = 25
 doubleprecision,parameter  :: dx = 5.0 ! [nm]
 integer , dimension(nx,ny) :: gindex ! converts local index (i,j) to global index
 integer :: i,j,k,p,ytip,xtip
 doubleprecision            :: lead_translation_vec(3),Ef
 doubleprecision            :: Vpot(nx*ny),Vtip,T1,T2
 complex*16,allocatable     :: H0(:,:),G0(:,:),waveR(:,:),waveRT(:,:),waveL(:,:),waveN(:,:)
 complex*16,allocatable     :: invGamma(:,:),dSkm(:,:),tmpSkm(:,:)
 complex*16 :: L(ny,ny),QdagUdag(ny,ny),tmpGU(ny,ny),Umodes(ny,ny),tmp2GU(ny,ny)
 type(qsmatrix),dimension(2,2) :: smatrix0
 ! Use atomic units in effective band model -> see modunit.f90 for more details
 call modunits_set_GaAs_params()
 a_dx = dx * L2LA ! convert it to atomic units
 a_Bz = BtoAtomicB(0.5D0) ! 1 Tesla to atomic units

 ! Initalize system
 call qt%init_system()
 QSYS_FORCE_SCHUR_DECOMPOSITION = .true.
 ! ----------------------------------------------------------
 ! 1. Create mesh - loop over width and height of the lattice
 ! ----------------------------------------------------------
 k      = 0
 gindex = 0
 do i = 1 , nx
 do j = 1 , ny

    call qt%qatom%init((/ (i-1) * dx , (j-1) * dx , 0.0 * dx /),no_in_states=1)
!    if( i == nx/2 .and. j == ny/3) then
!       gindex(i,j) = 1
!    else
    ! Add atom wavefuncto the system.
    call qt%qsystem%add_atom(qt%qatom)
        k           = k + 1
        gindex(i,j) = k
!    endif
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
QSYS_DEBUG_LEVEL = 2
p = 1

call qtools_fd_template(v_dnFdxn,0,1)

call qtools_fd_template(v_dnFdxn,1,p)
call qtools_fd_template(v_dnFdyn,1,p)

call qtools_fd_template(v_dnFdxn,2,p)
call qtools_fd_template(v_dnFdyn,2,p)

Vpot = 0

call qt%qsystem%make_lattice(qt%qnnbparam,c_matrix=coupling)



! ----------------------------------------------------------
! 4. Use generated mesh to calculate the band structure
! in the region of homogenous lead.
! ----------------------------------------------------------
! Setup shape and initialize lead
! Lead needs to know where it is (lead_shape) and using this information it will
! create propper matrices and connections using list of atoms
lead_translation_vec = (/  p*dx , 0.0D0 , 0.0D0 /)
call lead_shape%init_range_3d((/-0.5*dx,0.0D0,0.0D0/),lead_translation_vec)
call qt%add_lead(lead_shape,lead_translation_vec)

! Add second lead at the end of the system
lead_translation_vec = (/  -p*dx , 0.0D0 , 0.0D0 /)
call lead_shape%init_range_3d((/(nx-0.5)*dx,0.0D0,0.0D0/),lead_translation_vec)
call qt%add_lead(lead_shape,lead_translation_vec)


call qt%save_system(output_folder//"system.xml")

! Solve scattering problem for Ef=0.001
Ef = 10/E0/1000.0
QSYS_DEBUG_LEVEL = 2 ! show more info


allocate(H0(nx*ny,nx*ny))
allocate(G0(nx*ny,nx*ny))

call qt%qsystem%update_lattice(c_matrix=coupling)
call qt%calculate_modes(Ef)
call qt%solve(1,Ef)

allocate(waveL(qt%leads(1)%no_in_modes,qt%NO_VARIABLES))
allocate(waveRT(qt%leads(2)%no_in_modes,qt%NO_VARIABLES))
allocate(waveR(qt%leads(2)%no_in_modes,qt%NO_VARIABLES))
allocate(waveN(qt%leads(2)%no_in_modes,qt%NO_VARIABLES))
allocate(invGamma(qt%leads(2)%no_sites,qt%leads(2)%no_sites))
allocate(dSkm(qt%leads(2)%no_sites,qt%leads(2)%no_in_modes))
allocate(tmpSkm(qt%leads(2)%no_sites,qt%leads(2)%no_in_modes))

print*,"Hamiltonian size:",size(H0,1)," mem:",size(H0,1)**2/1024.0**2/8.0*2,"[MB]"
H0   = 0
do i = 1 , qt%NO_NON_ZERO_VALUES
    j = qt%ROWCOLID(i,1)
    k = qt%ROWCOLID(i,2)
    H0(j,k) = qt%MATHVALS(i)
!    print*,i,j,k,dble(H0(j,k))
enddo

G0    = -H0
print*,"Finding Greens function: G0"
call inverse_matrix(nx*ny,G0)
waveL = qt%wavefunc

Umodes = -transpose(qt%leads(2)%modes(M_OUT,:,:))
!call inverse_matrix(qt%leads(2)%no_sites,Umodes)

call qt%copy_smatrix(smatrix0)
a_Bz = -a_Bz
call qt%qsystem%update_lattice(c_matrix=coupling)
call qt%calculate_modes(Ef)
a_Bz = -a_Bz
call qt%solve(2,Ef)
waveR = qt%wavefunc

!call write_matrix(H0,output_folder//"H0.dat","(50000e12.3)")
!call write_matrix(G0,output_folder//"G0.dat","(50000e12.3)")

! Save calculated electron density to file
open(433,file=output_folder//"ldos.dat")
do i = 1 , nx
do j = 1 , ny
    k = gindex(i,j)
    write(433,"(2i,100e20.6)"),i,j,sum( abs(waveL(:,k))**2 ),  sum( abs(waveR(:,k))**2 ), &
                                  sum( abs(waveL(:,k))**2 ) + sum( abs(waveR(:,k))**2 ),abs(waveL(:,k))**2
enddo
    write(433,*),""
enddo
close(433)

invGamma = -qt%leads(2)%GammaVecs
print*,"sum(GammaVec)=",sum(invGamma)

call inverse_matrix(qt%leads(2)%no_sites,invGamma)
call write_matrix(qt%leads(2)%GammaVecs,output_folder//"GammaVec.dat","(500f10.4)")
print*,"cond(invGamma)=",alg_cond(invGamma)



print*,shape(-qt%leads(2)%GammaVecs)
print*,shape(L)
print*,shape(QdagUdag)


tmp2GU = -transpose(qt%leads(2)%modes(M_OUT,:,:))
call gemm(qt%leads(2)%GammaVecs,Umodes,tmpGU)



call zalg_SVD_QL(tmpGU,L,QdagUdag)
print*,"cond(L)       =",alg_cond(L)
print*,"cond(L(modes))=",alg_cond(L(1:qt%leads(2)%no_in_modes,1:qt%leads(2)%no_in_modes))
print*,"cond(QdagUdag)=",alg_cond(QdagUdag)
!
!call gemm(QdagUdag,L,tmp2GU,transa="C")

!call gemm(tmp2GU,tmpGU,QdagUdag)
!print*,"sum(GammaVec2)=",sum(QdagUdag)
!stop
!call write_matrix(L,output_folder//"L.dat","(500f10.4)")

call qt%qsystem%update_lattice(c_matrix=coupling)
call qt%calculate_modes(Ef)

call qt%solve(2,Ef)
waveRT = qt%wavefunc


print*,"T=",sum(abs(smatrix0(1,2)%Tnm)**2),"R=",sum(abs(smatrix0(1,1)%Tnm)**2)
print*,"W=",sum(abs(smatrix0(1,2)%Tnm)**2)+sum(abs(smatrix0(1,1)%Tnm)**2)


QSYS_DEBUG_LEVEL = 0
open(unit=111,file=output_folder//"sgmT.dat")
open(unit=112,file=output_folder//"sgmP1.dat")
xtip = nx/2-5
do ytip = 1 , ny
do xtip = nx/2-nx/4 , nx/2+nx/4
    Vpot = 0
    Vtip = Ef/1.0
!    T1 = get_clock()
    Vpot(gindex(xtip,ytip)) = Vtip
!    call qt%qsystem%update_lattice(c_matrix=coupling)
!    call qt%calculate_modes(Ef)
!    call qt%solve(1,Ef)
!    write(111,"(100e20.6)"),ytip*dx,sum(qt%Tn(:))
!    print*,"time full:",get_clock()-T1
    T1 = get_clock()
    call calc_sgm()
!    print*,"time fast:",get_clock()-T1
enddo
    write(112,*),""
enddo
close(111)



!open(unit=111,file=output_folder//"T.dat")
!print*,"Performing energy scan..."
!do Ef = 0.0 , 0.006 , 0.00005
!    call qt%qsystem%update_lattice(c_matrix=coupling)
!    call qt%calculate_modes(Ef)
!    call qt%solve(1,Ef)
!    write(111,"(100f20.6)"),Ef,sum(qt%Tn(:))
!enddo
!close(111)


! ----------------------------------------------------------
! X. Clean memory...
! ----------------------------------------------------------
call lead_shape%destroy_shape()
call qt%destroy_system()

contains
subroutine calc_sgm()
    doubleprecision :: total_Tn , total_Rn , T0 , ldos
    integer         :: i,k , m , j , leadID, la,ls,lg,p,no_modes
    complex*16      :: tmpT,tmpTn(500)
    leadID = 1
    p      = gindex(xtip,ytip)
    ldos   = 0.5*(sum(abs(waveL(:,p))**2)+sum(abs(waveRT(:,p))**2))

    tmpSkm = 0
    do m = 1 , qt%leads(1)%no_in_modes
    do k = 1 , qt%leads(2)%no_in_modes
        tmpSkm(k,m) = 0
        do j = 1 , qt%NO_VARIABLES
        tmpSkm(k,m) = tmpSkm(k,m) - &
!          (waveR(k,j)*Vpot(j)/(1+G0(p,p)*Vpot(j))*waveL(m,j))
          (waveR(k,j)*Vpot(j)/(1+cmplx(ldos,ldos)*Vpot(j))*waveL(m,j))
!          ((waveR(k,j))*Vpot(j)/(1+0*G0(p,p)*Vpot(j))*waveL(m,j))
        enddo
    enddo
    enddo
!    call gemm(invGamma,tmpSkm,dSkm)



!    do m = 1 , qt%leads(1)%no_in_modes
!    do k = 1 , qt%leads(2)%no_sites
!        la = qt%leads(2)%l2g(k,1)
!        ls = qt%leads(2)%l2g(k,2)
!        lg = qt%qsystem%atoms(la)%globalIDs(ls);
!        tmpSkm(k,m) = 0
!        do j = 1 , qt%NO_VARIABLES
!        tmpSkm(k,m) = tmpSkm(k,m) - &
!!          (G0(lg,j)*Vpot(j)*waveL(m,j))
!!          (G0(lg,j)*Vpot(j)/(1+cmplx(ldos,ldos)*Vpot(j))*waveL(m,j))
!          (G0(lg,j)*Vpot(j)/(1+G0(p,p)*Vpot(j))*waveL(m,j))
!        enddo
!    enddo
!    enddo

!    call gemm(invGamma,tmpSkm,dSkm)
!    print"(A,6f12.4)","G(p,p)=",G0(p,p),G0(p,p)/ldos

!
!    waveN = 0
!    do m = 1 , qt%leads(1)%no_in_modes
!        do i = 1 , qt%leads(2)%no_sites
!            la = qt%leads(2)%l2g(i,1)
!            ls = qt%leads(2)%l2g(i,2)
!            lg = qt%qsystem%atoms(la)%globalIDs(ls);
!            waveN(m,lg) = dSkm(i,m)
!!            waveN(m,lg) = tmpSkm(i,m)
!        enddo
!    enddo
!    waveN = waveN + waveL
!
!    total_Tn = 0
!    total_Rn = 0
!    qt%Tn    = 0
!    qt%Rn    = 0
!    do m = 1 , qt%leads(1)%no_in_modes
!
!    do i = 1 , qt%no_leads
!        qt%leads(i)%totalT = 0
!        ! Skip tranmission from input lead
!        if( i /= leadID) then
!            call qt%leads(i)%calculate_Tnm(qt%qsystem%atoms,waveN(m,:))
!            total_Tn = total_Tn  + qt%leads(i)%modeT
!        else
!            call qt%leads(i)%calculate_Tnm(qt%qsystem%atoms,waveN(m,:),m)
!            total_Rn = total_Rn  + qt%leads(i)%modeT
!        endif
!        qt%leads(i)%totalT = qt%leads(i)%totalT + qt%leads(i)%modeT
!    enddo
!        qt%Tn(m) = total_Tn
!        qt%Rn(m) = total_Rn
!    enddo
!!    write(112,"(100e20.6)"),ytip*dx,total_Tn
!!    print"(100e20.6)",ytip*dx,total_Tn

    ! calculate QdagUdag * dS
    no_modes = qt%leads(2)%no_in_modes
    call gemm(QdagUdag,tmpSkm,dSkm)

    total_Tn = 0
    do m = 1 , qt%leads(1)%no_in_modes

        tmpTn           = 0
        qt%leads(2)%Tnm = 0
        do i = 1 , qt%leads(2)%no_in_modes
            tmpT = 0
            do j = 1 , i-1
                tmpT = tmpT + tmpTn(j)*L(i,j)
            enddo
            tmpTn(i) =  (dSkm(i,m) - tmpT)/L(i,i)
            qt%leads(2)%Tnm(i) = tmpTn(i) + smatrix0(1,2)%Tnm(m,i)
!            print*,m,i,tmpTn(i)
        enddo

        total_Tn = total_Tn  + sum(abs(qt%leads(2)%Tnm(:))**2)
    enddo
!    print"(100e20.6)",ytip*dx,total_Tn
    write(112,"(100e20.6)"),xtip*dx,ytip*dx,total_Tn



end subroutine calc_sgm

logical function coupling(atomA,atomB,coupling_mat)
    use modcommons
    implicit none
    type(qatom) :: atomA,atomB
    complex*16  :: coupling_mat(:,:) ! you must overwrite this variable
    ! local variables
    integer         :: xdiff,ydiff,idx,idy,ix,iy
    doubleprecision :: dydiff,dxdiff,t0,y,rs,bz,Vext
    ! Calculate distance between atoms in units of dx.
    dxdiff = (atomA%atom_pos(1)-atomB%atom_pos(1))/dx
    dydiff = (atomA%atom_pos(2)-atomB%atom_pos(2))/dx
    ! Convert it to integers
    idx = NINT(dxdiff)
    idy = NINT(dydiff)
    xdiff = NINT(dxdiff)
    ydiff = NINT(dydiff)
    ix  = NINT(atomA%atom_pos(1)/dx)+1
    iy = NINT(atomA%atom_pos(2)/dx)+1
    ! default return value
    coupling       = .true.
    coupling_mat   =  0.00

    t0 = 1/(2*m_eff*a_dx**2)
    Vext = 0
    if(ix == nx/2 .and. iy == ny/3) Vext = 2*t0

!    coupling_mat = -t0*(dnFdxn(2,idx,idy,0) + dnFdyn(2,idx,idy,0)) &
!    + dnFdxn(0,idx,idy,0)*(Vpot(gindex(ix,iy)) + Vext )
!    if(sum(abs(coupling_mat)) < qsys_double_error) coupling = .false.


    coupling       = .false.
    if( xdiff == 0 .and. ydiff == 0 ) then
        coupling      = .true.
        coupling_mat = 4*t0 + Vpot(gindex(ix,iy)) + Vext
    else if( abs(xdiff) ==  1 .and. ydiff == 0 ) then
        coupling = .true.
        ! magnetic Field enters the hamiltonian by phase: EXP(i*DX*Bx*y)
        y = (atomA%atom_pos(2) - ny/2*dx ) * L2LA ! convert to atomic units
        coupling_mat = -t0 * EXP(II*xdiff*a_dx*a_Bz*y)
    else if( xdiff ==  0 .and. abs(ydiff) == 1 ) then
        coupling = .true.
        coupling_mat = -t0
    endif


end function coupling

subroutine write_matrix(matrix,filename, mformat)
    complex*16 ,intent(in) :: matrix(:,:)
    character(*) :: mformat,filename
    integer :: i,j,sx,sy

    sx = size(matrix,1)
    sy = size(matrix,2)

    open(unit=32123,file=filename)
    do i = 1 , sx
        write(32123, mformat),abs(matrix(i,:))
    enddo
    close(32123)

end subroutine write_matrix

end program transporter

