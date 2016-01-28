! ------------------------------------------------------ !
! Quantulaba - recursive_wfm.f90 - Krzysztof Kolasinski 2016
! This is an advanced example on how Quantulaba may be use to solve
! scattering problem using recursive formulation of WFM.
! This algorithm is very similar to well known RGF. Here
! we consider the channel with two barriers inside in presence
! of external magnetic field. In recursive version of WFM
! channel is divided into slices and scattering wave function
! is restored by iterative inversion of Green like matrix.
! After recursion process one may calculate Tranmission matrix
! and scattering wave functions (after inverse recursion loop).
! Here we compare standard WFM with recursive one.
! ------------------------------------------------------ !

program recursive_wfm

use modunits   ! unit conversion tools
use modscatter ! eigen values and transport
use modlead    ! bandgap structure
use modshape
use modalgs
implicit none
character(*),parameter     :: output_folder = "recursive_wfm_output/"
type(qscatter)             :: qt
type(qshape)               :: lead_shape

doubleprecision            :: a_dx,a_Bz ! dx and magnetic field in atomic units
integer ,parameter         :: nx = 150
integer ,parameter         :: ny = 50
doubleprecision,parameter  :: dx = 2.5  ! grid spacing [nm]
integer , dimension(nx,ny) :: gindex    ! converts local index (i,j) to global index
doubleprecision            :: Vbarriers(nx,ny) ! matrix with potential inside the channel

doubleprecision            :: lead_translation_vec(3),Ef
doubleprecision            :: T_WFM,T_RECURSIVE ! Transmissions obtained with difference methods
integer :: i,j,k


call qt%init_system()
QSYS_USE_ZGGEV_TO_FIND_MODES = .true.  ! in case of magnetic field some times one must force to
                                       ! solve generalized eigenvalue problem -> singulat hoping matrix

! Use atomic units in effective band model -> see modunit.f90 for more details
call modunits_set_GaAs_params()
a_dx = dx * L2LA         ! convert it to atomic units
a_Bz = BtoAtomicB(0.5D0) ! convert magnetic field from Tesla SI unit to atomic units

! ----------------------------------------------------------
! 1. Create mesh - loop over width and height of the lattice
! ----------------------------------------------------------
Vbarriers = 0
k      = 0
gindex = 0
do i = 1 , nx
do j = 1 , ny
    call qt%qatom%init((/ (i-1) * dx , (j-1) * dx , 0.0 * dx /))
    call qt%qsystem%add_atom(qt%qatom)
    k           = k + 1
    gindex(i,j) = k

    ! Create double potential barrier for electron
    if( i == nx/2-30 ) Vbarriers(i,j) = 3.0/E0/1000.0
    if( i == nx/2+30 ) Vbarriers(i,j) = 3.0/E0/1000.0
enddo
enddo
! ----------------------------------------------------------
! 2. Construct logical connections between sites on the mesh.
! ----------------------------------------------------------
qt%qnnbparam%box = (/2*dx,2*dx,0.0D0/)
call qt%qsystem%make_lattice(qt%qnnbparam,c_simple=connect)
! ----------------------------------------------------------
! 3. Use generated mesh to calculate the band structure
! in the region of homogenous lead.
! ----------------------------------------------------------
! Setup shape and initialize lead
call lead_shape%init_rect(SHAPE_RECTANGLE_XY,-0.5*dx,0.5*dx,-0.5*dx,(ny+1)*dx)
lead_translation_vec = (/  dx , 0.0D0 , 0.0D0 /)
call qt%add_lead(lead_shape,lead_translation_vec)

! Add second lead at the end of the system
call lead_shape%init_rect(SHAPE_RECTANGLE_XY,(-0.5+nx-1)*dx,(0.5+nx-1)*dx,-0.5*dx,(ny+1)*dx)
lead_translation_vec = (/  -dx , 0.0D0 , 0.0D0 /)
call qt%add_lead(lead_shape,lead_translation_vec)


QSYS_DEBUG_LEVEL = 0 ! do not show all debug messages

! Perform loop over Fermi energies
open(unit=222,file=output_folder//"T.dat")
do Ef = .2/E0/1000.0 ,  2.0/E0/1000.0 ,  2.0/E0/1000.0/40.0
    ! Calculate lead self energies
    call qt%calculate_modes(Ef)

    call calc_recursive(Ef) ! Calculate T from recursive algorithm
    call qt%solve(1,Ef)     ! Calculate T from standard WFM
    T_WFM = sum(qt%Tn)

    print*,Ef*E0*1000,"T=",T_WFM
    write(222,*),Ef*E0*1000,T_WFM,T_RECURSIVE
enddo

! ----------------------------------------------------------
! X. Clean memory...
! ----------------------------------------------------------
call lead_shape%destroy_shape()
call qt%destroy_system()
print*,"Generating plots..."
print*,"Plotting Transmission..."
call system("cd "//output_folder//"; ./plot_T.py")

contains


! ---------------------------------------------------------------------------
! This function decides if site A (called here atomA)  has hoping
! to atom B, and what is the value of the coupling constant.
! If there is no interaction between them returns false, otherwise true.
! ---------------------------------------------------------------------------
logical function connect(atomA,atomB,coupling_val)
    use modcommons
    implicit none
    type(qatom) :: atomA,atomB
    complex*16  :: coupling_val ! you must overwrite this variable

    ! local variables
    integer         :: xdiff,ydiff,ix,iy
    doubleprecision :: dydiff,dxdiff,t0,y,Vpot

    ! Calculate distance between atoms in units of dx.
    dxdiff = (atomA%atom_pos(1)-atomB%atom_pos(1))/dx
    dydiff = (atomA%atom_pos(2)-atomB%atom_pos(2))/dx

    ix = nint(atomA%atom_pos(1)/dx)+1
    iy = nint(atomA%atom_pos(2)/dx)+1
    ! Convert it to integers
    xdiff = NINT(dxdiff)
    ydiff = NINT(dydiff)

    ! default return value
    connect = .false.

    ! hoping parameter
    t0 = 1/(2*m_eff*a_dx**2)

    if( xdiff == 0 .and. ydiff == 0 ) then
        connect      = .true.
        coupling_val = 4*t0 + Vbarriers(ix,iy)
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

! ---------------------------------------------------------------------------
! This is simple application of recursive WFM method. Here we allocate all
! the arrays statically since considered system is small enough. Additionally,
! all multiplications are done with matmul() function which is considerably
! slow in comparision to optimized version of gemm from BLAS library. However
! it can be easily optimized. The system is divided into slices:
!   tau * c(k-1,:) + (Ef - H_k ) * c(k,:) + tau^Dagg * c(k,:) = 0
! knowing that
!   tau * c(Nx-1,:) + (Ef - H_Nx - SigmaN ) * c(Nx,:) = 0
! one may derive recurtion relation for other slices, connected by G(i,:,:) matrix
!   c(i,:) = G(i,:,:) * c(i-1,:)
! Once we have all G(i,:,:) matrices we can calculate scattering wave function Y(:,:)
! inside device. From this we calculate Transmission T, however in order to calculate
! T one does not have to have whole Y(:,:), same as in case RGF.
! ---------------------------------------------------------------------------
subroutine calc_recursive(pEf)
    doubleprecision :: pEf
    complex*16 :: H0(ny,ny) , Tau(ny,ny) , TauDag(ny,ny)  , Gk(nx,ny,ny)
    complex*16 :: ck(nx,ny) , Sigma1(ny,ny) , SigmaN(ny,ny) , ctmp(ny,ny) , Gamma1(ny)
    complex*16 :: vphi(nx*ny)
    doubleprecision :: Diag(ny,ny),Vpot,rho(nx,ny),modes(qt%leads(1)%no_in_modes,nx,ny),Tk
    doubleprecision :: total_Tn , total_Rn , Tn(qt%leads(1)%no_in_modes) , Rn(qt%leads(1)%no_in_modes)
    doubleprecision :: gamma(nx,ny)
    integer :: i,j,k,m,leadID

    Gk = 0
    Diag = 0
    do i = 1 , ny
        Diag(i,i) = 1
    enddo
    ! Initialize necessary matrices from lead
    Tau    = conjg(transpose(qt%leads(1)%valsTau)) ! coupling matrix
    TauDag = qt%leads(1)%valsTau    ! hermitian cojugate of coupling matrix
    H0     = qt%leads(1)%valsH0     ! Hamiltonian matrix of slice
    Sigma1 = qt%leads(1)%SigmaMat - qt%leads(1)%valsH0 ! Self energy, here we remove H0 since SigmaMat already contains it.
    SigmaN = qt%leads(2)%SigmaMat - qt%leads(2)%valsH0 ! Self energy of second lead.
    leadID = 1

    ! Calculation of block greens functions, we start from k=N to k=1
    ! k == 1 and k == N requires special treatement.

    ctmp = pEf*Diag - ( H0 + SigmaN )
    call inverse_matrix(ny,ctmp)
    Gk(Nx,:,:) = matmul(ctmp,Tau)

    do k = nx - 1 , 2 , -1
        ctmp = pEf*Diag - ( H0 + Diag * Vbarriers(k,1) + matmul(TauDag,Gk(k+1,:,:)) )
        call inverse_matrix(ny,ctmp)
        Gk(k,:,:) = matmul(ctmp,Tau)
    enddo

    ctmp = pEf*Diag - ( H0 + Sigma1 + matmul(TauDag,Gk(2,:,:)) )
    call inverse_matrix(ny,ctmp)
    Gk(1,:,:) = ctmp


    rho     = 0
    modes   = 0
    Tk      = 0
    ! Inverse recursion, for each incoming mode calculate scattering wave function
    do m = 1 , qt%leads(1)%no_in_modes
        Gamma1 = 0
        ck     = 0
        ! First slice k==1 is treated differently
        call calc_gamma(m,1,Gamma1)
        ck(1,:) = matmul(Gk(1,:,:),Gamma1)
        do k = 2 , nx
            ck(k,:) = matmul(Gk(k,:,:),ck(k-1,:))
        enddo
        ! Fill calculated scattering matrix for 1D array
        do i = 1 , nx
        do j = 1 , ny
            vphi(gindex(i,j)) = ck(i,j)
        enddo
        enddo
        ! ----------------------------------------------------
        !               Calculate transmission from
        !                      wave function
        ! ----------------------------------------------------
        total_Tn = 0
        total_Rn = 0
        do i = 1 , qt%no_leads
            qt%leads(i)%totalT = 0
            ! Skip tranmission from input lead
            if( i /= leadID) then
                call qt%leads(i)%calculate_Tnm(qt%qsystem%atoms,vphi)
                total_Tn = total_Tn  + qt%leads(i)%modeT
            else
                call qt%leads(i)%calculate_Tnm(qt%qsystem%atoms,vphi,m)
                total_Rn = total_Rn  + qt%leads(i)%modeT
            endif
            qt%leads(i)%totalT = qt%leads(i)%totalT + qt%leads(i)%modeT
        enddo

        ! Calculate transmission and reflection probabilities
        Tn(m) = total_Tn
        Rn(m) = total_Rn
        ! ----------------------------------------------------
        ! Scattering modes and LDOS is calculated
        modes(m,:,:) = abs(ck)**2!
        rho = rho + abs(ck)**2
    enddo
    T_RECURSIVE = sum(Tn)

end subroutine calc_recursive

! Caclulate the left hand side vector present in WFM method.
subroutine calc_gamma(modin,leadID,GammaL)
    integer :: modin,leadID
    complex*16 :: GammaL(ny),cval
    integer :: i,j

    cval    = qt%leads(leadID)%lambdas(M_IN,modin)**(-1)
    do i = 1 , ny
        do j = 1 , ny
            GammaL(i) = GammaL(i) + (-qt%leads(leadID)%LambdaMat(i,j) + cval*conjg(qt%leads(leadID)%valsTau(j,i)) ) &
                                    *qt%leads(leadID)%modes(M_IN,modin,j)
        enddo
    enddo
end subroutine calc_gamma
end program recursive_wfm
