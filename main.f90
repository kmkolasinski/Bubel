! ------------------------------------------------------ !
! Quantulaba - rashba.f90 - Krzysztof Kolasinski 2015
!
! ------------------------------------------------------ !
program rashba

use modunits    ! unit conversion tools
use modshape    ! leads shape definition
use modscatter  ! transport and eigen problem
use modinip     ! reading from ini file
implicit none
character(*),parameter     :: output_folder = "./"
type(qscatter)             :: qt
type(qshape)               :: lead_shape

! Parameters in SI units and atomic units
doubleprecision :: si_DX  , a_DX ! [nm]
doubleprecision :: si_Ef  , a_EF ! [meV]
doubleprecision :: si_Bz  , a_Bz ! [T]
doubleprecision :: si_Bx  , a_Bx ! [T]
doubleprecision :: si_By  , a_By ! [T]
doubleprecision :: si_a3D , a_a3D! [nm^2]
doubleprecision :: si_Fz  , a_Fz ! [meV/nm]
doubleprecision :: si_Rsb , a_Rsb! [meVnm] = si_Fz * si_a3D
doubleprecision :: si_Drs , a_Drs! [meVnm] => Dresselhaus
doubleprecision :: si_Vqpc1 , a_Vqpc1! [meV]
doubleprecision :: si_Vqpc2 , a_Vqpc2! [meV]

doubleprecision :: t0,tR,tD
doubleprecision,allocatable :: Apot(:,:,:) ! Potential Vector in atomic units
doubleprecision,allocatable :: Ugates(:,:)
doubleprecision            :: a_Emin,a_Emax
integer                    :: nx = 30
integer                    :: ny = 100
integer                    :: i,j,k
integer                    :: xin,xout,win,wout,qpcw
integer                    :: DETECTOR_ID
doubleprecision            :: lead_translation_vec(3),dx,xtip
character(300) :: line

! system definition 'hardcoded'
xin  = 80
xout = 20
win  = 25
wout = 100


!QSYS_FORCE_SCHUR_DECOMPOSITION = .true.
!QSYS_USE_ZGGEV_TO_FIND_MODES   = .true.
call read_config()
allocate(Apot  (nx,ny,2)) ! components of A = {Ax(x,y),Ay(x,y)}
allocate(Ugates(nx,ny))
Apot = 0
dx   = si_dx
qpcw = 100.0/dx
call convert_units()
call qt%init_system()
call calc_vector_potential()
call add_gates()
call create_system()

call qt%calculate_modes(a_Ef)

!a_Emin =-0.01 / E0 / 1000.0 ! converting from [meV] to atomic units
!a_Emax =1000.00 / E0 / 1000.0
!call qt%leads(1)%bands(output_folder//"bands.dat",-M_PI,+M_PI,M_PI/100.0,a_Emin,a_Emax)


!! Save lattice to file to see if constructed system is OK!
call qt%save_system(output_folder//"system.xml")

QSYS_DEBUG_LEVEL = 0 ! See more
open(unit=4321,file="T.txt")
!write(4321,"(A)"),"#Vqpc1 [meV]     Bz [T]     Rqpc1      Tqpc2        Ttotal (fake if used pseudo transparent gates)"
!do si_a3D = 0.0 , 1.0 , 0.01

!do xtip = 500.0 , 700.0 , 4.0

    call convert_units()
    call calc_vector_potential()
    call add_gates()
    call add_tip(xtip,400.0D0,50.0D0,50.0D0,1.5*a_EF)

    call qt%qsystem%update_lattice(c_matrix=coupling)
    call qt%calculate_modes(a_Ef)
    call qt%solve(1,a_Ef)

    print("(10e)"),si_Vqpc1,si_Bz,sum(qt%Rn(:)),qt%leads(1)%no_in_modes+0.0
    write(4321,"(10e)"),si_Vqpc1,si_Bz,xtip,sum(qt%Rn(:)),qt%leads(1)%no_in_modes+0.0
!enddo
close(4321)

! Save calculated electron density to file
do i = 1 , size(qt%qauxvec)
    qt%qauxvec(i) = sum(qt%densities(:,i))
enddo
call qt%qsystem%save_data(output_folder//"densities.xml",array2d=qt%densities,array1d=qt%qauxvec)
! Perform scan in function of Energy
! ----------------------------------------------------------
! X. Clean memory...
! ----------------------------------------------------------

call qt%destroy_system()
deallocate(Apot)
deallocate(Ugates)
!print*,"Generating plots..."
!print*,"Plotting band structure..."
!call system("cd "//output_folder//"; ./plot_bands.py")
!print*,"Plotting Transmission..."
!call system("cd "//output_folder//"; ./plot_T.py")
!print*,"Use Viewer program to see the structure and created leads."
contains

! -----------------------------------------------------------------
! Coupling function
! -----------------------------------------------------------------
logical function coupling(atomA,atomB,coupling_mat)
    use modcommons
    implicit none
    type(qatom) :: atomA,atomB
    complex*16  :: coupling_mat(:,:) ! you must overwrite this variable
    ! local variables
    integer         :: xdiff,ydiff,ix,iy
    doubleprecision :: dydiff,dxdiff,y
    complex*16      :: Cx,Cy

    ! Calculate distance between atoms in units of dx.
    dxdiff = (atomB%atom_pos(1)-atomA%atom_pos(1))/dx
    dydiff = (atomB%atom_pos(2)-atomA%atom_pos(2))/dx
    ! Convert it to integers
    xdiff = NINT(dxdiff)
    ydiff = NINT(dydiff)
    ix    = NINT(atomA%atom_pos(1)/dx)+1
    iy    = NINT(atomA%atom_pos(2)/dx)+1

    ! default return value
    coupling       = .false.
    coupling_mat   = 0.0

    Cx = EXP(II*a_DX*Apot(ix,iy,1))
    Cy = EXP(II*a_DX*Apot(ix,iy,2))

    if( xdiff == 0 .and. ydiff == 0 ) then
        coupling      = .true.
        coupling_mat = (4*t0 + Ugates(ix,iy) )* MAT_DIAG+ &
                    0.25*G_LAN*(a_Bx * MAT_SX + a_By * MAT_SY + a_Bz * MAT_SZ)
    else if( abs(xdiff) ==  1 .and. ydiff == 0 ) then
        coupling = .true.
        coupling_mat = -t0*Cx*MAT_DIAG + II*tR*Cx*MAT_SY - II*tD*Cx*MAT_SX
    else if( xdiff ==  0 .and. abs(ydiff) == 1 ) then
        coupling = .true.
        coupling_mat = -t0*Cy*MAT_DIAG - II*tR*Cy*MAT_SX + II*tD*Cy*MAT_SY
    endif
end function coupling

subroutine calc_vector_potential()
    integer :: i,j,n1
    doubleprecision :: x,y,w
    Apot = 0
    n1   = 30
    do i = 1 , nx
    do j = 1 , ny
        x = a_DX * i
        y = a_DX * j
        Apot(i,j,1) = -a_Bz*y
    enddo
    enddo
end subroutine calc_vector_potential

subroutine create_system()
integer :: i,j,nl,ny1,ny2
logical :: bSkip

! ----------------------------------------------------------
! 1. Create mesh - loop over width and height of the lattice
! ----------------------------------------------------------
do i = 1 , nx
do j = 1 , ny
    bSkip = .false.
    call qt%qatom%init((/ (i-1) * dx , (j-1) * dx + (i-1) * dx * 0.0 , 0.0 * dx /),no_in_states=2)
    !if(i <     xin .and. abs(j-ny/2)  > win ) bSkip = .true.
    !if(i >= nx-xout .and. abs(j-ny/2) > wout ) bSkip = .true.
!    if(i == nx-xout-1 .and. j == 1 )  bSkip = .true.
!    if(i == nx-xout-1 .and. j == ny ) bSkip = .true.
!
!    if(i == nx-xout-1 .and. j == ny/2-wout ) bSkip = .true.
!    if(i == nx-xout-1 .and. j == ny/2+wout ) bSkip = .true.



!    if(i == nx .and. j == 1  ) bSkip = .true.
!    if(i == nx .and. j == ny ) bSkip = .true.
    if(.not. bSkip) call qt%qsystem%add_atom(qt%qatom)
enddo
enddo
qt%qnnbparam%box = (/2*dx,2*dx,0.0D0/) ! do not search for the sites far than (-dx:+dx) direction
call qt%qsystem%make_lattice(qt%qnnbparam,c_matrix=coupling)
! ----------------------------------------------------------
! 2. Use generated mesh to calculate the band structure
! in the region of homogenous lead.
! ----------------------------------------------------------
call lead_shape%init_rect(SHAPE_RECTANGLE_XY,-0.5*dx,0.5*dx,-0.5*dx,(ny+1)*dx)
lead_translation_vec = (/  dx , 0.0D0 , 0.0D0 /)
call qt%add_lead(lead_shape,lead_translation_vec)


!call lead_shape%init_rect(SHAPE_RECTANGLE_XY,(xin-1.5)*dx,(-2.5+nx-xout)*dx,(-0.5)*dx,(+0.5)*dx)
call lead_shape%init_rect(SHAPE_RECTANGLE_XY,(xin-1.5)*dx,(-2.5+nx)*dx,(-0.5)*dx,(+0.5)*dx)
lead_translation_vec = (/  0.0D0, dx , 0.0D0 /)
!call qt%add_lead(lead_shape,lead_translation_vec)
!qt%leads(2)%lead_type         = QSYS_LEAD_TYPE_PSEUDO_TRANSPARENT
!qt%leads(2)%pseudo_lead_phase = exp(II*sqrt(2*a_Ef*M_EFF)*a_DX)

call lead_shape%init_rect(SHAPE_RECTANGLE_XY,(xin-1.5)*dx,(-2.5+nx)*dx,(ny-1.5)*dx,(ny-0.5)*dx)
!call lead_shape%init_rect(SHAPE_RECTANGLE_XY,(xin-1.5)*dx,(-2.5+nx-xout)*dx,(ny-1.5)*dx,(ny-0.5)*dx)
lead_translation_vec = (/  0.0D0, -dx , 0.0D0 /)
!call qt%add_lead(lead_shape,lead_translation_vec)
!qt%leads(3)%lead_type         = QSYS_LEAD_TYPE_PSEUDO_TRANSPARENT
!qt%leads(3)%pseudo_lead_phase = exp(II*sqrt(2*a_Ef*M_EFF)*a_DX)

call lead_shape%init_rect(SHAPE_RECTANGLE_XY,(nx-xout-2.5)*dx,(-1.5+nx-xout)*dx,(0.5)*dx,(ny/2-wout-1.5)*dx)
lead_translation_vec = (/  -dx , 0.0D0 , 0.0D0 /)
!call qt%add_lead(lead_shape,lead_translation_vec)
!
!qt%leads(4)%lead_type         = QSYS_LEAD_TYPE_PSEUDO_TRANSPARENT
!qt%leads(4)%pseudo_lead_phase = exp(II*sqrt(2*a_Ef*M_EFF)*a_DX)

call lead_shape%init_rect(SHAPE_RECTANGLE_XY,(nx-xout-2.5)*dx,(-1.5+nx-xout)*dx,(ny/2+wout-0.5)*dx,(ny-1.5)*dx)
lead_translation_vec = (/  -dx , 0.0D0 , 0.0D0 /)
!call qt%add_lead(lead_shape,lead_translation_vec)
!
!qt%leads(5)%lead_type         = QSYS_LEAD_TYPE_PSEUDO_TRANSPARENT
!qt%leads(5)%pseudo_lead_phase = exp(II*sqrt(2*a_Ef*M_EFF)*a_DX)

call lead_shape%init_rect(SHAPE_RECTANGLE_XY,(nx-1.5)*dx,(nx-0.5)*dx,-0.5*dx,(ny)*dx)
!call lead_shape%init_rect(SHAPE_RECTANGLE_XY,(nx-1.5)*dx,(nx-0.5)*dx,0.5*dx,(ny-1.5)*dx)
lead_translation_vec = (/  -dx , 0.0D0 , 0.0D0 /)
call qt%add_lead(lead_shape,lead_translation_vec)

DETECTOR_ID = qt%no_leads

end subroutine create_system


subroutine add_gates()
    integer :: i,j
    doubleprecision :: x,y
    Ugates = 0
    do i = 1 , nx
    do j = 1 , ny
        x = dx * i
        y = dx * j
        Ugates(i,j) = &
                    + finite_rect(x,y,(xin-qpcw/2.0-qpcw/2.0)*dx,(xin+qpcw/2.0-qpcw/2.0)*dx,-100*dx,(ny/2-15)*dx,10.0D0,a_Vqpc1) &
                    + finite_rect(x,y,(xin-qpcw/2.0-qpcw/2.0)*dx,(xin+qpcw/2.0-qpcw/2.0)*dx,(ny/2+15)*dx,(ny+100)*dx,10.0D0,a_Vqpc1)
!        write(444,*),i,j,Ugates(i,j)*E0*1000.0
    enddo
!        write(444,*),""
    enddo
    Ugates(1,:) = Ugates(2,:)

end subroutine add_gates


subroutine add_tip(xt,yt,dxt,dyt,Utip)
    doubleprecision :: xt,yt,dxt,dyt,Utip
    integer :: i,j
    doubleprecision :: x,y

    do i = 5 , nx-5
    do j = 5 , ny-5
        x = dx * i
        y = dx * j
        Ugates(i,j) = Ugates(i,j) + gauss_gate(Utip,xt,yt,dxt,dyt,x,y)

    enddo
    enddo

end subroutine add_tip
! --------------------------------------------------------
! Reads system parameters from config.ini file
! --------------------------------------------------------
subroutine read_config()
    call setIniFilename("config.ini")
    print*,"Reading config file:"
    call getDoubleValue("scattering","m_eff",M_EFF)
    call getDoubleValue("scattering","eps"  ,E_MAT)
    call getDoubleValue("scattering","g_lan",G_LAN)
    call getIntValue   ("scattering","nx",nx)
    call getIntValue   ("scattering","ny",ny)
    call getDoubleValue("scattering","dx",si_DX)
    call getDoubleValue("scattering","Ef",si_Ef)
    call getDoubleValue("scattering","Bz",si_Bz)
    call getDoubleValue("scattering","Bx",si_Bx)
    call getDoubleValue("scattering","By",si_By)
    call getDoubleValue("scattering","so_alpha3D",si_a3D)
    call getDoubleValue("scattering","so_Fz",si_Fz)
    call getDoubleValue("scattering","so_Drs",si_Drs)
    call getDoubleValue("scattering","V_qpc1",si_Vqpc1)
    call getDoubleValue("scattering","V_qpc2",si_Vqpc2)
    call modunits_set_semiconductor_effective_mass_units(M_EFF,E_MAT)
endsubroutine read_config

! ------------------------------------------------------
! Gate potential implementation from paper:
! Modeling the patterned two-dimensional electron gas: Electrostatics
! ------------------------------------------------------
doubleprecision function finite_rect(x,y,L,R,B,T,d,Vg)
    doubleprecision :: x,y,L,R,T,B,d,Vg
    finite_rect = Vg*(g_aux(x-L,y-B,d) + g_aux(x-L,T-y,d) + &
                  g_aux(R-x,y-B,d) + g_aux(R-x,T-y,d))
endfunction finite_rect
doubleprecision function g_aux(u,v,d)
    doubleprecision :: u,v,d
    doubleprecision :: a1,R
    R  = sqrt(u**2 + v**2 + d**2)
    a1 = atan(u*v/(R*d))
    g_aux = 1.0 / 2.0 / M_PI * a1
end function g_aux


double precision function gauss_gate(U0,xp,yp,sigmax,sigmay,x,y) result(rval)
        doubleprecision, intent(in):: U0,xp,yp,sigmax,sigmay,x,y

        rval =  U0 * exp( -(( x - xp)/(2*sigmax))**2 ) * exp( -(( y - yp)/(2*sigmay))**2 )

endfunction gauss_gate

! --------------------------------------------------------
! Converts parameters from SI {nm,meV,T} to atomic units
! --------------------------------------------------------
subroutine convert_units()
    a_DX    = si_DX * L2LA
    a_EF    = si_Ef / E0 / 1000.0
    a_Vqpc1 = si_Vqpc1 / E0 / 1000.0
    a_Vqpc2 = si_Vqpc2 / E0 / 1000.0
    a_Bz    = BtoAtomicB(si_Bz)
    a_By    = BtoAtomicB(si_By)
    a_Bx    = BtoAtomicB(si_Bx)
    a_a3D   = si_a3D * L2LA**2
    a_Fz    = si_Fz  / L2LA / E0 / 1000.0 ! meV/nm => 1/E0 * 1/L2LA
    a_Rsb   = a_a3D*a_Fz
    a_Drs   = si_Drs * L2LA / E0 / 1000.0 ! meVnm => 1/A0 * 1/E0
    ! other variables
    t0 = 1.0/(2*M_EFF*a_DX**2)
    tR = a_Rsb/(2.0*a_DX)
    tD = a_Drs/(2.0*a_DX)
!    print*,"Rashba      coupling:",tR
!    print*,"Dresselhaus coupling:",tD
end subroutine convert_units
end program rashba
