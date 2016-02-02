MODULE modunits
    USE, INTRINSIC :: ISO_C_BINDING
    use modcommons
    implicit none


    double precision,parameter  :: M_PI      = 3.1415926
    double precision,parameter  :: A0        = 0.0529177249 ! Bohr radius in [nm]
    double precision,parameter  :: E0        = 27.211384523 ! Atomic unit of energy in [eV]


    double precision  :: KbT
    double precision  :: Rd    ! For donors units
    double precision  :: L2LA  ! Conversion from nm to atomic unit of lenght
    double precision  :: LA2L  ! Conversion from atomic unit of lenght to nm
    double precision  :: L2LR  ! Same but for Rydberg units
    double precision  :: LR2L
    double precision  :: M_EFF ! effective mass for seminconductor systems
    double precision  :: E_MAT ! dielectric constant
    double precision  :: G_LAN ! lande factor
    double precision  :: so_alpha3D  ! SO rashba constant [nm^2]
    double precision  :: so_Fz       ! gradient of electric field in z plane [meV/nm ]
    double precision  :: so_rashba   ! rashba coeficient [nm*meV ]

    complex*16 ,parameter :: MAT_SX  (2,2) = (/ (/ 0.0  , 1.0 /) , (/ 1.0 , 0.0  /) /)
    complex*16 ,parameter :: MAT_SY  (2,2) = (/ (/ 0*II , -II /) , (/ II  , 0*II /) /)
    complex*16 ,parameter :: MAT_SZ  (2,2) = (/ (/ 1.0  , 0.0 /)    , (/ 0.0 , -1.0   /) /)
    doubleprecision ,parameter :: MAT_DIAG(2,2) = (/ (/ 1.0  , 0.0 /)    , (/ 0.0 , 1.0  /) /)



    contains
    subroutine modunits_set_GaAs_params()
        ! based on -- http://www.ioffe.ru/SVA/NSM/Semicond/GaInAs/
        ! lande factor -- http://arxiv.org/pdf/1406.2848v1.pdf
        ! http://journals.aps.org/prb/abstract/10.1103/PhysRevB.35.7729

        G_LAN      =  +0.8 ! very small value
        so_Fz      =   0.0 ! in [meV/nm]
        so_alpha3D =   0.0 ! in [nm^2]
        call modunits_set_semiconductor_effective_mass_units(0.067D0,12.0D0)
    end subroutine modunits_set_GaAs_params

    subroutine modunits_set_InGaAs_params()
        ! based on -- http://www.ioffe.ru/SVA/NSM/Semicond/GaInAs/
        ! lande factor -- http://arxiv.org/pdf/1406.2848v1.pdf
        ! http://journals.aps.org/prb/abstract/10.1103/PhysRevB.35.7729

        G_LAN      =  -8.97
        so_Fz      =  20.0   ! in [meV/nm]
        so_alpha3D =   0.572 ! in [nm^2]

        call modunits_set_semiconductor_effective_mass_units(0.0465D0,12.0D0)
    end subroutine modunits_set_InGaAs_params

    subroutine modunits_set_InSb_params()
        ! na podstawie pracy doktorskiej (colwiz zakladka SO)
        G_LAN = -50.0
        call modunits_set_semiconductor_effective_mass_units(0.014D0,16.0D0)
    end subroutine modunits_set_InSb_params

    subroutine modunits_set_semiconductor_effective_mass_units(pM_EFF,pE_MAT)
        double precision,intent(in) :: pM_EFF
        double precision,intent(in) :: pE_MAT

        E_MAT     = pE_MAT
        M_EFF     = pM_EFF
        Rd        = E0*M_EFF/E_MAT/E_MAT
        L2LR      = M_EFF/E_MAT/A0
        LR2L      = 1.0/L2LR
        L2LA      = 1.0/A0
        LA2L      = 1.0/L2LA

        so_rashba = so_alpha3D*so_Fz

    end subroutine modunits_set_semiconductor_effective_mass_units



    ! Convert magnetic field from [T] to [donor T] and from [donor T] to [T]
    double precision function BtoDonorB(inBZ) result(rval)
      double precision, intent(in) :: inBZ
      rval = 2*inBZ/m_eff*5.78838e-5/Rd
    endfunction BtoDonorB
    double precision function DonorBtoB(inBZ) result(rval)
      double precision, intent(in) :: inBZ
      rval = inBZ*m_eff/5.78838e-5*Rd/2
    endfunction DonorBtoB

    ! Convert magnetic field from [T] to [atomic T] and from [atomic T] to [T]
    double precision function BtoAtomicB(inBZ) result(rval)
      double precision, intent(in) :: inBZ
      rval = 2*inBZ*5.78838e-5/E0
    endfunction BtoAtomicB
    double precision function AtomicBtoB(inBZ) result(rval)
      double precision, intent(in) :: inBZ
      rval = inBZ/5.78838e-5*E0/2
    endfunction AtomicBtoB

END MODULE modunits
