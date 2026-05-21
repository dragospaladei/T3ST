!========================================================================================================
! File: Run_parameters.f90
! Module kept as constants for compatibility with the rest of the code.
! Owns shared constants, run parameters, derived normalization values,
! simulation-matrix import, and EFIT parameter post-processing.
!
! Purpose: define main physical and numerical constants of the model
!
! Record of revisions:
!    Date       | Programmer          | Description
!  -----------  | ------------------- | -----------------------------------------
!  --/--/2021   | D. I. Palade        | Initial version
!  --/--/2023   | L. M. Pomarjanschi  | Corrections and development
!  01/01/2024   | D. I. Palade        | Remake
!  --/--/2026   | D. I. Palade        | Continuous improvements
! Coordinate convention:
!   (x,y,z) = (q1,q2,q3) are field-aligned turbulence coordinates.
!   (X,Y,Z) are real-space coordinates, with Z the toroidal angle.
!========================================================================================================

MODULE constants
   USE sims
   IMPLICIT NONE
   SAVE   ! all data are preserved in different procedures

   !=====================================================================================================
   ! PART 1: MODULE DATA
   !=====================================================================================================

   !---------------------------------------------------------------------------------
   ! Universal constants
   !---------------------------------------------------------------------------------
   REAL(KIND=rp), PARAMETER :: pi   = 3.1415926545_rp         ! pi
   REAL(KIND=rp), PARAMETER :: E    = 2.7182818284_rp         ! Euler's number
   REAL(KIND=rp), PARAMETER :: q0   = 1.60217*10.0**(-19.0)   ! elementary charge
   REAL(KIND=rp), PARAMETER :: mi   = 1.67262*10.0**(-27.0)   ! proton mass
   REAL(KIND=rp), PARAMETER :: me   = 9.10938*10.0**(-31.0)   ! electron mass
   REAL(KIND=rp), PARAMETER :: eps  = 8.85418*10.0**(-12.0)   ! vacuum dielectric constant
   REAL(KIND=rp), PARAMETER :: mu0  = 1.25663*10.0**(-6.0)    ! vacuum magnetic permeability
   REAL(KIND=rp), PARAMETER :: kB   = 1.38064*10.0**(-23.0)   ! Boltzmann constant
   REAL(KIND=rp), PARAMETER :: logL = 17.0                    ! Coulomb logarithm
   INTEGER, PARAMETER       :: ON = 1
   INTEGER, PARAMETER       :: OFF = 0
   REAL(KIND=rp), PARAMETER :: eps_xi    = 1.0e-15_rp
   REAL(KIND=rp), PARAMETER :: eps_rr    = 1.0e-15_rp
   REAL(KIND=rp), PARAMETER :: eps_root  = 1.0e-15_rp
   REAL(KIND=rp), PARAMETER :: eps_B     = 1.0e-15_rp
   REAL(KIND=rp), PARAMETER :: eps_Bsp   = 1.0e-15_rp
   REAL(KIND=rp), PARAMETER :: eps_omega = 1.0e-15_rp
   REAL(KIND=rp), PARAMETER :: eps_large = 1.0e-6_rp


   !---------------------------------------------------------------------------------
   ! Switches
   !---------------------------------------------------------------------------------
   ! INTEGER :: Method                              ! legacy full-f/delta-f switch
   INTEGER       :: USE_larmor                     ! FLR effects: 0=OFF, 1=ON
   INTEGER       :: USE_coll                       ! collisions: 0=OFF, 1=ON
   INTEGER       :: USE_turb                       ! turbulence: 0=OFF, 1=ON
   INTEGER       :: USE_magnturb                   ! RMP: 0=OFF, 1=ON
   INTEGER       :: USE_freq                       ! mode freqs: 0=frozen, 1=real frequencies
   INTEGER       :: USE_polar                      ! polarization drift: 0=OFF, 1=ON
   INTEGER       :: USE_PC                         ! initial condition constancy psi/Pc
   INTEGER       :: USE_real                       ! same/different turbulence realizations across particles
   INTEGER       :: USE_corr                       ! compute/export Lagrangian correlations
   INTEGER       :: USE_balloon                    ! 0=homogeneous kz, 1=ballooning g(z)
   INTEGER       :: USE_tilt                       ! tilting switch
   INTEGER       :: USE_testing                    ! testing switch

   !---------------------------------------------------------------------------------
   ! Numerical parameters
   !---------------------------------------------------------------------------------
   REAL(KIND=rp) :: t0                             ! simulation start time        [R0/vth]
   REAL(KIND=rp) :: tc                             ! collisions start time        [R0/vth]
   REAL(KIND=rp) :: tt                             ! turbulence start time        [R0/vth]
   REAL(KIND=rp) :: tmax                           ! simulation end time          [R0/vth]

   INTEGER       :: Nt                             ! number of time steps
   INTEGER       :: Nreal                          ! realizations of turbulent field
   INTEGER       :: Nc                             ! number of waves / modes
   INTEGER       :: Nloop                          ! simulation loops (super-ensemble)
   INTEGER       :: Np                             ! particles per loop
   INTEGER       :: ntraj                          ! number of trajectories exported

   INTEGER       :: Nci                            ! number of ITG partial-waves
   INTEGER       :: Nce                            ! number of TEM partial-waves

   !---------------------------------------------------------------------------------
   ! Plasma / normalization
   !---------------------------------------------------------------------------------
   REAL(KIND=rp) :: Ti                             ! ion temperature @ r00
   REAL(KIND=rp) :: Te                             ! electron temperature @ r00
   REAL(KIND=rp) :: Zeff                           ! effective charge number
   REAL(KIND=rp) :: Aeff                           ! effective mass number
   REAL(KIND=rp) :: ndens                          ! ion density @ r00             [m^-3] (imported in 1e19 m^-3)

   REAL(KIND=rp) :: vth                            ! thermal velocity (H ions)     [m/s] @ r00
   REAL(KIND=rp) :: rhoi                           ! Larmor radius (H ions)        [m]   @ r00 (B0)
   REAL(KIND=rp) :: wi                             ! Larmor frequency              [Hz]
   REAL(KIND=rp) :: c0                             ! collision-related constant
   REAL(KIND=rp) :: delta                          ! collision-related constant

   !---------------------------------------------------------------------------------
   ! Gradients / drives
   !---------------------------------------------------------------------------------
   REAL(KIND=rp) :: Ln                             ! density gradient length       [R0*d(lnn)/dr]
   REAL(KIND=rp) :: Li                             ! ion temperature gradient      [R0*d(lnTi)/dr]
   REAL(KIND=rp) :: Le                             ! electron temperature gradient [R0*d(lnTe)/dr]

   !---------------------------------------------------------------------------------
   ! Geometry / equilibrium
   !---------------------------------------------------------------------------------
   INTEGER       :: magnetic_model                 ! equilibrium: 1-EFIT, 2-EFIT-sa, 3-Solovev, 4-Circular

   REAL(KIND=rp) :: B0                             ! magnetic field @ axis         [T]
   REAL(KIND=rp) :: R0                             ! major radius                  [m]
   REAL(KIND=rp) :: a0                             ! minor radius                  [m]

   REAL(KIND=rp) :: s1                             ! circular q-profile coefficient
   REAL(KIND=rp) :: s2                             ! circular q-profile coefficient
   REAL(KIND=rp) :: s3                             ! circular q-profile coefficient

   REAL(KIND=rp) :: amp                            ! Solovev flux amplitude
   REAL(KIND=rp) :: elong                          ! elongation

   INTEGER       :: device                         ! (1=MAST, 2=TCV, 3=WEST, 4=JTSA)
   INTEGER       :: shot                           ! shot number
   INTEGER       :: shotslice                      ! time (ms) for equilibrium slice

   REAL(KIND=rp) :: Omgt0                          ! toroidal angular frequency    [1e3 Hz] @ r00
   REAL(KIND=rp) :: Omgtprim                       ! d ln Omgt / d rhot            [1] @ r00

   REAL(KIND=rp) :: q00                            ! safety factor @ r00
   REAL(KIND=rp) :: r00                            ! reference surface (for global evaluation)
   REAL(KIND=rp) :: q10                            ! reference surface (for global evaluation)
   REAL(KIND=rp) :: psi0                           ! poloidal flux @ r0

   ! Extra equilibrium / imported-data helpers
   REAL(KIND=rp) :: alfa                           ! Solovev parameter (separatrix/ellipticity-related)
   REAL(KIND=rp) :: gama                           ! Solovev parameter (sqrt(gama)=separatrix position)
   CHARACTER(LEN=:), ALLOCATABLE            :: efit_file
   REAL(KIND=rp), ALLOCATABLE, DIMENSION(:) :: Efit_data
   INTEGER                                  :: NgridR, NgridZ, Nqua, Ngrid
   REAL(KIND=rp)                            :: minR, maxR, minZ, maxZ
   REAL(KIND=rp)                            :: stepR, stepZ
   REAL(KIND=rp)                            :: triang

   !---------------------------------------------------------------------------------
   ! Turbulence model
   !---------------------------------------------------------------------------------
   INTEGER       :: turb_model                     ! type of turbulence model: 1=old, 2=new
   INTEGER       :: x_corr                         ! x correlation: 1-gaussian, 2-exponential, 3-GENE-like
   INTEGER       :: y_corr                         ! y correlation: 1-Madi, 2-GENE
   INTEGER       :: t_corr                         ! t correlation: 1-gaussian, 2-exponential

   REAL(KIND=rp) :: Phi                            ! turbulence amplitude [ePhi/Ti]
   REAL(KIND=rp) :: turbprof                       ! exponent for synthetic turb. profile
   REAL(KIND=rp) :: Ai                             ! ITG fraction

   REAL(KIND=rp) :: lambdax                        ! x correlation length [rhoi units]
   REAL(KIND=rp) :: lambday                        ! y correlation length [rhoi units]
   REAL(KIND=rp) :: lambdaz                        ! z correlation length due to decoherence [q0*R0 units]
   REAL(KIND=rp) :: lbalonz                        ! z correlation length due to ballooning envelope [q0*R0 units]
   REAL(KIND=rp) :: tauc                           ! time correlation [R0/vth]
   REAL(KIND=rp) :: k0i                            ! ITG dominant ky [1/rhoi]
   REAL(KIND=rp) :: k0e                            ! TEM dominant ky [1/rhoi]
   REAL(KIND=rp) :: gamma_ZF                       ! zonal-flow shearing rate
   REAL(KIND=rp) :: gamma_E                        ! rotational shearing rate

   REAL(KIND=rp) :: Ae                             ! TEM fraction
   INTEGER       :: dmmax                          ! maximal parallel number

   !---------------------------------------------------------------------------------
   ! Particle initial conditions / markers
   !---------------------------------------------------------------------------------
   REAL(KIND=rp) :: X0                             ! starting position R
   REAL(KIND=rp) :: Y0                             ! starting position Z
   REAL(KIND=rp) :: Z0                             ! starting position phi

   REAL(KIND=rp) :: Ts                             ! energy average [Ti units] @ r00
   REAL(KIND=rp) :: Es                             ! average ion energy [Ti units]
   REAL(KIND=rp) :: pitch                          ! average pitch angle
   REAL(KIND=rp) :: As                             ! ion mass number
   REAL(KIND=rp) :: Zs                             ! ionization state
   REAL(KIND=rp) :: taucc                          ! inverse collision frequency [R0/vth]
   REAL(KIND=rp) :: Lns                            ! density gradient for the species
   REAL(KIND=rp) :: Lts                            ! temperature gradient of the species
   REAL(KIND=rp) :: annulus_width                  ! full radial width for local-annulus initialization

   INTEGER       :: position_type                  ! init position: 1-fixed point, ...
   INTEGER       :: pitch_type                     ! init pitch: 1-fixed, 2-random uniform
   INTEGER       :: energy_type                    ! init energy: 1-fixed, 2-Boltzmann

   REAL(KIND=rp) :: usetilt, balloon, noballoon, norm ! helpers

   !---------------------------------------------------------------------------------
   ! Misc. derived helpers
   !---------------------------------------------------------------------------------
   REAL(KIND=rp) :: C1                             ! q1 = C1*Q1
   REAL(KIND=rp) :: C2                             ! q2 = C2*Q2
   REAL(KIND=rp) :: C3                             ! q3 = C3*Q3

CONTAINS

   !=====================================================================================================
   ! PART 2: SUBROUTINE load_run_parameters(run)
   !
   ! Load one row from array1 into named module variables, then compute all derived
   ! normalization constants used by the rest of the simulation.
   !
   ! This should be called once at the beginning of every run.
   !=====================================================================================================

   SUBROUTINE load_run_parameters(run)
      USE sims
      IMPLICIT NONE

      INTEGER, INTENT(IN)             :: run      ! row index in Sim_XXX.dat / DB_XXX.dat
      REAL(KIND=rp), DIMENSION(cols1) :: pp       ! parameter vector for this run

      pp = array1(run, :)

      !---------------------------------------------------------------------------------
      ! Import raw parameters from the selected run row.
      !
      ! Name-based indexing keeps this code robust if the GUI reorders columns.
      ! Integer parameters are converted explicitly; physical quantities remain real.
      !---------------------------------------------------------------------------------
      ! Method         = int(pp(parameter_index("Method")))
      USE_larmor     = int(pp(parameter_index("USE_larmor")))
      USE_coll       = int(pp(parameter_index("USE_coll")))
      USE_turb       = int(pp(parameter_index("USE_turb")))
      USE_magnturb   = int(pp(parameter_index("USE_magnturb")))
      USE_freq       = int(pp(parameter_index("USE_freq")))
      USE_polar      = int(pp(parameter_index("USE_polar")))
      USE_PC         = int(pp(parameter_index("USE_PC")))
      USE_real       = int(pp(parameter_index("USE_real")))
      USE_corr       = int(pp(parameter_index("USE_corr")))
      USE_balloon    = int(pp(parameter_index("USE_balloon")))
      USE_tilt       = int(pp(parameter_index("USE_tilt")))
      USE_testing    = int(pp(parameter_index("USE_testing")))

      t0             = pp(parameter_index("t0"))
      tc             = pp(parameter_index("tc"))
      tt             = pp(parameter_index("tt"))
      tmax           = pp(parameter_index("tmax"))
      Nt             = int(pp(parameter_index("Nt")))
      Nreal          = int(pp(parameter_index("Nreal")))
      Nc             = int(pp(parameter_index("Nc")))
      Nloop          = int(pp(parameter_index("Nloop")))
      Np             = int(pp(parameter_index("Np")))
      ntraj          = int(pp(parameter_index("ntraj")))

      Ti             = pp(parameter_index("Ti"))
      Te             = pp(parameter_index("Te"))
      Zeff           = pp(parameter_index("Zeff"))
      Aeff           = pp(parameter_index("Aeff"))
      ndens          = pp(parameter_index("ndens"))

      Ln             = pp(parameter_index("Ln"))
      Li             = pp(parameter_index("Li"))
      Le             = pp(parameter_index("Le"))

      magnetic_model = int(pp(parameter_index("magnetic_model")))
      B0             = pp(parameter_index("B0"))
      R0             = pp(parameter_index("R0"))
      a0             = pp(parameter_index("a0"))
      s1             = pp(parameter_index("s1"))
      s2             = pp(parameter_index("s2"))
      s3             = pp(parameter_index("s3"))
      amp            = pp(parameter_index("amp"))
      elong          = pp(parameter_index("elong"))
      device         = int(pp(parameter_index("device")))
      shot           = int(pp(parameter_index("shot")))
      shotslice      = int(pp(parameter_index("shotslice")))

      Omgt0          = pp(parameter_index("Omgt0"))
      Omgtprim       = pp(parameter_index("Omgtprim"))

      turb_model     = int(pp(parameter_index("turb_model")))
      x_corr         = int(pp(parameter_index("x_corr")))
      y_corr         = int(pp(parameter_index("y_corr")))
      t_corr         = int(pp(parameter_index("t_corr")))
      Phi            = pp(parameter_index("Phi"))
      turbprof       = pp(parameter_index("turbprof"))
      Ai             = pp(parameter_index("Ai"))
      lambdax        = pp(parameter_index("lambdax"))
      lambday        = pp(parameter_index("lambday"))
      lambdaz        = pp(parameter_index("lambdaz"))
      lbalonz        = pp(parameter_index("lbalonz"))
      tauc           = pp(parameter_index("tauc"))
      k0i            = pp(parameter_index("k0i"))
      k0e            = pp(parameter_index("k0e"))
      gamma_ZF       = pp(parameter_index("gamma_ZF"))

      X0             = pp(parameter_index("X0"))
      Y0             = pp(parameter_index("Y0"))
      Z0             = pp(parameter_index("Z0"))
      Ts             = pp(parameter_index("Ts"))
      Es             = pp(parameter_index("Es"))
      pitch          = pp(parameter_index("pitch"))
      As             = pp(parameter_index("As"))
      Zs             = pp(parameter_index("Zs"))
      position_type  = int(pp(parameter_index("position_type")))

      annulus_width = pp(parameter_index("annulus_width"))
      pitch_type     = int(pp(parameter_index("pitch_type")))
      energy_type    = int(pp(parameter_index("energy_type")))
      Lns            = pp(parameter_index("Lns"))
      Lts            = pp(parameter_index("Lts"))

      CALL validate_raw_run_parameters

      
      Nqua = 16   ! number of quantities stored on the EFIT-like grid
      CALL load_efit_equilibrium()

      !---------------------------------------------------------------------------------
      ! Convert imported units and evaluate derived normalization quantities.
      !---------------------------------------------------------------------------------
      Ti    = Ti*1000.0_rp                     ! imported in keV
      Te    = Te*1000.0_rp                     ! imported in keV
      ndens = ndens*10.0_rp**19                ! imported in 1e19 m^-3

      vth  = sqrt(Ti*q0/mi)
      rhoi = mi*vth/q0/abs(B0)
      wi   = q0*abs(B0)/mi

      r00 = sqrt((X0 - R0)**2 + Y0**2)

      gama = 1.0_rp - 2.0_rp*sqrt((a0/R0)**2 - (a0/R0)**4)
      alfa = elong*(sqrt(2.0_rp - gama) - sqrt(gama)) / sqrt(2.0_rp - 2.0_rp*gama)
      amp  = 0.5_rp*(1.0_rp + alfa**2) / (1.0_rp - gama) / s1

      a0    = a0/R0
      Omgt0 = Omgt0*10.0_rp**3 / (vth/R0)

      r00 = r00/R0
      X0  = X0/R0
      Y0  = Y0/R0

      CALL OnlyQ(X0, Y0, q00, psi0)
      CALL Onlycoord_new(X0, Y0, Z0, q10) ! reference q1 for delta_q1

      C1 = a0*R0/rhoi
      C2 = -r00/q00 * R0/rhoi
      C3 = 1.0_rp

      Nci = int(Nc*Ai**2)
      Nce = Nc - Nci

      Ae = sqrt(1.0_rp - Ai**2)

      c0 = logL*ndens/8.0_rp/pi * (Zeff*Zs*q0**2/(eps*mi*As))**2

      delta = 5.0_rp*4.0_rp/3.0_rp/sqrt(pi) * (tmax - t0)/Nt * R0*2.0_rp*logL*ndens*(Zeff*Zs*q0**2)**2 / &
              (pi*eps**2*mi**2*As**2*vth**4*sqrt(2.0_rp/Aeff)**3)
      taucc = 0.0_rp

      CALL validate_derived_run_parameters

      gamma_E = gamma_ZF + Omgtprim          ! total shearing rate used by turbulence
      usetilt = real(USE_tilt, rp)
      balloon = real(USE_balloon, rp)        ! integer switch -> 0.0 or 1.0
      noballoon = 1.0_rp - balloon
      norm = 1.0_rp / sqrt(real(Nc, rp)/2.0_rp)

      ! dmmax = 2
      dmmax = 2 + int(sqrt(-2.0_rp*log(0.05_rp)/lbalonz))

   END SUBROUTINE load_run_parameters

   SUBROUTINE load_efit_equilibrium()
      USE sims, ONLY: project_root
      IMPLICIT NONE

      CHARACTER(LEN=35) :: aux
      INTEGER           :: i

      IF ((magnetic_model /= 1) .AND. (magnetic_model /= 2)) RETURN

      WRITE(aux, '(I0)') shot

      ! Device-specific EFIT folder.
      IF (device == 1) THEN
         efit_file = project_root//'/G_EQDSK/MAST-U/'//TRIM(aux)//'/QQ_'
      ELSEIF (device == 2) THEN
         efit_file = project_root//'/G_EQDSK/TCV/'//TRIM(aux)//'/QQ_'
      ELSEIF (device == 3) THEN
         efit_file = project_root//'/G_EQDSK/WEST/'//TRIM(aux)//'/QQ_'
      ELSEIF (device == 4) THEN
         efit_file = project_root//'/G_EQDSK/JTSA/'//TRIM(aux)//'/QQ_'
      ELSE
         WRITE(*,*) 'the specified tokamak does not exist'
      END IF

      WRITE(aux, '(I0)') shotslice

      ! WEST files use zero-padded shot-slice names.
      IF ((device == 3) .AND. (shotslice < 10)) THEN
         aux = '00'//TRIM(aux)
      ELSEIF ((device == 3) .AND. (shotslice < 100) .AND. (shotslice >= 10)) THEN
         aux = '0'//TRIM(aux)
      END IF

      efit_file = efit_file//TRIM(aux)//'.dat'

      OPEN (UNIT=100, FILE=TRIM(efit_file), STATUS='old', ACTION='read')

      ! Sign convention note:
      ! The magnetic-field direction and field-line twisting are controlled by B0
      ! and q(r). For analytical equilibria these signs must be chosen carefully in
      ! the input parameters; EFIT files already carry the intended convention.

      READ (100, *) R0
      READ (100, *) B0
      READ (100, *) NgridR
      READ (100, *) NgridZ
      READ (100, *) minR
      READ (100, *) maxR
      READ (100, *) minZ
      READ (100, *) maxZ
      READ (100, *) a0
      READ (100, *) triang
      READ (100, *) elong
      READ (100, *) Nqua

      Ngrid = NgridR*NgridZ
      ALLOCATE (Efit_data(Ngrid*Nqua))

      DO i = 1, Ngrid*Nqua
         READ (100, *) Efit_data(i)
      END DO

      CLOSE (100)

      stepR = (maxR - minR) / (NgridR - 1)
      stepZ = (maxZ - minZ) / (NgridZ - 1)

      ! Scale EFIT arrays to the dimensionless units used by the drift routines.
      minR  = minR/R0
      maxR  = maxR/R0
      minZ  = minZ/R0
      maxZ  = maxZ/R0
      stepR = stepR/R0
      stepZ = stepZ/R0

      ! Efit_data block order:
      ! {psi, psir, psiz, psirr, psirz, psizz, F, Fprim,
      !  qpsi, qprim, rhot, rhotR, rhotZ, chi, chiR, chiZ}
      Efit_data(1 +  0*Ngrid :  1*Ngrid) = Efit_data(1 +  0*Ngrid :  1*Ngrid) / (B0*R0**2)  ! psi
      Efit_data(1 +  1*Ngrid :  2*Ngrid) = Efit_data(1 +  1*Ngrid :  2*Ngrid) / (B0*R0)     ! psir
      Efit_data(1 +  2*Ngrid :  3*Ngrid) = Efit_data(1 +  2*Ngrid :  3*Ngrid) / (B0*R0)     ! psiz
      Efit_data(1 +  3*Ngrid :  4*Ngrid) = Efit_data(1 +  3*Ngrid :  4*Ngrid) / (B0)        ! psirr
      Efit_data(1 +  4*Ngrid :  5*Ngrid) = Efit_data(1 +  4*Ngrid :  5*Ngrid) / (B0)        ! psirz
      Efit_data(1 +  5*Ngrid :  6*Ngrid) = Efit_data(1 +  5*Ngrid :  6*Ngrid) / (B0)        ! psizz
      Efit_data(1 +  6*Ngrid :  7*Ngrid) = Efit_data(1 +  6*Ngrid :  7*Ngrid) / (B0*R0)     ! F
      Efit_data(1 +  7*Ngrid :  8*Ngrid) = Efit_data(1 +  7*Ngrid :  8*Ngrid) * R0          ! Fprim
      Efit_data(1 +  8*Ngrid :  9*Ngrid) = Efit_data(1 +  8*Ngrid :  9*Ngrid)               ! qpsi
      Efit_data(1 +  9*Ngrid : 10*Ngrid) = Efit_data(1 +  9*Ngrid : 10*Ngrid)               ! qprim (dq/drhot)
      Efit_data(1 + 10*Ngrid : 11*Ngrid) = Efit_data(1 + 10*Ngrid : 11*Ngrid)               ! rhot
      Efit_data(1 + 11*Ngrid : 12*Ngrid) = Efit_data(1 + 11*Ngrid : 12*Ngrid) * R0          ! rhotR
      Efit_data(1 + 12*Ngrid : 13*Ngrid) = Efit_data(1 + 12*Ngrid : 13*Ngrid) * R0          ! rhotZ
      Efit_data(1 + 13*Ngrid : 14*Ngrid) = Efit_data(1 + 13*Ngrid : 14*Ngrid)               ! chi
      Efit_data(1 + 14*Ngrid : 15*Ngrid) = Efit_data(1 + 14*Ngrid : 15*Ngrid) * R0          ! chiR
      Efit_data(1 + 15*Ngrid : 16*Ngrid) = Efit_data(1 + 15*Ngrid : 16*Ngrid) * R0          ! chiZ

   END SUBROUTINE load_efit_equilibrium

END MODULE constants
