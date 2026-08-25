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
!   (X,Y,Z) = (R,Z,varphi) are real-space coordinates, the cylindrical coordinates
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
   REAL(KIND=wp), PARAMETER :: pi   = 3.1415926545_wp         ! pi
   REAL(KIND=wp), PARAMETER :: E    = 2.7182818284_wp         ! Euler's number
   REAL(KIND=wp), PARAMETER :: q0   = 1.60217*10.0**(-19.0)   ! elementary charge
   REAL(KIND=wp), PARAMETER :: mi   = 1.67262*10.0**(-27.0)   ! proton mass
   REAL(KIND=wp), PARAMETER :: me   = 9.10938*10.0**(-31.0)   ! electron mass
   REAL(KIND=wp), PARAMETER :: eps  = 8.85418*10.0**(-12.0)   ! vacuum dielectric constant
   REAL(KIND=wp), PARAMETER :: mu0  = 1.25663*10.0**(-6.0)    ! vacuum magnetic permeability
   REAL(KIND=wp), PARAMETER :: kB   = 1.38064*10.0**(-23.0)   ! Boltzmann constant
   REAL(KIND=wp), PARAMETER :: logL = 17.0                    ! Coulomb logarithm
   INTEGER, PARAMETER       :: ON = 1
   INTEGER, PARAMETER       :: OFF = 0
   REAL(KIND=wp), PARAMETER :: eps_xi    = 1.0e-15_wp
   REAL(KIND=wp), PARAMETER :: eps_rr    = 1.0e-15_wp
   REAL(KIND=wp), PARAMETER :: eps_root  = 1.0e-15_wp
   REAL(KIND=wp), PARAMETER :: eps_B     = 1.0e-15_wp
   REAL(KIND=wp), PARAMETER :: eps_Bsp   = 1.0e-15_wp
   REAL(KIND=wp), PARAMETER :: eps_omega = 1.0e-15_wp
   REAL(KIND=wp), PARAMETER :: eps_large = 1.0e-6_wp


   !---------------------------------------------------------------------------------
   ! Switches
   !---------------------------------------------------------------------------------
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
   REAL(KIND=wp) :: t0                             ! simulation start time        [R0/vth]
   REAL(KIND=wp) :: tc                             ! collisions start time        [R0/vth]
   REAL(KIND=wp) :: tt                             ! turbulence start time        [R0/vth]
   REAL(KIND=wp) :: tmax                           ! simulation end time          [R0/vth]

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
   REAL(KIND=wp) :: Ti                             ! ion temperature @ r00
   REAL(KIND=wp) :: Te                             ! electron temperature @ r00
   REAL(KIND=wp) :: Zeff                           ! effective charge number
   REAL(KIND=wp) :: Aeff                           ! effective mass number
   REAL(KIND=wp) :: ndens                          ! ion density @ r00             [m^-3] (imported in 1e19 m^-3)

   REAL(KIND=wp) :: vth                            ! thermal velocity (H ions)     [m/s] @ r00
   REAL(KIND=wp) :: rhoi                           ! Larmor radius (H ions)        [m]   @ r00 (B0)
   REAL(KIND=wp) :: wi                             ! Larmor frequency              [Hz]
   REAL(KIND=wp) :: c0                             ! collision-related constant
   REAL(KIND=wp) :: delta                          ! collision-related constant

   !---------------------------------------------------------------------------------
   ! Gradients / drives
   !---------------------------------------------------------------------------------
   REAL(KIND=wp) :: Ln                             ! density gradient length       [R0*d(lnn)/dr]
   REAL(KIND=wp) :: Li                             ! ion temperature gradient      [R0*d(lnTi)/dr]
   REAL(KIND=wp) :: Le                             ! electron temperature gradient [R0*d(lnTe)/dr]

   !---------------------------------------------------------------------------------
   ! Geometry / equilibrium
   !---------------------------------------------------------------------------------
   INTEGER       :: magnetic_model                 ! equilibrium: 1-EFIT, 2-Solovev, 3-Circular

   REAL(KIND=wp) :: B0                             ! magnetic field @ axis         [T]
   REAL(KIND=wp) :: R0                             ! major radius                  [m]
   REAL(KIND=wp) :: a0                             ! minor radius                  [m]
   REAL(KIND=wp) :: safety                         ! reference analytical safety factor
   REAL(KIND=wp) :: shear_mag                      ! analytical magnetic shear

   REAL(KIND=wp) :: s1                             ! circular q-profile coefficient
   REAL(KIND=wp) :: s2                             ! circular q-profile coefficient
   REAL(KIND=wp) :: s3                             ! circular q-profile coefficient

   REAL(KIND=wp) :: amp                            ! Solovev flux amplitude
   REAL(KIND=wp) :: elong                          ! elongation

   INTEGER       :: device                         ! (1=MAST, 2=TCV, 3=WEST, 4=JTSA)
   INTEGER       :: shot                           ! shot number
   INTEGER       :: shotslice                      ! time (ms) for equilibrium slice

   REAL(KIND=wp) :: Omgt0                          ! toroidal angular frequency    [1e3 Hz] @ r00
   REAL(KIND=wp) :: Omgtprim                       ! d ln Omgt / d rhot            [1] @ r00

   REAL(KIND=wp) :: q00                            ! safety factor @ r00
   REAL(KIND=wp) :: r00                            ! reference surface (for global evaluation)
   REAL(KIND=wp) :: rhot0                          ! toroidal normalized radius @ r00
   REAL(KIND=wp) :: psi0                           ! poloidal flux @ r00

   !---------------------------------------------------------------------------------
   ! Solovev precomputed straight-field-line geometry
   !---------------------------------------------------------------------------------
   INTEGER, PARAMETER :: Nsol     = 257     ! radial table points; must be odd
   INTEGER, PARAMETER :: Ksol     = 16      ! retained cosine harmonics
   INTEGER, PARAMETER :: NetaSol  = 1024    ! periodic eta quadrature points
   INTEGER, PARAMETER :: NedgeSol = 1024    ! quadrature points for Phi_t at separatrix

   ! The straight-field-line Fourier representation is deliberately stopped
   ! inside the separatrix, where q and chi_eta become singular.
   ! rpsi itself remains the physical coordinate sqrt(psi/psi_sep).
   REAL(KIND=wp), PARAMETER :: SOL_RMAX = 0.95_wp

   ! Flat Sol_data layout.
   INTEGER, PARAMETER :: SOL_W0_OFF    = 0
   INTEGER, PARAMETER :: SOL_W0R_OFF   = Nsol
   INTEGER, PARAMETER :: SOL_T_OFF     = 2*Nsol
   INTEGER, PARAMETER :: SOL_B_OFF     = 3*Nsol                 ! b_k blocks
   INTEGER, PARAMETER :: SOL_BR_OFF    = (3 + Ksol)*Nsol        ! db_k/drpsi blocks
   INTEGER, PARAMETER :: SOL_JEDGE_IDX = (3 + 2*Ksol)*Nsol + 1  ! true separatrix Jedge
   INTEGER, PARAMETER :: SOL_NDATA     = SOL_JEDGE_IDX

   ! Uniform table spacing in the physical coordinate rpsi=sqrt(psi/psi_sep).
   REAL(KIND=wp), PARAMETER :: SOL_DR     = SOL_RMAX/REAL(Nsol - 1, wp)
   REAL(KIND=wp), PARAMETER :: SOL_INV_DR = REAL(Nsol - 1, wp)/SOL_RMAX

   REAL(KIND=wp), ALLOCATABLE :: Sol_data(:)

   ! Extra equilibrium / imported-data helpers
   REAL(KIND=wp) :: alfa                           ! Solovev parameter (separatrix/ellipticity-related)
   REAL(KIND=wp) :: gama                           ! Solovev parameter (sqrt(gama)=separatrix position)
   CHARACTER(LEN=:), ALLOCATABLE            :: efit_file
   REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: Efit_data
   INTEGER                                  :: NgridR, NgridZ, Nqua, Ngrid
   REAL(KIND=wp)                            :: minR, maxR, minZ, maxZ
   REAL(KIND=wp)                            :: stepR, stepZ
   REAL(KIND=wp)                            :: triang

   !---------------------------------------------------------------------------------
   ! Turbulence model
   !---------------------------------------------------------------------------------
   INTEGER       :: turb_model                     ! type of turbulence model: 1=old, 2=new
   INTEGER       :: x_corr                         ! x correlation: 1-gaussian, 2-exponential, 3-GENE-like
   INTEGER       :: y_corr                         ! y correlation: 1-Madi, 2-GENE
   INTEGER       :: t_corr                         ! t correlation: 1-gaussian, 2-exponential

   REAL(KIND=wp) :: Phi                            ! turbulence amplitude [ePhi/Ti]
   REAL(KIND=wp) :: Beta                           ! magnetic turbulence amplitude [e vth A_parallel/Ti]
   REAL(KIND=wp) :: turbprof                       ! exponent for synthetic turb. profile
   REAL(KIND=wp) :: Ai                             ! ITG fraction

   REAL(KIND=wp) :: lambdax                        ! x correlation length [rhoi units]
   REAL(KIND=wp) :: lambday                        ! y correlation length [rhoi units]
   REAL(KIND=wp) :: lambdaz                        ! z correlation length due to decoherence [q0*R0 units]
   REAL(KIND=wp) :: lbalonz                        ! z correlation length due to ballooning envelope [q0*R0 units]
   REAL(KIND=wp) :: tauc                           ! time correlation [R0/vth]
   REAL(KIND=wp) :: k0i                            ! ITG dominant ky [1/rhoi]
   REAL(KIND=wp) :: k0e                            ! TEM dominant ky [1/rhoi]
   REAL(KIND=wp) :: gamma_ZF                       ! zonal-flow shearing rate
   REAL(KIND=wp) :: gamma_E                        ! rotational shearing rate

   REAL(KIND=wp) :: Ae                             ! TEM fraction
   INTEGER       :: dmmax                          ! maximal parallel number

   !---------------------------------------------------------------------------------
   ! Particle initial conditions / markers
   !---------------------------------------------------------------------------------
   REAL(KIND=wp) :: X0                             ! starting position R
   REAL(KIND=wp) :: Y0                             ! starting position Z
   REAL(KIND=wp) :: Z0                             ! starting position phi

   REAL(KIND=wp) :: Ts                             ! energy average [Ti units] @ r00
   REAL(KIND=wp) :: Es                             ! average ion energy [Ti units]
   REAL(KIND=wp) :: pitch                          ! average pitch angle
   REAL(KIND=wp) :: As                             ! ion mass number
   REAL(KIND=wp) :: Zs                             ! ionization state
   REAL(KIND=wp) :: taucc                          ! inverse collision frequency [R0/vth]
   REAL(KIND=wp) :: Lns                            ! density gradient for the species
   REAL(KIND=wp) :: Lts                            ! temperature gradient of the species
   REAL(KIND=wp) :: annulus_width                  ! full radial width for local-annulus initialization

   INTEGER       :: position_type                  ! init position: 1-fixed point, ...
   INTEGER       :: pitch_type                     ! init pitch: 1-fixed, 2-random uniform
   INTEGER       :: energy_type                    ! init energy: 1-fixed, 2-Boltzmann

   REAL(KIND=wp) :: usetilt, balloon, noballoon, norm ! helpers

   !---------------------------------------------------------------------------------
   ! Misc. derived helpers
   !---------------------------------------------------------------------------------
   REAL(KIND=wp) :: C1                             ! q1 = C1*f(psi)
   REAL(KIND=wp) :: C2                             ! q2 = C2*(zeta - q(psi)*chi)
   REAL(KIND=wp) :: C3                             ! q3 = C3*chi

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
      REAL(KIND=wp), DIMENSION(cols1) :: pp       ! parameter vector for this run
      INTEGER                         :: ierr_sol

      pp = array1(run, :)

      !---------------------------------------------------------------------------------
      ! Import raw parameters from the selected run row.
      !
      ! Name-based indexing keeps this code robust if the GUI reorders columns.
      ! Integer parameters are converted explicitly; physical quantities remain real.
      !---------------------------------------------------------------------------------
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
      safety         = pp(parameter_index("safety"))
      shear_mag      = pp(parameter_index("shear_mag"))
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
      Beta           = pp(parameter_index("Beta"))
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
      
      Nqua = 16   ! number of quantities stored on the EFIT-like grid
      IF (magnetic_model == 1) CALL load_efit_equilibrium()
      CALL validate_raw_run_parameters

      !---------------------------------------------------------------------------------
      ! Convert imported units and evaluate derived normalization quantities.
      !---------------------------------------------------------------------------------
      Ti    = Ti*1000.0_wp                     ! imported in keV
      Te    = Te*1000.0_wp                     ! imported in keV
      ndens = ndens*10.0_wp**19                ! imported in 1e19 m^-3

      vth  = sqrt(Ti*q0/mi)
      rhoi = mi*vth/q0/abs(B0)
      wi   = q0*abs(B0)/mi

      r00 = sqrt((X0 - R0)**2 + Y0**2)

      IF (.NOT. ALLOCATED(Sol_data)) ALLOCATE(Sol_data(SOL_NDATA))
      IF (magnetic_model == 3) THEN
           s1 = safety*(1.0_wp - shear_mag) + eps_large
           s2 = safety*shear_mag*a0/r00 + eps_large
           s3 = eps_large
      ENDIF

      IF (magnetic_model == 2) THEN
	epsa = a0/R0
	gama = 1.0_wp - 2.0_wp*epsa*sqrt(1.0_wp-epsa**2)
	alfa = elong*(sqrt(2.0_wp-gama)-sqrt(gama))/sqrt(2.0_wp*(1.0_wp-gama))

	xref = X0/R0
	zref = Y0/R0

	dsep = 1.0_wp-gama
	rpsi_ref = sqrt((xref*xref-1.0_wp)**2 + 4.0_wp*zref*zref*(xref*xref-gama)/(alfa*alfa))/dsep
	w0_ref = solovev_w0_at_rpsi(rpsi_ref,gama)
	disc = safety*safety - 0.25_wp*gama*dsep*dsep*rpsi_ref*rpsi_ref*w0_ref*w0_ref

	IF (disc <= 0.0_wp) THEN
	   ERROR STOP 'Requested Solovev q is incompatible with geometry'
	END IF

	amp = SIGN(1.0_wp,safety)*(1.0_wp+alfa*alfa)*w0_ref/alfa/SQRT(disc)

	CALL build_solovev_data(Sol_data,ierr_sol)

         IF (ierr_sol /= 0) THEN
            WRITE(*, *) 'ERROR building Solovev table. ierr = ', ierr_sol
            ERROR STOP 'Failed to build Solovev geometry data'
         END IF
      ELSE
         Sol_data = 0.0_wp 
      END IF

      a0  = a0/R0
      r00 = r00/R0
      X0  = X0/R0
      Y0  = Y0/R0
      Omgt0 = Omgt0*10.0_wp**3 / (vth/R0)
                 
      CALL equilibrium_q_rhot_psi(X0, Y0, q00, rhot0, psi0)
      C1 = a0*R0/rhoi
      C2 = -r00/q00 * R0/rhoi
      C3 = 1.0_wp

      Nci = int(Nc*Ai**2)
      Nce = Nc - Nci

      Ae = sqrt(1.0_wp - Ai**2)

      c0 = logL*ndens/8.0_wp/pi * (Zeff*Zs*q0**2/(eps*mi*As))**2

      delta = 5.0_wp*4.0_wp/3.0_wp/sqrt(pi) * (tmax - t0)/Nt * R0*2.0_wp*logL*ndens*(Zeff*Zs*q0**2)**2 / &
              (pi*eps**2*mi**2*As**2*vth**4*sqrt(2.0_wp/Aeff)**3)
      taucc = 0.0_wp

      CALL validate_derived_run_parameters

      gamma_E = gamma_ZF + Omgtprim          ! total shearing rate used by turbulence
      usetilt = real(USE_tilt, wp)
      balloon = real(USE_balloon, wp)        ! integer switch -> 0.0 or 1.0
      noballoon = 1.0_wp - balloon
      norm = 1.0_wp / sqrt(real(Nc, wp)/2.0_wp)

      ! dmmax = 2
      dmmax = 2 + int(sqrt(-2.0_wp*log(0.05_wp)/lbalonz))

   END SUBROUTINE load_run_parameters

   SUBROUTINE load_efit_equilibrium()
      USE sims, ONLY: project_root
      IMPLICIT NONE

      CHARACTER(LEN=35) :: aux
      INTEGER           :: i

      IF (magnetic_model /= 1) RETURN

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
   
     
   
   SUBROUTINE build_solovev_data(data_out, ierr)
      IMPLICIT NONE

      !-------------------------------------------------------------------
      ! Outputs
      !-------------------------------------------------------------------
      REAL(KIND=wp), INTENT(OUT) :: data_out(SOL_NDATA)
      INTEGER, INTENT(OUT)       :: ierr

      !-------------------------------------------------------------------
      ! Local arrays
      !-------------------------------------------------------------------
      REAL(KIND=wp) :: qnode(Nsol)
      REAL(KIND=wp) :: fJ(Nsol)
      REAL(KIND=wp) :: Jnode(Nsol)
      REAL(KIND=wp) :: asum(Ksol)
      REAL(KIND=wp) :: arsum(Ksol)
      REAL(KIND=wp) :: coseta(NetaSol)
      REAL(KIND=wp) :: coseta_mid(NetaSol)

      !-------------------------------------------------------------------
      ! Scalars
      !-------------------------------------------------------------------
      REAL(KIND=wp) :: delta_sep, delta
      REAL(KIND=wp) :: psi_sep, psi
      REAL(KIND=wp) :: rpsi
      REAL(KIND=wp) :: eta, deta
      REAL(KIND=wp) :: ce
      REAL(KIND=wp) :: d1, d2
      REAL(KIND=wp) :: w, wr
      REAL(KIND=wp) :: w0, w0r
      REAL(KIND=wp) :: ak, akr
      REAL(KIND=wp) :: bk, bkr
      REAL(KIND=wp) :: sum0, sum0r
      REAL(KIND=wp) :: ck, ckm1, ckp1
      REAL(KIND=wp) :: Fpsi, Farg
      REAL(KIND=wp) :: Cqsol
      REAL(KIND=wp) :: dr
      REAL(KIND=wp) :: Jedge
      REAL(KIND=wp) :: tedge, redge, dtedge
      REAL(KIND=wp) :: qedge_int

      INTEGER :: i, j, k
      INTEGER :: ib, ibr

      !===================================================================
      ! Initialization / sanity checks
      !===================================================================
      ierr = 0
      data_out = 0.0_wp

      IF (alfa <= 0.0_wp) THEN
         ierr = 1
         RETURN
      END IF

      IF (MOD(Nsol - 1, 2) /= 0) THEN
         ! Cumulative Simpson construction below assumes an even number
         ! of radial intervals.
         ierr = 2
         RETURN
      END IF

      IF (NetaSol <= 2*Ksol) THEN
         ierr = 3
         RETURN
      END IF

      IF ((SOL_RMAX <= 0.0_wp) .OR. (SOL_RMAX >= 1.0_wp)) THEN
         ierr = 4
         RETURN
      END IF

      !===================================================================
      ! Solovev separatrix normalization
      !
      ! For this parametrization the separatrix is delta_sep = 1-gamma.
      ! Its midplane intersections are R^2=gamma and R^2=2-gamma.
      ! Therefore
      !
      !   rpsi = sqrt(psi/psi_sep) = delta/delta_sep,
      !
      ! and the straight-angle table is built only up to SOL_RMAX<1.
      !===================================================================
      delta_sep = 1.0_wp - gama

      IF ((delta_sep <= 0.0_wp) .OR. (gama <= 0.0_wp)) THEN
         ierr = 5
         RETURN
      END IF

      psi_sep = amp*alfa**2*delta_sep**2 &
              /(8.0_wp*(1.0_wp + alfa**2))

      deta  = 2.0_wp*Pi/REAL(NetaSol, wp)
      dr    = SOL_DR
      Cqsol = (1.0_wp + alfa**2)/(alfa*amp)

      ! Periodic nodes for the Fourier table, and half-shifted nodes for
      ! the separate separatrix-flux quadrature.
      DO j = 1, NetaSol
         eta = deta*REAL(j - 1, wp)
         coseta(j) = COS(eta)

         eta = deta*(REAL(j - 1, wp) + 0.5_wp)
         coseta_mid(j) = COS(eta)
      END DO

      !===================================================================
      ! Radial table: 0 <= rpsi <= SOL_RMAX
      !===================================================================
      DO i = 1, Nsol

         rpsi = REAL(i - 1, wp)*dr
         delta = delta_sep*rpsi
         psi = psi_sep*rpsi*rpsi

         !---------------------------------------------------------------
         ! F(psi)
         !---------------------------------------------------------------
         Farg = 1.0_wp &
              + 2.0_wp*amp*gama*psi/(1.0_wp + alfa**2)

         IF (Farg <= 0.0_wp) THEN
            ierr = 6
            RETURN
         END IF

         Fpsi = SQRT(Farg)

         !---------------------------------------------------------------
         ! Fourier integrals of
         !
         !   w = 1 / [(1+delta cos eta)
         !            sqrt(1-gamma+delta cos eta)]
         !
         ! and its exact derivative with respect to physical rpsi:
         !
         !   dw/drpsi = -delta_sep cos(eta) w *
         !      [1/(1+delta cos eta)
         !       + 1/(2(1-gamma+delta cos eta))].
         !---------------------------------------------------------------
         sum0  = 0.0_wp
         sum0r = 0.0_wp
         asum  = 0.0_wp
         arsum = 0.0_wp

         DO j = 1, NetaSol

            ce = coseta(j)
            d1 = 1.0_wp - 0.0_wp + delta*ce
            d2 = delta_sep + delta*ce

            IF ((d1 <= 0.0_wp) .OR. (d2 <= 0.0_wp)) THEN
               ierr = 7
               RETURN
            END IF

            w = 1.0_wp/(d1*SQRT(d2))

            wr = -delta_sep*ce*w &
               *(1.0_wp/d1 + 0.5_wp/d2)

            sum0  = sum0  + w
            sum0r = sum0r + wr

            ! cos(k eta) recurrence: no Ksol transcendental calls.
            ckm1 = 1.0_wp
            ck   = ce

            DO k = 1, Ksol
               asum(k)  = asum(k)  + w*ck
               arsum(k) = arsum(k) + wr*ck

               IF (k < Ksol) THEN
                  ckp1 = 2.0_wp*ce*ck - ckm1
                  ckm1 = ck
                  ck   = ckp1
               END IF
            END DO

         END DO

         !---------------------------------------------------------------
         ! w = w0 + Sum_k a_k cos(k eta)
         !---------------------------------------------------------------
         w0  = sum0 /REAL(NetaSol, wp)
         w0r = sum0r/REAL(NetaSol, wp)

         data_out(SOL_W0_OFF  + i) = w0
         data_out(SOL_W0R_OFF + i) = w0r

         !---------------------------------------------------------------
         ! Normalized straight-angle Fourier coefficients
         !
         !   b_k = a_k/w0,
         !   b'_k = (a'_k w0 - a_k w'_0)/w0^2.
         !---------------------------------------------------------------
         DO k = 1, Ksol

            ak  = 2.0_wp*asum(k) /REAL(NetaSol, wp)
            akr = 2.0_wp*arsum(k)/REAL(NetaSol, wp)

            bk  = ak/w0
            bkr = (akr*w0 - ak*w0r)/(w0*w0)

            ib  = SOL_B_OFF  + (k - 1)*Nsol + i
            ibr = SOL_BR_OFF + (k - 1)*Nsol + i

            data_out(ib)  = bk
            data_out(ibr) = bkr

         END DO

         ! Exact symmetry values on the magnetic axis.
         IF (i == 1) THEN
            data_out(SOL_W0R_OFF + i) = 0.0_wp

            DO k = 1, Ksol
               ib = SOL_B_OFF + (k - 1)*Nsol + i
               data_out(ib) = 0.0_wp

               ! b_1'(0) is generally non-zero; higher first derivatives vanish.
               IF (k >= 2) THEN
                  ibr = SOL_BR_OFF + (k - 1)*Nsol + i
                  data_out(ibr) = 0.0_wp
               END IF
            END DO
         END IF

         !---------------------------------------------------------------
         ! q(rpsi) and toroidal-flux radial integrand
         !---------------------------------------------------------------
         qnode(i) = Cqsol*Fpsi*w0
         fJ(i) = 2.0_wp*rpsi*qnode(i)

      END DO

      !===================================================================
      ! Cumulative J(rpsi)=Integral_0^r 2 r' q(r') dr' on the table domain.
      ! Simpson integration in pairs gives J at every radial table node.
      !===================================================================
      Jnode = 0.0_wp
      Jnode(1) = 0.0_wp

      DO i = 1, Nsol - 2, 2

         Jnode(i + 1) = Jnode(i) &
            + dr*(5.0_wp*fJ(i) + 8.0_wp*fJ(i + 1) - fJ(i + 2)) &
            /12.0_wp

         Jnode(i + 2) = Jnode(i) &
            + dr*(fJ(i) + 4.0_wp*fJ(i + 1) + fJ(i + 2)) &
            /3.0_wp

      END DO

      !===================================================================
      ! True separatrix toroidal-flux normalization Jedge.
      !
      ! q(rpsi) diverges logarithmically as rpsi->1, but Jedge is finite.
      ! We therefore use the endpoint-removing substitution
      !
      !      rpsi = 1 - t^2,    0<t<1,
      !
      ! so that
      !
      ! Jedge = Integral_0^1 4 t rpsi q(rpsi) dt.
      !
      ! Midpoint nodes avoid evaluating q exactly at the separatrix.  The
      ! eta quadrature is half-shifted as well, which resolves the X-point
      ! peak much better than sampling eta=pi directly.
      !===================================================================
      Jedge = 0.0_wp
      dtedge = 1.0_wp/REAL(NedgeSol, wp)

      DO i = 1, NedgeSol

         tedge = (REAL(i, wp) - 0.5_wp)*dtedge
         redge = 1.0_wp - tedge*tedge
         delta = delta_sep*redge
         psi   = psi_sep*redge*redge

         Farg = 1.0_wp &
              + 2.0_wp*amp*gama*psi/(1.0_wp + alfa**2)

         IF (Farg <= 0.0_wp) THEN
            ierr = 8
            RETURN
         END IF

         Fpsi = SQRT(Farg)
         sum0 = 0.0_wp

         DO j = 1, NetaSol
            ce = coseta_mid(j)
            d1 = 1.0_wp + delta*ce
            d2 = delta_sep + delta*ce

            IF ((d1 <= 0.0_wp) .OR. (d2 <= 0.0_wp)) THEN
               ierr = 9
               RETURN
            END IF

            sum0 = sum0 + 1.0_wp/(d1*SQRT(d2))
         END DO

         w0 = sum0/REAL(NetaSol, wp)
         qedge_int = Cqsol*Fpsi*w0

         Jedge = Jedge &
               + dtedge*4.0_wp*tedge*redge*qedge_int

      END DO

      IF (ABS(Jedge) <= TINY(1.0_wp)) THEN
         ierr = 10
         RETURN
      END IF

      !===================================================================
      ! T(rpsi)=rho_t^2 = J(rpsi)/Jedge.
      !
      ! Important: the last table point is at SOL_RMAX<1, therefore
      ! T(Nsol)<1.  rho_t=1 belongs to the true separatrix, not the table edge.
      !===================================================================
      DO i = 1, Nsol
         data_out(SOL_T_OFF + i) = Jnode(i)/Jedge
      END DO

      data_out(SOL_T_OFF + 1) = 0.0_wp
      data_out(SOL_JEDGE_IDX) = Jedge

   END SUBROUTINE build_solovev_data

PURE FUNCTION solovev_w0_at_rpsi(rpsi, gama) RESULT(w0)
   IMPLICIT NONE

   REAL(KIND=wp), INTENT(IN) :: rpsi, gama
   REAL(KIND=wp) :: w0
   REAL(KIND=wp) :: d, delta_loc, eta, ce
   REAL(KIND=wp) :: d1, d2, deta
   INTEGER :: j

   d         = 1.0_wp - gama
   delta_loc = d*rpsi
   deta      = 2.0_wp*pi/REAL(NetaSol,wp)

   w0 = 0.0_wp

   DO j = 1, NetaSol
      eta = deta*REAL(j-1,wp)
      ce  = COS(eta)

      d1 = 1.0_wp + delta_loc*ce
      d2 = d        + delta_loc*ce

      w0 = w0 + 1.0_wp/(d1*SQRT(d2))
   END DO

   w0 = w0/REAL(NetaSol,wp)

END FUNCTION solovev_w0_at_rpsi

END MODULE constants
