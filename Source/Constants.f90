!========================================================================================================
! Module: constants
! Shared constants / parameters across the project.
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
!========================================================================================================

MODULE constants
   USE sims
   IMPLICIT NONE
   SAVE   ! all data are preserved in different procedures

!========================================================================================================
! PART 1: MAIN MODULE BODY  variable declarations
!========================================================================================================

   !---------------------------------------------------------------------------------
   ! Universal constants
   !---------------------------------------------------------------------------------
   REAL(KIND=dp), PARAMETER :: pi   = 3.1415926545_dp         ! pi
   REAL(KIND=dp), PARAMETER :: E    = 2.7182818284_dp         ! Euler's number
   REAL(KIND=dp), PARAMETER :: q0   = 1.60217*10.0**(-19.0)   ! elementary charge
   REAL(KIND=dp), PARAMETER :: mi   = 1.67262*10.0**(-27.0)   ! proton mass
   REAL(KIND=dp), PARAMETER :: me   = 9.10938*10.0**(-31.0)   ! electron mass
   REAL(KIND=dp), PARAMETER :: eps  = 8.85418*10.0**(-12.0)   ! vacuum dielectric constant
   REAL(KIND=dp), PARAMETER :: mu0  = 1.25663*10.0**(-6.0)    ! vacuum magnetic permeability
   REAL(KIND=dp), PARAMETER :: kB   = 1.38064*10.0**(-23.0)   ! Boltzmann constant
   REAL(KIND=dp), PARAMETER :: logL = 17.0                    ! Coulomb logarithm
   INTEGER, PARAMETER       :: ON = 1 
   INTEGER, PARAMETER       :: OFF = 0 
   REAL(KIND=dp), parameter :: eps_xi   = 1.0e-15_dp
   REAL(KIND=dp), parameter :: eps_rr   = 1.0e-15_dp
   REAL(KIND=dp), parameter :: eps_root = 1.0e-15_dp
   REAL(KIND=dp), parameter :: eps_B    = 1.0e-15_dp
   REAL(KIND=dp), parameter :: eps_Bsp  = 1.0e-15_dp   
   REAL(KIND=dp), parameter :: eps_omega= 1.0e-15_dp   

   !---------------------------------------------------------------------------------
   ! Numerical parameters
   !---------------------------------------------------------------------------------
   REAL(KIND=dp) :: t0                             ! simulation start time        [R0/vth]
   REAL(KIND=dp) :: tc                             ! collisions start time        [R0/vth]
   REAL(KIND=dp) :: tt                             ! turbulence start time        [R0/vth]
   REAL(KIND=dp) :: tmax                           ! simulation end time          [R0/vth]

   INTEGER       :: Np                             ! particles per loop (total: Nloop*Np)
   INTEGER       :: Nloop                          ! simulation loops (super-ensemble)
   INTEGER       :: Nreal                          ! realizations of turbulent field
   INTEGER       :: Nc                             ! number of waves (typically ~1e3)
   INTEGER       :: Nt                             ! number of time steps
   INTEGER       :: ntraj                          ! number of trajectories exported
   INTEGER       :: Nci                            ! number of ITG partial-waves
   INTEGER       :: Nce                            ! number of TEM partial-waves

   !---------------------------------------------------------------------------------
   ! Equilibrium
   !---------------------------------------------------------------------------------
   REAL(KIND=dp) :: Ti                             ! ion temperature @ r00         [eV]
   REAL(KIND=dp) :: Te                             ! electron temperature @ r00    [eV]
   REAL(KIND=dp) :: B0                             ! magnetic field @ axis         [T]
   REAL(KIND=dp) :: vth                            ! thermal velocity (H ions)     [m/s] @ r00
   REAL(KIND=dp) :: rhoi                           ! Larmor radius (H ions)        [m]   @ r00 (B0)
   REAL(KIND=dp) :: R0                             ! major radius                  [m]
   REAL(KIND=dp) :: Ln                             ! ITG grad-length               [R0*d(lnn)/dr]
   REAL(KIND=dp) :: Li                             ! ITG grad-length               [R0*d(lnTi)/dr]
   REAL(KIND=dp) :: Le                             ! TEM grad-length               [R0*d(lnTe)/dr]
   REAL(KIND=dp) :: Aeff                           ! effective mass number
   REAL(KIND=dp) :: Zeff                           ! effective charge number
   REAL(KIND=dp) :: ndens                          ! ion density @ r00             [m^-3] (imported in 1e19 m^-3)
   REAL(KIND=dp) :: a0                             ! minor radius                  [m]
   REAL(KIND=dp) :: Omgt0                          ! toroidal angular frequency    [1e3 Hz] @ r00
   REAL(KIND=dp) :: Omgtprim                       ! d ln Omgt / d rhot            [adimensional]     @ r00
   REAL(KIND=dp) :: wi                             ! Larmor frequency              [Hz]
   REAL(KIND=dp) :: q00                            ! safety factor @ r00           [1]

   !---------------------------------------------------------------------------------
   ! Circular model parameters (indexmodel = 4): q = safe1 + safe2*(r/a) + safe3*(r/a)^2
   !---------------------------------------------------------------------------------
   REAL(KIND=dp) :: safe1
   REAL(KIND=dp) :: safe2
   REAL(KIND=dp) :: safe3

   !---------------------------------------------------------------------------------
   ! Solovev model parameters (indexmodel = 3)
   !---------------------------------------------------------------------------------
   REAL(KIND=dp) :: amp                            ! Solovev flux amplitude
   REAL(KIND=dp) :: alfa                           ! Solovev parameter (separatrix/ellipticity-related)
   REAL(KIND=dp) :: elong                          ! elongation
   REAL(KIND=dp) :: gama                           ! Solovev parameter (sqrt(gama)=separatrix position)

   !---------------------------------------------------------------------------------
   ! EFIT input (indexmodel = 1,2)
   !---------------------------------------------------------------------------------
   INTEGER                                  :: device       ! (1=MAST, 2=TCV, 3=WEST, 4=JTSA)
   INTEGER                                  :: shot         ! shot number
   INTEGER                                  :: shotslice    ! time (ms) for equilibrium slice
   CHARACTER(LEN=:), ALLOCATABLE            :: efit_file    ! EFIT file name
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Efit_data    ! equilibrium grid values

   INTEGER                                  :: NgridR, NgridZ, Nqua, Ngrid     ! EFIT grid sizes / number of quantities
   REAL(KIND=dp)                            :: minR, maxR, minZ, maxZ          ! grid extents
   REAL(KIND=dp)                            :: stepR, stepZ                    ! grid steps (computed)
   REAL(KIND=dp)                            :: triang                          ! triangularity

   !---------------------------------------------------------------------------------
   ! Turbulence parameters
   !---------------------------------------------------------------------------------
   REAL(KIND=dp) :: Phi                             ! turbulence amplitude [ePhi/Ti]
   REAL(KIND=dp) :: turbprof                        ! exponent for synthetic turb. profile
   REAL(KIND=dp) :: Ai                              ! ITG fraction
   REAL(KIND=dp) :: Ae                              ! TEM fraction
   REAL(KIND=dp) :: lambdax                         ! x correlation length [rhoi units]
   REAL(KIND=dp) :: lambday                         ! y correlation length [rhoi units]
   REAL(KIND=dp) :: lambdaz                         ! z correlation length due to decoherence [q0*R0 units]
   REAL(KIND=dp) :: lbalonz                         ! z correlation length due to ballooning envelope [q0*R0 units]
! both lambdaz, and lbalonz are basically "angle" correlation lengths since they are given in qR units;
   REAL(KIND=dp) :: tauc                            ! time correlation [R0/vth]
   REAL(KIND=dp) :: k0i                             ! ITG dominant ky [1/rhoi]
   REAL(KIND=dp) :: k0e                             ! TEM dominant ky [1/rhoi]
   REAL(KIND=dp) :: gamma_E                         ! rotational sheering [1/t0]
   REAL(KIND=dp) :: gamma_ZF                        ! Zonal flow rotational sheering [1/t0]

   !---------------------------------------------------------------------------------
   ! Particle parameters
   !---------------------------------------------------------------------------------
   REAL(KIND=dp) :: X0, Y0, Z0                      ! starting positions
   REAL(KIND=dp) :: r00                             ! reference surface (for global evaluation)
   REAL(KIND=dp) :: q10                             ! reference surface (for global evaluation)
   REAL(KIND=dp) :: Tw                              ! energy average [Ti units] @ r00
   REAL(KIND=dp) :: Ew                              ! average ion energy [Ti units]
   REAL(KIND=dp) :: pitch                           ! average pitch angle
   REAL(KIND=dp) :: Aw                              ! ion mass number
   REAL(KIND=dp) :: Zw                              ! ionization state
   REAL(KIND=dp) :: c0                              ! collision-related constant
   REAL(KIND=dp) :: delta                           ! collision-related constant
   REAL(KIND=dp) :: taucc                           ! inverse collision frequency (simplified model)
   REAL(KIND=dp) :: psi0                            ! poloidal flux @ r0

   !---------------------------------------------------------------------------------
   ! Choices / switches
   !---------------------------------------------------------------------------------
   INTEGER :: magnetic_model                        ! equilibrium: 1-EFIT, 2-EFIT-sa, 3-Solovev, 4-Circular
   INTEGER :: USE_larmor                            ! FLR effects: 0 = OFF, 1 = ON, else = ERROR
   INTEGER :: position_type                         ! init position: 1-all same (X0,Y0,Z0)
   INTEGER :: pitch_type                            ! init pitch: 1-fixed, 2-random uniform (2pi)
   INTEGER :: energy_type                           ! init energy: 1-fixed (Ew), 2-Boltzmann
   INTEGER :: USE_coll                              ! collisions:   0 = OFF, 1 = ON, else = ERROR
   INTEGER :: USE_turb                              ! turbulence:   0 = OFF, 1 = ON, else = ERROR
   INTEGER :: USE_magnturb                          ! RMP:          0 = OFF, 1 = ON, else = ERROR
   INTEGER :: USE_freq                              ! mode freqs:   0 = frozen turbulence, 1 = real frequencies, else = ERROR
   INTEGER :: USE_polar                             ! polar. drift: 0 = OFF, 1 = ON, else = ERROR
   INTEGER :: USE_PC                                ! initial condition constancy psi/Pc:  0 = NO, 1 = YES, else = ERROR
   INTEGER :: USE_real                              ! turb realizations differ/same across particles: 0 = NO, 1 = YES, else = ERROR
   INTEGER :: USE_corr                              ! compute & export Lagrangian correlations: 0 = NO, 1 = YES, else = ERROR
   INTEGER :: USE_balloon                           ! 0 = homogeneous parallel structure (kz); 1 = ballooning parallel structure g(z); else = ERROR
   INTEGER :: USE_tilt                              ! ????
   INTEGER :: USE_testing
   REAL(KIND=dp) :: usetilt, balloon, noballoon, norm ! helpers

   REAL(KIND=dp) :: C1                              ! q1 = C1*Q1,  Q1 = rho_tor(psi)!(not anymore)psi-psi0
   REAL(KIND=dp) :: C2                              ! q2 = C2*Q2,  Q2 = zeta-nu
   REAL(KIND=dp) :: C3                              ! q3 = C3*Q3,  Q3 = chi, nu=qpsi*chi

CONTAINS

!========================================================================================================
! PART 2: SUBROUTINE parameters(run)
! Sets values for the parameters in the current run.
! Should be called in main at the start of every run.
!========================================================================================================

   SUBROUTINE parameters(run)
      USE sims
      IMPLICIT NONE

      INTEGER, INTENT(IN)                 :: run      ! simulation index in values.dat
      REAL(KIND=dp), DIMENSION(cols1)     :: pp       ! parameter vector for this run
      CHARACTER(LEN=10000)                :: CWD
      CHARACTER(LEN=35)                   :: aux
      INTEGER                             :: lenCWD, i

      pp = array1(run, :)

      !---------------------------------------------------------------------------------
      ! Import parameters
      !---------------------------------------------------------------------------------
t0       = pp(param_index("t0"))
tc       = pp(param_index("tc"))
tt       = pp(param_index("tt"))
tmax     = pp(param_index("tmax"))

Np       = INT(pp(param_index("Np")))
Nc       = INT(pp(param_index("Nc")))
Nloop    = INT(pp(param_index("Nloop")))
Nreal    = INT(pp(param_index("Nreal")))
Nt       = INT(pp(param_index("Nt")))
ntraj    = INT(pp(param_index("ntraj")))

Ti       = pp(param_index("Ti"))
Te       = pp(param_index("Te"))
B0       = pp(param_index("B0"))
R0       = pp(param_index("R0"))
Ln       = pp(param_index("Ln"))
Li       = pp(param_index("Li"))
Le       = pp(param_index("Le"))
Zeff     = pp(param_index("Zeff"))
Aeff     = pp(param_index("Aeff"))
ndens    = pp(param_index("ndens"))
a0       = pp(param_index("a0"))

Omgt0    = pp(param_index("Omgt0"))
Omgtprim = pp(param_index("Omgtprim"))

safe1    = pp(param_index("safe1"))
safe2    = pp(param_index("safe2"))
safe3    = pp(param_index("safe3"))

amp      = pp(param_index("amp"))
elong    = pp(param_index("elong"))

device    = INT(pp(param_index("device")))
shot      = INT(pp(param_index("shot")))
shotslice = INT(pp(param_index("shotslice")))

Phi      = pp(param_index("Phi"))
turbprof = pp(param_index("turbprof"))
Ai       = pp(param_index("Ai"))
lambdax  = pp(param_index("lambdax"))
lambday  = pp(param_index("lambday"))
lambdaz  = pp(param_index("lambdaz"))
lbalonz  = pp(param_index("lbalonz"))
tauc     = pp(param_index("tauc"))
k0i      = pp(param_index("k0i"))
k0e      = pp(param_index("k0e"))
gamma_ZF = pp(param_index("gamma_ZF"))

X0       = pp(param_index("X0"))
Y0       = pp(param_index("Y0"))
Z0       = pp(param_index("Z0"))
Tw       = pp(param_index("Tw"))
Ew       = pp(param_index("Ew"))
pitch    = pp(param_index("pitch"))
Aw       = pp(param_index("Aw"))
Zw       = pp(param_index("Zw"))

magnetic_model = INT(pp(param_index("magnetic_model")))
USE_larmor     = INT(pp(param_index("USE_larmor")))
position_type  = INT(pp(param_index("position_type")))
pitch_type     = INT(pp(param_index("pitch_type")))
energy_type    = INT(pp(param_index("energy_type")))
USE_coll       = INT(pp(param_index("USE_coll")))
USE_turb       = INT(pp(param_index("USE_turb")))
USE_magnturb   = INT(pp(param_index("USE_magnturb")))
USE_freq       = INT(pp(param_index("USE_freq")))
USE_polar      = INT(pp(param_index("USE_polar")))
USE_PC         = INT(pp(param_index("USE_PC")))
USE_real       = INT(pp(param_index("USE_real")))
USE_corr       = INT(pp(param_index("USE_corr")))
USE_balloon    = INT(pp(param_index("USE_balloon")))
USE_tilt       = INT(pp(param_index("USE_tilt")))
USE_testing    = INT(pp(param_index("USE_testing")))


      ! End of import parameters
      Nqua = 16   ! number of quantities stored on the EFIT-like grid

      !---------------------------------------------------------------------------------
      ! EFIT file & data import
      !---------------------------------------------------------------------------------
      IF ((magnetic_model == 1) .OR. (magnetic_model == 2)) THEN

         CALL GetCWD(CWD)
         efit_file = TRIM(ADJUSTL(CWD))
         lenCWD    = LEN(efit_file)

         WRITE(aux, '(I0)') shot

         ! Device path
         if (device == 1) then
            efit_file = efit_file(1:lenCWD - 7)//'/G_EQDSK/MAST-U/'//trim(aux)//'/QQ_'
         elseif (device == 2) then
            efit_file = efit_file(1:lenCWD - 7)//'/G_EQDSK/TCV/'//trim(aux)//'/QQ_'
         elseif (device == 3) then
            efit_file = efit_file(1:lenCWD - 7)//'/G_EQDSK/WEST/'//trim(aux)//'/QQ_'
         elseif (device == 4) then
            efit_file = efit_file(1:lenCWD - 7)//'/G_EQDSK/JTSA/'//trim(aux)//'/QQ_'
         else
            WRITE(*,*) 'the specified tokamak does not exist'
         end if

         WRITE(aux, '(I0)') shotslice

         ! Correct naming for WEST
         if ((device == 3) .and. (shotslice < 10)) then
            aux = '00'//trim(aux)
         elseif ((device == 3) .and. (shotslice < 100) .and. (shotslice >= 10)) then
            aux = '0'//trim(aux)
         end if

         efit_file = efit_file//trim(aux)//'.dat'

         OPEN (UNIT=100, FILE=trim(efit_file), STATUS='old', ACTION='read')

! a thought here: the direction of the magnetic field and the twisting of the field lines is controlled by B0 and q(r); if you set them +/- you can controll eveything, so 
! you must be carefull at the input parameters; from the EFIT files there is no problem;


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
      END IF

      !---------------------------------------------------------------------------------
      ! Processing and evaluation of derived quantities
      !---------------------------------------------------------------------------------
      Ti    = Ti*1000.0_dp                     ! imported in keV
      Te    = Te*1000.0_dp                     ! imported in keV
      ndens = ndens*10.0_dp**19                ! imported in 1e19 m^-3

      vth  = sqrt(Ti*q0/mi)
      rhoi = mi*vth/q0/abs(B0)
      wi   = q0*abs(B0)/mi

      stepR = (maxR - minR) / (NgridR - 1)
      stepZ = (maxZ - minZ) / (NgridZ - 1)

      r00 = sqrt((X0 - R0)**2 + Y0**2)

      !---------------------------------------------------------------------------------
      ! Scaling quantities
      !---------------------------------------------------------------------------------
      IF ((magnetic_model == 1) .OR. (magnetic_model == 2)) THEN
         minR  = minR/R0
         maxR  = maxR/R0
         minZ  = minZ/R0
         maxZ  = maxZ/R0
         stepR = stepR/R0
         stepZ = stepZ/R0

         ! Efit_data ::
         ! {psi,psir,psiz,psirr,psirz,psizz, F,Fprim, qpsi,qprim, rhot,rhotR,rhotZ, chi,chiR,chiZ}
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
      END IF

      gama = 1.0_dp - 2.0_dp*sqrt((a0/R0)**2 - (a0/R0)**4)
      alfa = elong*(sqrt(2.0_dp - gama) - sqrt(gama)) / sqrt(2.0_dp - 2.0_dp*gama)
      amp  = 0.5_dp*(1.0_dp + alfa**2) / (1.0_dp - gama) / safe1

      a0    = a0/R0
      Omgt0 = Omgt0*10.0_dp**3 / (vth/R0)

      r00 = r00/R0
      X0  = X0/R0
      Y0  = Y0/R0

      CALL OnlyQ(X0, Y0, q00, psi0)
      CALL Onlycoord_new(X0, Y0, Z0, q10) !asta e pentru deltax

      C1 = a0*R0/rhoi
      C2 = -r00/q00 * R0/rhoi
      C3 = 1.0_dp

      Nci = int(Nc*Ai**2)
      Nce = Nc - Nci

      Ae = sqrt(1.0_dp - Ai**2)

      c0 = logL*ndens/8.0_dp/pi * (Zeff*Zw*q0**2/(eps*mi*Aw))**2

      delta = 5.0_dp*4.0_dp/3.0_dp/sqrt(pi) * (tmax - t0)/Nt * R0*2.0_dp*logL*ndens*(Zeff*Zw*q0**2)**2 / &
              (pi*eps**2*mi**2*Aw**2*vth**4*sqrt(2.0_dp/Aeff)**3)

      CALL error_signal

      gamma_E = gamma_ZF + 0.0_dp*Omgtprim    ! sheering
      usetilt = real(USE_tilt, dp)
      balloon = real(USE_balloon, dp)                ! integer -> 0.0 or 1.0
      noballoon = 1.0_dp - balloon
      norm = 1.0_dp / sqrt(real(Nc, dp)/2.0_dp)

   END SUBROUTINE parameters

END MODULE constants

