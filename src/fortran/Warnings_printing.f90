!========================================================================================================
! File: Warnings_printing.f90
! Console/header output and run-parameter validation.
!========================================================================================================

!--------------------------------------------------------------------------------------------------------
! Routine: print_run_header
! Purpose: Print concise run metadata to stdout and to the parameter-output unit.
!--------------------------------------------------------------------------------------------------------
SUBROUTINE print_run_header(runs, nsim, date_time, folder, lene, simf, simi, iono)
   USE sims
   USE constants
   IMPLICIT NONE

   ! I/O variables
   CHARACTER(LEN=35), INTENT(IN) :: runs, nsim
   INTEGER,            INTENT(IN) :: lene
   CHARACTER(LEN=lene),INTENT(IN) :: folder
   CHARACTER(LEN=18),  INTENT(IN) :: date_time
   INTEGER,            INTENT(IN) :: iono, simf, simi

   real :: timpnufigrabit1, timpnufigrabit2

   WRITE (*, *) 'Sim number  🔥 :: ', noofsimstr
   WRITE (*, *) 'Run number  ⚡ :: ', TRIM(ADJUSTl(runs)), '/', TRIM(ADJUSTl(nsim))
   WRITE (*, *) 'Start time  🕒 :: ', date_time(7:8), '/', date_time(5:6), '/', date_time(1:4), &
                                 &'  ', date_time(9:10), ':', date_time(11:12), ':', date_time(13:14)

   timpnufigrabit1 = real(Np,wp)*real(Nloop,wp)*real(Nreal,wp)*real(Nt,wp)*(1.2_wp + 0.7_wp*real(USE_coll,wp) + 3.9_wp*real(USE_turb,wp) + 1.5_wp*real(USE_coll,wp)*real(USE_turb,wp))/10.0**8
   timpnufigrabit2 = (simf-simi+1)*timpnufigrabit1
   WRITE (*, *) 'EsCompTime  ⏳ ::', timpnufigrabit1, 's  out of a total of ', timpnufigrabit2, ' s'
   WRITE (*, '(A)') ' Folder      📁 :: '//TRIM(folder)

   WRITE (iono, *) 'Run number :: ', TRIM(ADJUSTl(runs)), '/', TRIM(ADJUSTl(nsim))
   WRITE (iono, *) 'Start time :: ', date_time(7:8), '/', date_time(5:6), '/', date_time(1:4), &
                                &'  ', date_time(9:10), ':', date_time(11:12), ':', date_time(13:14)
   WRITE (iono, *) 'Folder     :: ', folder


END SUBROUTINE print_run_header


!--------------------------------------------------------------------------------------------------------
! Routine: validate_raw_run_parameters
! Purpose: Validate imported parameters before unit conversion, EFIT import, and
!          derived quantities that may divide by them.
!--------------------------------------------------------------------------------------------------------
SUBROUTINE validate_raw_run_parameters
   USE constants
   USE, INTRINSIC :: iso_fortran_env, ONLY: error_unit, output_unit
   IMPLICIT NONE

   INTEGER :: nerr
   REAL(KIND=wp), PARAMETER :: tiny_dp = 1.0e-30_wp

   nerr = 0

   ! --- time / stepping ---
   IF (t0 >= tmax) THEN
      WRITE(error_unit,'(A,2(1X,ES12.4))') 'Error :: require t0 < tmax, got (t0,tmax)=', t0, tmax
      nerr = nerr! + 1
   END IF

   IF (tc > tmax) THEN
      WRITE(error_unit,'(A,3(1X,ES12.4))') 'Error :: tc must be <= tmax, got (t0,tc,tmax)=', t0, tc, tmax
      nerr = nerr! + 1
   END IF

   IF (tt > tmax) THEN
      WRITE(error_unit,'(A,3(1X,ES12.4))') 'Error :: tt must be <= tmax, got (t0,tt,tmax)=', t0, tt, tmax
      nerr = nerr! + 1
   END IF

   IF (Nt <= 0) THEN
      WRITE(error_unit,'(A,I0)') 'Error :: Nt must be > 0, got Nt = ', Nt
      nerr = nerr + 1
   END IF

   ! --- sizes / counts ---
   IF (Np <= 0) THEN
      WRITE(error_unit,'(A,I0)') 'Error :: Np must be > 0, got Np = ', Np
      nerr = nerr + 1
   END IF

   IF (Nc <= 0) THEN
      WRITE(error_unit,'(A,I0)') 'Error :: Nc must be > 0, got Nc = ', Nc
      nerr = nerr + 1
   END IF

   IF (Nloop <= 0) THEN
      WRITE(error_unit,'(A,I0)') 'Error :: Nloop must be > 0, got Nloop = ', Nloop
      nerr = nerr + 1
   END IF

   IF (Nreal <= 0) THEN
      WRITE(error_unit,'(A,I0)') 'Error :: Nreal must be > 0, got Nreal = ', Nreal
      nerr = nerr + 1
   END IF

   IF (ntraj < 0 .OR. ntraj > Np) THEN
      WRITE(error_unit,'(A,2(1X,I0))') 'Error :: ntraj must satisfy 0 <= ntraj <= Np, got (ntraj,Np)=', ntraj, Np
      nerr = nerr + 1
   END IF

   ! --- physical parameters before normalization ---
   IF (magnetic_model == 3 .OR. magnetic_model == 4) THEN
      IF (R0 <= tiny_dp) THEN
         WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: R0 must be > 0 before normalization, got R0 =', R0
         nerr = nerr + 1
      END IF

      IF (a0 <= 0.0_wp) THEN
         WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: a0 must be > 0 before normalization, got a0 =', a0
         nerr = nerr + 1
      END IF

      IF (a0 >= R0) THEN
         WRITE(error_unit,'(A,2(1X,ES12.4))') 'Error :: require a0 < R0 before Solovev/circular normalization, got (a0,R0)=', a0, R0
         nerr = nerr + 1
      END IF

      IF (ABS(B0) <= tiny_dp) THEN
         WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: |B0| must be > 0, got B0 =', B0
         nerr = nerr + 1
      END IF
   END IF

   IF (Ti <= 0.0_wp .OR. Te <= 0.0_wp) THEN
      WRITE(error_unit,'(A,2(1X,ES12.4))') 'Error :: Ti and Te must be > 0, got (Ti,Te)=', Ti, Te
      nerr = nerr + 1
   END IF

   IF (ndens <= 0.0_wp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: ndens must be > 0, got ndens =', ndens
      nerr = nerr + 1
   END IF

   IF (As <= 0.0_wp .OR. Aeff <= 0.0_wp) THEN
      WRITE(error_unit,'(A,2(1X,ES12.4))') 'Error :: As and Aeff must be > 0, got (As,Aeff)=', As, Aeff
      nerr = nerr + 1
   END IF

   IF (Zeff <= 0.0_wp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: Zeff must be > 0, got Zeff =', Zeff
      nerr = nerr + 1
   END IF

   IF (USE_PC == ON .AND. ABS(Zs) <= tiny_dp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: Zs must be nonzero when USE_PC=ON, got Zs =', Zs
      nerr = nerr + 1
   END IF

   IF (ABS(s1) <= tiny_dp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: s1 must be nonzero; amp is derived using 1/s1, got s1 =', s1
      nerr = nerr + 1
   END IF

   IF (elong <= 0.0_wp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: elong must be > 0, got elong =', elong
      nerr = nerr + 1
   END IF

   ! --- model choices ---
   SELECT CASE (magnetic_model)
   CASE (1, 2, 3, 4)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0)') 'Error :: magnetic_model not implemented, got magnetic_model = ', magnetic_model
      nerr = nerr + 1
   END SELECT

   SELECT CASE (device)
   CASE (1, 2, 3, 4)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0,A)') 'Error :: device = ', device, ' not implemented (1=MAST,2=TCV,3=WEST,4=JTSA).'
      nerr = nerr + 1
   END SELECT

   SELECT CASE (turb_model)
   CASE (1, 2)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0)') 'Error :: turb_model unknown, got turb_model = ', turb_model
      nerr = nerr + 1
   END SELECT

   SELECT CASE (x_corr)
   CASE (1, 2, 3)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0)') 'Error :: x_corr unknown, got x_corr = ', x_corr
      nerr = nerr + 1
   END SELECT

   SELECT CASE (y_corr)
   CASE (1, 2)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0)') 'Error :: y_corr unknown, got y_corr = ', y_corr
      nerr = nerr + 1
   END SELECT

   SELECT CASE (t_corr)
   CASE (1, 2)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0)') 'Error :: t_corr unknown, got t_corr = ', t_corr
      nerr = nerr + 1
   END SELECT

   ! --- turbulence and marker parameters ---
   IF (Ai < 0.0_wp .OR. Ai > 1.0_wp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: Ai must be in [0,1], got Ai =', Ai
      nerr = nerr + 1
   END IF

   IF (Phi < 0.0_wp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: Phi must be >= 0, got Phi =', Phi
      nerr = nerr + 1
   END IF

   IF (lbalonz <= 0.0_wp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: lbalonz must be > 0; dmmax divides by it, got lbalonz =', lbalonz
      nerr = nerr + 1
   END IF

   IF (USE_turb == ON) THEN
      IF (lambdax <= 0.0_wp .OR. lambday <= 0.0_wp .OR. lambdaz <= 0.0_wp) THEN
         WRITE(error_unit,'(A,3(1X,ES12.4))') 'Error :: turbulence lengths must be > 0, got (lambdax,lambday,lambdaz)=', &
              lambdax, lambday, lambdaz
         nerr = nerr + 1
      END IF
      IF (tauc <= 0.0_wp) THEN
         WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: tauc must be > 0 when turbulence is enabled, got tauc =', tauc
         nerr = nerr + 1
      END IF
      IF (k0i <= 0.0_wp .OR. k0e <= 0.0_wp) THEN
         WRITE(error_unit,'(A,2(1X,ES12.4))') 'Error :: k0i and k0e must be > 0 when turbulence is enabled, got (k0i,k0e)=', k0i, k0e
         nerr = nerr + 1
      END IF
   END IF

   IF (Ts <= 0.0_wp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: Ts must be > 0, got Ts =', Ts
      nerr = nerr + 1
   END IF

   IF (energy_type == 1 .AND. Es < 0.0_wp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: Es must be >= 0 for fixed-energy initialization, got Es =', Es
      nerr = nerr + 1
   END IF

   IF (pitch_type == 1 .AND. ABS(pitch) > 1.0_wp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: fixed pitch must be in [-1,1], got pitch =', pitch
      nerr = nerr + 1
   END IF

   IF (annulus_width < 0.0_wp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: annulus_width must be >= 0, got annulus_width =', annulus_width
      nerr = nerr + 1
   END IF

   CALL validate_switch(USE_larmor,   'USE_larmor',   nerr)
   CALL validate_switch(USE_coll,     'USE_coll',     nerr)
   CALL validate_switch(USE_turb,     'USE_turb',     nerr)
   CALL validate_switch(USE_magnturb, 'USE_magnturb', nerr)
   CALL validate_switch(USE_freq,     'USE_freq',     nerr)
   CALL validate_switch(USE_polar,    'USE_polar',    nerr)
   CALL validate_switch(USE_PC,       'USE_PC',       nerr)
   CALL validate_switch(USE_real,     'USE_real',     nerr)
   CALL validate_switch(USE_corr,     'USE_corr',     nerr)
   CALL validate_switch(USE_balloon,  'USE_balloon',  nerr)
   CALL validate_switch(USE_tilt,     'USE_tilt',     nerr)
   CALL validate_switch(USE_testing,  'USE_testing',  nerr)

   SELECT CASE (position_type)
   CASE (1, 2, 3)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0,A)') 'Error :: position_type unknown or unimplemented, got ', position_type, ' (implemented: 1,2,3).'
      nerr = nerr + 1
   END SELECT

   SELECT CASE (pitch_type)
   CASE (1, 2)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0)') 'Error :: pitch_type unknown, got pitch_type = ', pitch_type
      nerr = nerr + 1
   END SELECT

   SELECT CASE (energy_type)
   CASE (1, 2)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0)') 'Error :: energy_type unknown, got energy_type = ', energy_type
      nerr = nerr + 1
   END SELECT

   CALL stop_if_errors(nerr, 'raw input')

   IF (Np < 100) THEN
      WRITE(output_unit,'(A,I0,A)') 'Warning :: Np is small (Np = ', Np, '); results may be noisy.'
   END IF

   IF (Nc < 50) THEN
      WRITE(output_unit,'(A,I0,A)') 'Warning :: Nc is small (Nc = ', Nc, '); spectrum may be under-resolved.'
   END IF

   IF (USE_turb == ON .AND. Phi == 0.0_wp) THEN
      WRITE(output_unit,'(A)') 'Warning :: USE_turb=ON but Phi=0; turbulent potential amplitude is zero.'
   END IF

END SUBROUTINE validate_raw_run_parameters

!--------------------------------------------------------------------------------------------------------
! Routine: validate_derived_run_parameters
! Purpose: Validate values derived after unit conversion, EFIT import, and reference
!          equilibrium queries.
!--------------------------------------------------------------------------------------------------------
SUBROUTINE validate_derived_run_parameters
   USE constants
   USE, INTRINSIC :: iso_fortran_env, ONLY: error_unit, output_unit
   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
   IMPLICIT NONE

   INTEGER :: nerr
   REAL(KIND=wp) :: dt
   REAL(KIND=wp), PARAMETER :: tiny_dp = 1.0e-30_wp

   nerr = 0
   dt = (tmax - t0)/REAL(Nt,wp)

   ! --- normalized geometry and reference position ---
   IF (R0 <= tiny_dp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: derived/imported R0 must be > 0, got R0 =', R0
      nerr = nerr + 1
   END IF

   IF (ABS(B0) <= tiny_dp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: derived/imported |B0| must be > 0, got B0 =', B0
      nerr = nerr + 1
   END IF

   IF (a0 <= 0.0_wp .OR. a0 >= 1.0_wp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: normalized a0 must be inside (0,1), got a0 =', a0
      nerr = nerr + 1
   END IF

   IF (r00 > a0) THEN
      WRITE(error_unit,'(A,2(1X,ES12.4))') 'Error :: reference position outside plasma (r00 > a0). (r00,a0)=', r00, a0
      nerr = nerr + 1
   END IF

   IF (ABS(q00) <= tiny_dp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: q00 must be nonzero; C2 divides by q00, got q00 =', q00
      nerr = nerr + 1
   END IF

   IF (vth <= 0.0_wp .OR. rhoi <= 0.0_wp .OR. wi <= 0.0_wp) THEN
      WRITE(error_unit,'(A,3(1X,ES12.4))') 'Error :: derived (vth,rhoi,wi) must be > 0, got ', vth, rhoi, wi
      nerr = nerr + 1
   END IF

   IF (.NOT. ieee_is_finite(vth) .OR. .NOT. ieee_is_finite(rhoi) .OR. .NOT. ieee_is_finite(q00) .OR. &
       .NOT. ieee_is_finite(psi0) .OR. .NOT. ieee_is_finite(C1) .OR. .NOT. ieee_is_finite(C2)) THEN
      WRITE(error_unit,'(A)') 'Error :: one or more derived normalization/reference quantities are NaN or Inf.'
      nerr = nerr + 1
   END IF

   ! --- EFIT-derived grid sanity ---
   IF (magnetic_model == 1 .OR. magnetic_model == 2) THEN
      IF (NgridR <= 1) THEN
         WRITE(error_unit,'(A,I0)') 'Error :: NgridR must be > 1 for EFIT, got NgridR = ', NgridR
         nerr = nerr + 1
      END IF
      IF (NgridZ <= 1) THEN
         WRITE(error_unit,'(A,I0)') 'Error :: NgridZ must be > 1 for EFIT, got NgridZ = ', NgridZ
         nerr = nerr + 1
      END IF
      IF (Nqua /= 16) THEN
         WRITE(error_unit,'(A,I0)') 'Error :: EFIT import expects Nqua=16, got Nqua = ', Nqua
         nerr = nerr + 1
      END IF
      IF (maxR <= minR) THEN
         WRITE(error_unit,'(A,2(1X,ES12.4))') 'Error :: EFIT grid invalid: maxR <= minR; (minR,maxR)=', minR, maxR
         nerr = nerr + 1
      END IF
      IF (maxZ <= minZ) THEN
         WRITE(error_unit,'(A,2(1X,ES12.4))') 'Error :: EFIT grid invalid: maxZ <= minZ; (minZ,maxZ)=', minZ, maxZ
         nerr = nerr + 1
      END IF
      IF (stepR <= 0.0_wp .OR. stepZ <= 0.0_wp) THEN
         WRITE(error_unit,'(A,2(1X,ES12.4))') 'Error :: EFIT steps must be > 0, got (stepR,stepZ)=', stepR, stepZ
         nerr = nerr + 1
      END IF
      IF (.NOT. ALLOCATED(Efit_data)) THEN
         WRITE(error_unit,'(A)') 'Error :: EFIT model selected but Efit_data is not allocated.'
         nerr = nerr + 1
      ELSEIF (SIZE(Efit_data) /= NgridR*NgridZ*Nqua) THEN
         WRITE(error_unit,'(A,2(1X,I0))') 'Error :: Efit_data size mismatch, got/expected =', SIZE(Efit_data), NgridR*NgridZ*Nqua
         nerr = nerr + 1
      END IF
   END IF

   ! --- derived turbulence/species quantities ---
   IF (Nci < 0 .OR. Nce < 0 .OR. Nci + Nce /= Nc) THEN
      WRITE(error_unit,'(A,3(1X,I0))') 'Error :: invalid derived mode counts (Nci,Nce,Nc)=', Nci, Nce, Nc
      nerr = nerr + 1
   END IF

   IF (.NOT. ieee_is_finite(Ae) .OR. Ae < 0.0_wp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: invalid derived Ae =', Ae
      nerr = nerr + 1
   END IF

   CALL stop_if_errors(nerr, 'derived')

   IF (dt > 1.0_wp) THEN
      WRITE(output_unit,'(A,1X,ES12.4)') 'Warning :: dt is large (dt = (tmax-t0)/Nt) =', dt
   END IF

   IF (USE_turb == ON .AND. tauc > 0.0_wp .AND. dt > 0.5_wp*tauc) THEN
      WRITE(output_unit,'(A,2(1X,ES12.4))') 'Warning :: dt is large relative to tauc, got (dt,tauc)=', dt, tauc
   END IF

   IF (10.0_wp < ABS(B0)) THEN
      WRITE(output_unit,'(A,1X,ES12.4)') 'Warning :: |B0| is large (B0 =', B0, '). Check normalization/units.'
   END IF

   IF (Phi > 0.2_wp) THEN
      WRITE(output_unit,'(A,1X,ES12.4)') 'Warning :: turbulence amplitude Phi is high (Phi =', Phi, ').'
   END IF

END SUBROUTINE validate_derived_run_parameters

!--------------------------------------------------------------------------------------------------------
! Routine: validate_switch
! Purpose: Validate an integer OFF/ON switch and append to the caller's error count.
!--------------------------------------------------------------------------------------------------------
SUBROUTINE validate_switch(value, name, nerr)
   USE constants
   USE, INTRINSIC :: iso_fortran_env, ONLY: error_unit
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: value
   CHARACTER(LEN=*), INTENT(IN) :: name
   INTEGER, INTENT(INOUT) :: nerr

   SELECT CASE (value)
   CASE (OFF, ON)
   CASE DEFAULT
      WRITE(error_unit,'(A,A,A,I0,A)') 'Error :: ', TRIM(name), ' must be OFF/ON (0/1), got ', value, '.'
      nerr = nerr + 1
   END SELECT
END SUBROUTINE validate_switch

!--------------------------------------------------------------------------------------------------------
! Routine: stop_if_errors
! Purpose: Emit a stage-specific validation summary and stop when fatal errors exist.
!--------------------------------------------------------------------------------------------------------
SUBROUTINE stop_if_errors(nerr, stage)
   USE, INTRINSIC :: iso_fortran_env, ONLY: error_unit
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: nerr
   CHARACTER(LEN=*), INTENT(IN) :: stage

   IF (nerr > 0) THEN
      WRITE(error_unit,'(A,I0,A,A,A)') 'Aborting: found ', nerr, ' fatal ', TRIM(stage), ' parameter error(s).'
      ERROR STOP "Invalid parameters"
   END IF
END SUBROUTINE stop_if_errors
