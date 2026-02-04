!      ! Subroutine :: Print
!      ! Purpose    :: Printing some info
!  ========================================================================================================================================================
!      ! Record of revisions:
!      !    Date       |    Programmer          |   Description of change
!      !  =======      |   =============        |   =========================
!      !  01/01/2024   |     D. I. Palade       |   Remake
!      !  --/--/2026   |     D. I. Palade       | Continuous improvements
!  ========================================================================================================================================================

SUBROUTINE printing(runs, nsim, date_time, folder, lene, simf, simi, iono)
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

   timpnufigrabit1 = real(Np,dp)*real(Nloop,dp)*real(Nreal,dp)*real(Nt,dp)*(1.2_dp + 0.7_dp*real(USE_coll,dp) + 3.9_dp*real(USE_turb,dp) + 1.5_dp*real(USE_coll,dp)*real(USE_turb,dp))/10.0**8
   timpnufigrabit2 = (simf-simi+1)*timpnufigrabit1
   WRITE (*, *) 'EsCompTime  ⏳ ::', timpnufigrabit1, 's  out of a total of ', timpnufigrabit2, ' s'
   WRITE (*, '(A)') ' Folder      📁 :: '//TRIM(folder)

   WRITE (iono, *) 'Run number :: ', TRIM(ADJUSTl(runs)), '/', TRIM(ADJUSTl(nsim))
   WRITE (iono, *) 'Start time :: ', date_time(7:8), '/', date_time(5:6), '/', date_time(1:4), &
                                &'  ', date_time(9:10), ':', date_time(11:12), ':', date_time(13:14)
   WRITE (iono, *) 'Folder     :: ', folder


END SUBROUTINE printing


! ======      errors   ==============
! ====== errors / warnings ==========================================
SUBROUTINE error_signal
   USE constants
   USE, INTRINSIC :: iso_fortran_env, ONLY: error_unit, output_unit
   IMPLICIT NONE

   INTEGER :: nerr
   REAL(KIND=dp), PARAMETER :: tiny_dp = 1.0e-30_dp

   nerr = 0

   ! ================================================================
   ! FATAL ERRORS (collect all, then stop once)
   ! ================================================================

   ! --- time / stepping ---
   IF (t0 > tmax) THEN
      WRITE(error_unit,'(A)') 'Error :: t0 > tmax (initial time greater than final time).'
      nerr = nerr + 1
   END IF

   IF (Nt <= 0) THEN
      WRITE(error_unit,'(A,I0)') 'Error :: Nt must be > 0, got Nt = ', Nt
      nerr = nerr + 1
   END IF

   IF (tmax <= t0) THEN
      WRITE(error_unit,'(A,2(1X,ES12.4))') 'Error :: non-positive simulation interval (t0,tmax)=', t0, tmax
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

   ! --- physical sanity that can cause NaNs/div-by-zero in derived quantities ---
   IF (R0 <= tiny_dp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: R0 must be > 0, got R0 =', R0
      nerr = nerr + 1
   END IF

   IF (a0 <= 0.0_dp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: a0 must be > 0 (normalized minor radius), got a0 =', a0
      nerr = nerr + 1
   END IF

   IF (ABS(B0) <= tiny_dp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: |B0| must be > 0, got B0 =', B0
      nerr = nerr + 1
   END IF

   IF (Ti <= 0.0_dp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: Ti must be > 0, got Ti =', Ti
      nerr = nerr + 1
   END IF

   IF (Te <= 0.0_dp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: Te must be > 0, got Te =', Te
      nerr = nerr + 1
   END IF

   IF (ndens <= 0.0_dp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: ndens must be > 0, got ndens =', ndens
      nerr = nerr + 1
   END IF

   IF (Aw <= 0.0_dp .OR. Aeff <= 0.0_dp) THEN
      WRITE(error_unit,'(A,2(1X,ES12.4))') 'Error :: Aw and Aeff must be > 0, got (Aw,Aeff)=', Aw, Aeff
      nerr = nerr + 1
   END IF

   IF (Zeff <= 0.0_dp .OR. Zw <= 0.0_dp) THEN
      WRITE(error_unit,'(A,2(1X,ES12.4))') 'Error :: Zeff and Zw must be > 0, got (Zeff,Zw)=', Zeff, Zw
      nerr = nerr + 1
   END IF

   ! --- EFIT-derived grid sanity (only if model uses EFIT data) ---
   SELECT CASE (magnetic_model)
   CASE (1, 2)
      IF (NgridR <= 1) THEN
         WRITE(error_unit,'(A,I0)') 'Error :: NgridR must be > 1 for EFIT, got NgridR = ', NgridR
         nerr = nerr + 1
      END IF
      IF (NgridZ <= 1) THEN
         WRITE(error_unit,'(A,I0)') 'Error :: NgridZ must be > 1 for EFIT, got NgridZ = ', NgridZ
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
      IF (.NOT. ALLOCATED(Efit_data)) THEN
         WRITE(error_unit,'(A)') 'Error :: EFIT model selected but Efit_data is not allocated.'
         nerr = nerr + 1
      END IF
   CASE (3, 4)
      ! Solovev / Circular: no EFIT grid requirements here
   CASE DEFAULT
      WRITE(error_unit,'(A,I0)') 'Error :: magnetic_model not implemented, got magnetic_model = ', magnetic_model
      nerr = nerr + 1
   END SELECT

   ! --- device choice ---
   SELECT CASE (device)
   CASE (1, 2, 3, 4)
      ! ok
   CASE DEFAULT
      WRITE(error_unit,'(A,I0,A)') 'Error :: device = ', device, ' not implemented (1=MAST,2=TCV,3=WEST,4=JTSA).'
      nerr = nerr + 1
   END SELECT

   ! --- turbulence fractions / amplitudes ---
   IF (Ai < 0.0_dp .OR. Ai > 1.0_dp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: Ai must be in [0,1], got Ai =', Ai
      nerr = nerr + 1
   END IF

   IF (Phi < 0.0_dp) THEN
      WRITE(error_unit,'(A,1X,ES12.4)') 'Error :: Phi must be >= 0, got Phi =', Phi
      nerr = nerr + 1
   END IF

   ! --- switches / enums ---
   SELECT CASE (USE_larmor)
   CASE (OFF, ON)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0,A)') 'Error :: USE_larmor must be OFF/ON (0/1), got ', USE_larmor, '.'
      nerr = nerr + 1
   END SELECT

   SELECT CASE (USE_freq)
   CASE (OFF, ON)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0,A)') 'Error :: USE_freq must be OFF/ON (0/1), got ', USE_freq, '.'
      nerr = nerr + 1
   END SELECT

   SELECT CASE (USE_coll)
   CASE (OFF, ON)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0,A)') 'Error :: USE_coll must be OFF/ON (0/1), got ', USE_coll, '.'
      nerr = nerr + 1
   END SELECT

   SELECT CASE (USE_balloon)
   CASE (OFF, ON)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0,A)') 'Error :: USE_balloon must be OFF/ON (0/1), got ', USE_balloon, '.'
      nerr = nerr + 1
   END SELECT

   SELECT CASE (USE_tilt)
   CASE (OFF, ON)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0,A)') 'Error :: USE_tilt must be OFF/ON (0/1), got ', USE_tilt, '.'
      nerr = nerr + 1
   END SELECT

SELECT CASE (USE_turb)
   CASE (OFF, ON)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0,A)') 'Error :: USE_turb must be OFF/ON (0/1), got ', USE_turb, '.'
      nerr = nerr + 1
   END SELECT

   SELECT CASE (USE_magnturb)
   CASE (OFF, ON)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0,A)') 'Error :: USE_magnturb must be OFF/ON (0/1), got ', USE_magnturb, '.'
      nerr = nerr + 1
   END SELECT

   SELECT CASE (USE_polar)
   CASE (OFF, ON)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0,A)') 'Error :: USE_polar must be OFF/ON (0/1), got ', USE_polar, '.'
      nerr = nerr + 1
   END SELECT

   SELECT CASE (USE_PC)
   CASE (OFF, ON)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0,A)') 'Error :: USE_PC must be OFF/ON (0/1), got ', USE_PC, '.'
      nerr = nerr + 1
   END SELECT

   SELECT CASE (USE_real)
   CASE (OFF, ON)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0,A)') 'Error :: USE_real must be OFF/ON (0/1), got ', USE_real, '.'
      nerr = nerr + 1
   END SELECT

   SELECT CASE (USE_corr)
   CASE (OFF, ON)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0,A)') 'Error :: USE_corr must be OFF/ON (0/1), got ', USE_corr, '.'
      nerr = nerr + 1
   END SELECT

   SELECT CASE (position_type)
   CASE (1, 2, 3, 4, 5)
   CASE DEFAULT
      WRITE(error_unit,'(A,I0)') 'Error :: position_type unknown, got position_type = ', position_type
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

   ! --- basic placement sanity in normalized units ---
   IF (r00 > a0) THEN
      WRITE(error_unit,'(A,2(1X,ES12.4))') 'Error :: particles placed outside plasma (r00 > a0). (r00,a0)=', r00, a0
      nerr = nerr + 1
   END IF

   ! Stop once if any fatal errors were found.
   IF (nerr > 0) THEN
      WRITE(error_unit,'(A,I0,A)') 'Aborting: found ', nerr, ' fatal input/parameter error(s).'
      ERROR STOP "Invalid parameters"
   END IF

   ! ================================================================
   ! WARNINGS (print all)
   ! ================================================================

   IF (Np < 100) THEN
      WRITE(output_unit,'(A,I0,A)') 'Warning :: Np is small (Np = ', Np, '); results may be noisy.'
   END IF

   IF (Nc < 50) THEN
      WRITE(output_unit,'(A,I0,A)') 'Warning :: Nc is small (Nc = ', Nc, '); spectrum may be under-resolved.'
   END IF

   IF (Nt > 0) THEN
      IF ((tmax - t0)/REAL(Nt,dp) > 1.0_dp) THEN
         WRITE(output_unit,'(A,1X,ES12.4)') 'Warning :: dt is large (dt = (tmax-t0)/Nt) =', (tmax - t0)/REAL(Nt,dp)
      END IF
   END IF

   IF (10.0_dp < ABS(B0)) THEN
      WRITE(output_unit,'(A,1X,ES12.4)') 'Warning :: |B0| is large (B0 =', B0, '). Check normalization/units.'
   END IF

   IF (10.0_dp < ABS(R0)) THEN
      WRITE(output_unit,'(A,1X,ES12.4)') 'Warning :: |R0| is large (R0 =', R0, '). Check normalization/units.'
   END IF

   IF (a0 > 1.0_dp) THEN
      WRITE(output_unit,'(A,1X,ES12.4,A)') 'Warning :: a0 > 1 (a0 =', a0, &
           '); in normalized units this implies minor radius exceeds major radius.'
   END IF

   IF (ABS(safe2) > 0.0_dp) THEN
      WRITE(output_unit,'(A,1X,ES12.4)') 'Warning :: safe2 is not 0 (safe2 =', safe2, ').'
   END IF

   IF (Phi > 0.2_dp) THEN
      WRITE(output_unit,'(A,1X,ES12.4)') 'Warning :: turbulence amplitude Phi is high (Phi =', Phi, ').'
   END IF

END SUBROUTINE error_signal

SUBROUTINE printall
   USE constants

write(*,*) 't0', t0
write(*,*) 'tc', tc
write(*,*) 'tt', tt
write(*,*) 'tmax', tmax
write(*,*) 'Np', Np
write(*,*) 'Nc', Nc
write(*,*) 'Nloop', Nloop
write(*,*) 'Nreal', Nreal
write(*,*) 'Nt', Nt
write(*,*) 'ntraj', ntraj
write(*,*) 'Ti', Ti
write(*,*) 'Te', Te
write(*,*) 'Zeff', Zeff
write(*,*) 'Aeff', Aeff
write(*,*) 'ndens', ndens
write(*,*) 'Ln', Ln
write(*,*) 'Li', Li
write(*,*) 'Le', Le
write(*,*) 'magnetic_model', magnetic_model
write(*,*) 'B0', B0
write(*,*) 'R0', R0
write(*,*) 'a0', a0
write(*,*) 'safe1', safe1
write(*,*) 'safe2', safe2
write(*,*) 'safe3', safe3
write(*,*) 'amp', amp
write(*,*) 'elong', elong
write(*,*) 'device', device
write(*,*) 'shot', shot
write(*,*) 'shotslice', shotslice
write(*,*) 'Omgt0', Omgt0
write(*,*) 'Omgp0', Omgp0
write(*,*) 'Omgtprim', Omgtprim
write(*,*) 'Phi', Phi
write(*,*) 'turbprof', turbprof
write(*,*) 'Ai', Ai
write(*,*) 'lambdax', lambdax
write(*,*) 'lambday', lambday
write(*,*) 'lambdaz', lambdaz
write(*,*) 'lbalonz', lbalonz
write(*,*) 'tauc', tauc
write(*,*) 'k0i', k0i
write(*,*) 'k0e', k0e
write(*,*) 'gamma_ZF', gamma_ZF
write(*,*) 'X0', X0
write(*,*) 'Y0', Y0
write(*,*) 'Z0', Z0
write(*,*) 'Tw', Tw
write(*,*) 'Ew', Ew
write(*,*) 'pitch', pitch
write(*,*) 'Aw', Aw
write(*,*) 'Zw', Zw
write(*,*) 'position_type', position_type
write(*,*) 'pitch_type', pitch_type
write(*,*) 'energy_type', energy_type
write(*,*) 'USE_larmor', USE_larmor
write(*,*) 'USE_coll', USE_coll
write(*,*) 'USE_turb', USE_turb
write(*,*) 'USE_magnturb', USE_magnturb
write(*,*) 'USE_freq', USE_freq
write(*,*) 'USE_polar', USE_polar
write(*,*) 'USE_PC', USE_PC
write(*,*) 'USE_real', USE_real
write(*,*) 'USE_corr', USE_corr
write(*,*) 'USE_balloon', USE_balloon
write(*,*) 'USE_tilt', USE_tilt
write(*,*) 'USE_testing', USE_testing
pause
stop


END SUBROUTINE printall

