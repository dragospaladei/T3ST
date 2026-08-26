MODULE testing_T3ST_mod
   USE constants
   USE random_numbers              ! provides kx,ky,kz,w,ph,L and kxs,kys,kzs,ws,phs,Ls and USE_real
   USE drift_kernels
   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_nan
   USE, INTRINSIC :: iso_fortran_env, ONLY: OUTPUT_UNIT, ERROR_UNIT
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: testing_T3ST
   PUBLIC :: diagnose_saved_trajectories

   ! -------------------------
   ! Unified formatting knobs
   ! -------------------------
   INTEGER, PARAMETER :: WLBL = 52     ! label column width (bigger to avoid truncation)
   INTEGER, PARAMETER :: WTAG = 7      ! width for tags like "[ OK ]"
   INTEGER, PARAMETER :: WSEC = 70     ! banner width
   INTEGER, PARAMETER :: MAX_ALERTS = 256
   INTEGER, PARAMETER :: MSG_LEN = 220

   ! -------------------------
   ! Logging (tee: terminal + file)
   ! -------------------------
   INTEGER, SAVE :: LOGU = -1
   LOGICAL, SAVE :: LOG_TO_FILE = .FALSE.
   CHARACTER(len=512), SAVE :: LOGFILE = "testing_T3ST.log"

   ! Shared scratch buffer for formatted lines (available to all contained procs)
   CHARACTER(len=1024), SAVE :: TMP

   ! More legible scientific formats
   CHARACTER(len=*), PARAMETER :: FVAL = 'es10.3'  ! values
   CHARACTER(len=*), PARAMETER :: FRAT = 'es10.3'  ! ratios/tolerances
   CHARACTER(len=*), PARAMETER :: FINT = 'i0'      ! integers

   ! -------------------------
   ! Alert collection (so we can recap at the end)
   ! -------------------------
   INTEGER :: n_alerts = 0
   CHARACTER(len=MSG_LEN) :: alert_msgs(MAX_ALERTS)

CONTAINS
   SUBROUTINE diagnose_saved_trajectories(F0, FMtraj, betaFtraj, Htraj, Pctraj, &
                                          Obs_FF, Obs_FD, Obs_DF, Obs_DD, no_errors, fail_counts)
      REAL(wp), INTENT(in) :: F0(:)
      REAL(wp), INTENT(in) :: FMtraj(:, :), betaFtraj(:, :)
      REAL(wp), INTENT(in) :: Htraj(:, :), Pctraj(:, :)
      REAL(wp), INTENT(in) :: Obs_FF(:, :), Obs_FD(:, :), Obs_DF(:, :), Obs_DD(:, :)
      INTEGER, INTENT(out) :: no_errors
      INTEGER, INTENT(out), OPTIONAL :: fail_counts(:)

      REAL(wp), PARAMETER :: tol_identity = 1.0E-3_WP
      REAL(wp), PARAMETER :: tol_obs_abs = 1.0E-6_WP
      REAL(wp), PARAMETER :: tol_obs_rel = 1.0E-5_WP
      REAL(wp), PARAMETER :: tol_H1 = 1.0E-5_WP
      REAL(wp), PARAMETER :: tol_H2 = 2.0E-2_WP
      REAL(wp), PARAMETER :: tol_Pc1 = 1.0E-4_WP
      REAL(wp), PARAMETER :: eps_ref = 1.0E-12_WP

      INTEGER :: ntime, nsaved
      INTEGER :: i, k
      INTEGER :: before_errors
      REAL(wp) :: max_identity_err
      REAL(wp) :: VALUE

      no_errors = 0
      IF (present(fail_counts)) fail_counts = 0
      ntime = size(FMtraj, 1)
      nsaved = size(FMtraj, 2)

      IF (ntime < 1 .OR. nsaved < 1) THEN
         no_errors = no_errors + 1
         IF (present(fail_counts)) fail_counts(1) = fail_counts(1) + 1
         RETURN
      END IF

      IF (size(betaFtraj, 1) /= ntime .OR. size(betaFtraj, 2) /= nsaved .OR. size(F0) < nsaved) THEN
         no_errors = no_errors + 1
         IF (present(fail_counts)) fail_counts(1) = fail_counts(1) + 1
         RETURN
      END IF

      max_identity_err = 0.0_WP
      DO i = 1, nsaved
         IF (F0(i) <= 0.0_WP) THEN
            no_errors = no_errors + 1
            IF (present(fail_counts)) fail_counts(2) = fail_counts(2) + 1
            RETURN
         END IF

         DO k = 1, ntime
            VALUE = FMtraj(k, i)*exp(betaFtraj(k, i))/F0(i)
            max_identity_err = max(max_identity_err, abs(VALUE - 1.0_WP))
         END DO
      END DO

      before_errors = no_errors
      CALL check_scalar_tolerance("FM*exp(betaF)/F0 over saved trajectories", max_identity_err, tol_identity, no_errors)
      IF (present(fail_counts)) fail_counts(3) = fail_counts(3) + no_errors - before_errors

      IF (size(Htraj, 2) >= 1) THEN
         before_errors = no_errors
         CALL check_time_constancy("Htraj(:,1)", Htraj(:, 1), tol_H1, eps_ref, no_errors)
         IF (present(fail_counts)) fail_counts(4) = fail_counts(4) + no_errors - before_errors
      ELSE
         no_errors = no_errors + 1
         IF (present(fail_counts)) fail_counts(4) = fail_counts(4) + 1
      END IF

      IF (size(Htraj, 2) >= 2) THEN
         before_errors = no_errors
         CALL check_time_constancy("Htraj(:,2)", Htraj(:, 2), tol_H2, eps_ref, no_errors)
         IF (present(fail_counts)) fail_counts(5) = fail_counts(5) + no_errors - before_errors
      ELSE
         no_errors = no_errors + 1
         IF (present(fail_counts)) fail_counts(5) = fail_counts(5) + 1
      END IF

      IF (size(Pctraj, 2) >= 1) THEN
         before_errors = no_errors
         CALL check_time_constancy("Pctraj(:,1)", Pctraj(:, 1), tol_Pc1, eps_ref, no_errors)
         IF (present(fail_counts)) fail_counts(6) = fail_counts(6) + no_errors - before_errors
      ELSE
         no_errors = no_errors + 1
         IF (present(fail_counts)) fail_counts(6) = fail_counts(6) + 1
      END IF

      before_errors = no_errors
      CALL check_observable_identity("Obs_FF - Obs_FD = Obs_DF - Obs_DD", &
                                     Obs_FF - Obs_FD, Obs_DF - Obs_DD, &
                                     tol_obs_abs, tol_obs_rel, eps_ref, no_errors)
      IF (present(fail_counts)) fail_counts(7) = fail_counts(7) + no_errors - before_errors

      before_errors = no_errors
      CALL check_observable_identity("Obs_FD - Obs_DD = Obs_FF - Obs_DF", &
                                     Obs_FD - Obs_DD, Obs_FF - Obs_DF, &
                                     tol_obs_abs, tol_obs_rel, eps_ref, no_errors)
      IF (present(fail_counts)) fail_counts(8) = fail_counts(8) + no_errors - before_errors
   END SUBROUTINE diagnose_saved_trajectories

   !=============================================================================
   ! Main entry point
   !=============================================================================
   SUBROUTINE testing_T3ST(X, Y, Z, Vp, mu, F0, G0, FM, betaF, betaD, Vstar1, Vstar2, Vstar3, no_errors)
      REAL(wp), INTENT(in) :: X(:), Y(:), Z(:), Vp(:), mu(:)
      REAL(wp), INTENT(in) :: F0(:), G0(:), FM(:), betaF(:), betaD(:)
      REAL(wp), INTENT(in) :: Vstar1, Vstar2, Vstar3
      INTEGER, INTENT(out) :: no_errors

      INTEGER :: Np_partial, Nt_partial, np_avail
      REAL(wp) :: dt, ratio32

      REAL(wp), ALLOCATABLE :: Xh(:, :), Yh(:, :), Zh(:, :), Vph(:, :), muh(:, :)
      REAL(wp), ALLOCATABLE :: Hh(:, :), Pch(:, :), c1h(:, :), c2h(:, :), c3h(:, :), Vtx(:, :), VFx(:, :)

      INTEGER :: nstep_eff
      REAL(wp) :: T, compare
      REAL(wp), PARAMETER :: epsH = 1.0E-300_WP

      no_errors = 0
      ! -------------------------
      ! Controls (hard-coded here)
      ! -------------------------
      Np_partial = 10000
      Nt_partial = 200
      np_avail = size(X)
      Np_partial = min(Np_partial, np_avail)
      dt = REAL(tmax/Nt, wp)   ! base dt used by rk4_propagation below

      ALLOCATE (Xh(1:Nt_partial, 1:Np_partial), &
                Yh(1:Nt_partial, 1:Np_partial), &
                Zh(1:Nt_partial, 1:Np_partial), &
                Vph(1:Nt_partial, 1:Np_partial), &
                muh(1:Nt_partial, 1:Np_partial), &
                Hh(1:Nt_partial, 1:Np_partial), &
                Pch(1:Nt_partial, 1:Np_partial), &
                c1h(1:Nt_partial, 1:Np_partial), &
                c2h(1:Nt_partial, 1:Np_partial), &
                c3h(1:Nt_partial, 1:Np_partial), &
                Vtx(1:Nt_partial, 1:Np_partial), &
                VFx(1:Nt_partial, 1:Np_partial))

      CALL alerts_reset()
      CALL banner_begin()
      CALL log_open("testing_T3ST.log")   ! or any filename you want

      ! -------------------------
      ! Inputs: NaN checks
      ! -------------------------
      CALL section("Inputs: NaN checks")
      CALL check_nan("X", X, no_errors)
      CALL check_nan("Y", Y, no_errors)
      CALL check_nan("Z", Z, no_errors)
      CALL check_nan("Vp", Vp, no_errors)
      CALL check_nan("mu", mu, no_errors)
      CALL check_nan("F0", F0, no_errors)
      CALL check_nan("G0", G0, no_errors)
      CALL check_nan("FM", FM, no_errors)
      CALL check_nan("betaF", betaF, no_errors)
      CALL check_nan("betaD", betaD, no_errors)
      CALL check_nan("Vstar1", Vstar1, no_errors)
      CALL check_nan("Vstar2", Vstar2, no_errors)
      CALL check_nan("Vstar3", Vstar3, no_errors)

      ! -------------------------
      ! NaN checks: waves (choose set based on USE_real)
      ! (Skip L/Ls entirely)
      ! -------------------------
      CALL section("Waves: NaN checks")
      IF (USE_real == ON) THEN
         CALL log_info("Wave-set selected: s (USE_real = ON)")
         CALL check_nan("kxs", kxs, no_errors)
         CALL check_nan("kys", kys, no_errors)
         CALL check_nan("kzs", kzs, no_errors)
         CALL check_nan("ws", ws, no_errors)
         CALL check_nan("phs", phs, no_errors)
      ELSE
         CALL log_info("Wave-set selected: base (USE_real = OFF)")
         CALL check_nan("kx", kx, no_errors)
         CALL check_nan("ky", ky, no_errors)
         CALL check_nan("kz", kz, no_errors)
         CALL check_nan("w", w, no_errors)
         CALL check_nan("ph", ph, no_errors)
      END IF

      ! -------------------------
      ! Vstar checks
      ! -------------------------
      CALL section("Vstar checks")
      CALL log_vals3("[INFO ]", "Vstar1 Vstar2 Vstar3", Vstar1, Vstar2, Vstar3)

      IF (abs(Vstar1) > 1.0E-9_WP) THEN
         no_errors = no_errors + 1
         CALL log_alert("Vstar1 not ~ 0 (radial drift present)")
      ELSE
         CALL log_ok("Vstar1 ~ 0 (radial drift OK)")
      END IF

      IF (abs(Vstar2) <= tiny(1.0_WP)) THEN
         no_errors = no_errors + 1
         CALL log_alert("Vstar2 ~ 0 (cannot evaluate Vstar3/Vstar2)")
      ELSE
         ratio32 = Vstar3/Vstar2
         CALL log_ratio("[INFO ]", "Vstar3/Vstar2", abs(ratio32), 0.1_WP)

         IF (abs(ratio32) > 0.1_WP) THEN
            no_errors = no_errors + 1
            CALL log_alert("Large parallel drift (|Vstar3/Vstar2| > 0.1)")
         ELSE
            CALL log_ok("Parallel drift OK (|Vstar3/Vstar2| <= 0.1)")
         END IF
      END IF

      ! -------------------------
      ! Input consistency checks
      ! -------------------------
      CALL section("Input consistency checks")
      CALL check_inputs(X, Y, Z, Vp, mu, F0, G0, FM, betaF, betaD, a0, 1.0E-12_WP, no_errors)

      ! -------------------------
      ! Wave moment/symmetry checks (same switch)
      ! -------------------------
      CALL section("Wave moment checks")
      IF (USE_real == ON) THEN
         CALL check_waves_any(kxs, kys, kzs, ws, phs, Vstar2, Vstar3, no_errors)
      ELSE
         CALL check_waves_any(kx, ky, kz, w, ph, Vstar2, Vstar3, no_errors)
      END IF

!!!!
!!!let's do something for eulerian correlation
!!!!

      ! =====================================================================
      ! Drift / RK4 propagation
      ! =====================================================================
      CALL section("RK4 propagation")
      CALL log_val("[INFO ]", "Np_partial", REAL(Np_partial, wp))
      CALL log_val("[INFO ]", "Nt_partial", REAL(Nt_partial, wp))
      CALL log_val("[INFO ]", "dt", dt)

      CALL rk4_propagation(X, Y, Z, Vp, mu, betaF, betaD, dt, Nt_partial, Np_partial, &
                           Xh, Yh, Zh, Vph, muh, Hh, Pch, c1h, c2h, c3h, Vtx, VFx)

      nstep_eff = max(1, Nt_partial - 1)
      T = dt*REAL(nstep_eff, wp)

      ! =====================================================================
      ! Energy proxy diagnostics (d ln H / dt)
      ! =====================================================================
      IF (USE_freq /= OFF) THEN
         CALL log_line("[ALERT]", "Non-zero frequencies (USE_freq /= OFF): forget conservation")
      END IF

      CALL section("Energy conservation proxy ")
      CALL print_dlnH_stats(Hh, Np_partial, Nt_partial, T, epsH, no_errors)
      CALL section("Toroidal(canonical) momentul conservation proxy")
      CALL print_dlnPc_stats(Pch, Np_partial, Nt_partial, T, epsH, no_errors)

      ! =====================================================================
      ! Gaussianity / moment checks: whole (time x particles)
      ! =====================================================================
      CALL section("Gaussianity summary: all times x particles (kurt/3 = 1 for Gaussian)")
      compare = 1.0_WP

      CALL print_gaussianity_block("phi", c1h, Np_partial, Nt_partial, no_errors, compare)
      compare = 1.0_WP
      CALL print_gaussianity_block("phix", c2h, Np_partial, Nt_partial, no_errors, compare)
      compare = 1.0_WP
      CALL print_gaussianity_block("phiy", c3h, Np_partial, Nt_partial, no_errors, compare)

      CALL section("Gaussianity drivers (qualitative)")
      CALL log_info("Larmor & Frequencies -> NO effects on Gaussianity of fields")
      CALL log_info("Balloning breaks phix, phiz Gaussianity as I_0[4/Lbalon^2]/I_0[2/Lbalon^2]^2/2")
      CALL log_info("Balloning breaks amplitude as: I_0(1/lbalonz^2)*exp(-lbalonz^2)")
      CALL log_info("Tilting breaks phix/phiz Gaussianit when particles are z distributey; k_x^eff = kx+ky a --> ??")

      ! =====================================================================
      ! Gaussianity / moment checks: time slice (k = 1)
      ! =====================================================================
      CALL section("Gaussianity summary: time slice k=1 (kurt/3 = 1 for Gaussian)")
      compare = 1.0_WP
      IF (USE_balloon .EQ. ON) THEN
         compare = compare**exp(-1.0_WP/lbalonz**2)*(1 + 1.0_WP/lbalonz**2/4.0_WP + 1.0_WP/lbalonz**4/64.0_WP)
      END IF
      CALL print_gaussianity_slice("phi", c1h, 1, Np_partial, no_errors, compare)
      compare = 1.0_WP/lambdax**4 &
                + (pi**4/5.0_WP)*(3.0_WP/lambday**2 + k0i**2)**2 &
                *(2.0_WP*C2/(C1*C3))**4*9.0_WP*(r00/a0**2)**4 &
                + (2.0_WP*pi**2/3.0_WP)*(1.0_WP/lambdax**2) &
                *(3.0_WP/lambday**2 + k0i**2) &
                *(2.0_WP*C2/(C1*C3))**2*3.0_WP*(r00/a0**2)**2
      CALL print_gaussianity_slice("phix", c2h, 1, Np_partial, no_errors, compare)

      CALL print_gaussianity_slice("phiy", c3h, 1, Np_partial, no_errors, compare)

      IF (USE_real .EQ. OFF) THEN
         compare = ((R0/rhoi*Phi)**2)*sum(ky**2)/REAL(Np*Nc, wp)
      ELSEIF (USE_real .EQ. ON) THEN
         compare = ((R0/rhoi*Phi)**2)*sum(ky**2)/REAL(Np*Nc, wp)
      END IF
      IF (USE_balloon .EQ. ON) THEN
         compare = compare**exp(-1.0_WP/lbalonz**2)*(1 + 1.0_WP/lbalonz**2/4.0_WP + 1.0_WP/lbalonz**4/64.0_WP)
      END IF
      ! Larmor

      CALL log_raw(' --------------------------------------------------------------')
      CALL log_raw(' ------ A ballistic/QL estimation of turbulent transport ------')
      WRITE (TMP, '(A,'//FVAL//',A,'//FVAL//',A,'//FVAL//')') &
         '  E[Vtx] =', sum(Vtx(1, :))/Np_partial, ' E[Vtx^2] =', sum(Vtx(1, :)**2)/Np_partial, ' compared to ', compare
      CALL log_raw(TMP)
      CALL log_raw(' --------------------------------------------------------------')
      CALL log_raw(' ------ A ballistic/QL estimation of neoclassical transport ------')
      WRITE (TMP, '(A,'//FVAL//',A,'//FVAL//',A,'//FVAL//')') &
         '  E[Vnx] =', sum(VFx(1, :) - Vtx(1, :))/Np_partial, &
         ' E[Vnx^2] =', sum((VFx(1, :) - Vtx(1, :))**2)/Np_partial, &
         ' compared to (only if theta=pi/2)', compare
      CALL log_raw(TMP)

      ! -------------------------
      ! Final recap + banner end
      ! -------------------------
      CALL alerts_recap()
      CALL log_close()
      CALL banner_end(no_errors)

   END SUBROUTINE testing_T3ST

   !=============================================================================
   ! RK4 propagation (unchanged numerics; formatting/indent only)
   !=============================================================================
   SUBROUTINE rk4_propagation(X, Y, Z, Vp, mu, betaF, betaD, dt, Nt_partial, Np_partial, &
                              Xh, Yh, Zh, Vph, muh, Hh, Pch, c1h, c2h, c3h, Vtxall, VFxall)

      REAL(wp), INTENT(in) :: X(:), Y(:), Z(:), Vp(:), mu(:), betaF(:), betaD(:)
      INTEGER, INTENT(in) :: Np_partial, Nt_partial
      REAL(wp), INTENT(in) :: dt

      REAL(wp), INTENT(out) :: Xh(:, :), Yh(:, :), Zh(:, :), Vph(:, :), muh(:, :)
      REAL(wp), INTENT(out) :: Hh(:, :), Pch(:, :), c1h(:, :), c2h(:, :), c3h(:, :), Vtxall(:, :), VFxall(:, :)

      REAL(wp), POINTER, CONTIGUOUS :: Qx(:), Qy(:), Qz(:), Qw(:), Qph(:)

      REAL(wp) :: xi, yi, zi, mui, vpi, betaFi, betaDi
      REAL(wp) :: vx, vy, vz, vm, ap, vbF, vbD, kin
      REAL(wp) :: Wx, Wy, Wz, Wm, Wvp, WbF, WbD
      REAL(wp) :: q1, q2, q3
      REAL(wp) :: Hi, FMi, B, Vtx, VFx, Pc
      REAL(wp) :: check_1, check_2, check_3, vs_1, vs_2, vs_3, time
      INTEGER :: ias, k

      DO ias = 1, Np_partial

         ! Select per-particle or single-spectrum arrays
         IF (USE_real == OFF) THEN
            Qx => kx(:, ias)
            Qy => ky(:, ias)
            Qz => kz(:, ias)
            Qw => w(:, ias)
            Qph => ph(:, ias)
         ELSE
            Qx => kxs
            Qy => kys
            Qz => kzs
            Qw => ws
            Qph => phs
         END IF

         xi = X(ias)
         yi = Y(ias)
         zi = Z(ias)
         mui = mu(ias)
         vpi = Vp(ias)
         betaFi = betaF(ias)
         betaDi = betaD(ias)

         DO k = 1, Nt_partial

            time = dt*REAL(k - 1, wp)

            CALL gyrocenter_drifts(xi, yi, zi, vpi, mui, q1, q2, q3, time, &
                                   vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx, &
                                   check_1, check_2, check_3, vs_1, vs_2, vs_3, Qx, Qy, Qz, Qw, Qph, Sol_data)

            Xh(k, ias) = xi
            Yh(k, ias) = yi
            Zh(k, ias) = zi
            Vph(k, ias) = vpi
            muh(k, ias) = mui
            Hh(k, ias) = Hi
            Pch(k, ias) = Pc
            c1h(k, ias) = check_1
            c2h(k, ias) = check_2
            c3h(k, ias) = check_3
            Vtxall(k, ias) = Vtx
            VFxall(k, ias) = VFx

            Wx = vx/6.0_WP
            Wy = vy/6.0_WP
            Wz = vz/6.0_WP
            Wm = vm/6.0_WP
            WbF = vbF/6.0_WP
            WbD = vbD/6.0_WP
            Wvp = ap/6.0_WP

            CALL gyrocenter_drifts(xi + vx*dt/2.0_WP, yi + vy*dt/2.0_WP, zi + vz*dt/2.0_WP, &
                                   vpi + ap*dt/2.0_WP, mui + vm*dt/2.0_WP, q1, q2, q3, time + dt/2.0_WP, &
                                   vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx, &
                                   check_1, check_2, check_3, vs_1, vs_2, vs_3, Qx, Qy, Qz, Qw, Qph, Sol_data)

            Wx = Wx + vx/3.0_WP
            Wy = Wy + vy/3.0_WP
            Wz = Wz + vz/3.0_WP
            Wm = Wm + vm/3.0_WP
            WbF = WbF + vbF/3.0_WP
            WbD = WbD + vbD/3.0_WP
            Wvp = Wvp + ap/3.0_WP

! RK2b
            CALL gyrocenter_drifts(xi + vx*dt/2.0_WP, yi + vy*dt/2.0_WP, zi + vz*dt/2.0_WP, &
                                   vpi + ap*dt/2.0_WP, mui + vm*dt/2.0_WP, q1, q2, q3, time + dt/2.0_WP, &
                                   vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx, &
                                   check_1, check_2, check_3, vs_1, vs_2, vs_3, Qx, Qy, Qz, Qw, Qph, Sol_data)

            Wx = Wx + vx/3.0_WP
            Wy = Wy + vy/3.0_WP
            Wz = Wz + vz/3.0_WP
            Wm = Wm + vm/3.0_WP
            WbF = WbF + vbF/3.0_WP
            WbD = WbD + vbD/3.0_WP
            Wvp = Wvp + ap/3.0_WP

! RK4
            CALL gyrocenter_drifts(xi + vx*dt, yi + vy*dt, zi + vz*dt, &
                                   vpi + ap*dt, mui + vm*dt, q1, q2, q3, time + dt, &
                                   vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx, &
                                   check_1, check_2, check_3, vs_1, vs_2, vs_3, Qx, Qy, Qz, Qw, Qph, Sol_data)

            Wx = Wx + vx/6.0_WP
            Wy = Wy + vy/6.0_WP
            Wz = Wz + vz/6.0_WP
            Wm = Wm + vm/6.0_WP
            WbF = WbF + vbF/6.0_WP
            WbD = WbD + vbD/6.0_WP
            Wvp = Wvp + ap/6.0_WP

            xi = xi + dt*Wx
            yi = yi + dt*Wy
            zi = zi + dt*Wz
            vpi = vpi + dt*Wvp
            mui = mui + dt*Wm
            betaFi = betaFi + dt*WbF
            betaDi = betaDi + dt*WbD

         END DO
      END DO
   END SUBROUTINE rk4_propagation

   !=============================================================================
   ! Unified NaN checker (assumed-rank)
   !=============================================================================
   SUBROUTINE check_nan(name, x, no_errors)
      CHARACTER(len=*), INTENT(in) :: name
      REAL(wp), INTENT(in) :: x(..)
      INTEGER, INTENT(inout) :: no_errors

      LOGICAL :: has_nan

      has_nan = .FALSE.
      SELECT rank (x)
      rank (0)
         has_nan = ieee_is_nan(x)
      rank (1)
         has_nan = any(ieee_is_nan(x))
      rank (2)
         has_nan = any(ieee_is_nan(x))
      rank (3)
         has_nan = any(ieee_is_nan(x))
      rank (4)
         has_nan = any(ieee_is_nan(x))
      rank default
         CALL log_warn("check_nan: unsupported rank for "//trim(name))
         RETURN
      END SELECT

      IF (has_nan) THEN
         no_errors = no_errors + 1
         CALL log_alert(trim(name)//" contains NaN")
      ELSE
         CALL log_ok("No NaNs in "//trim(name))
      END IF
   END SUBROUTINE check_nan

   !=============================================================================
   ! Section header
   !=============================================================================
   SUBROUTINE section(title)
      CHARACTER(len=*), INTENT(in) :: title
      CALL log_raw(" ")
      CALL log_line("[INFO ]", repeat("-", 1)//" "//trim(title)//" "//repeat("-", max(0, WSEC - 9 - len_trim(title))))
   END SUBROUTINE section

   !=============================================================================
   ! Banners
   !=============================================================================
   SUBROUTINE banner_begin()
      CALL log_raw(" ")
      CALL log_raw(repeat("=", WSEC))
      CALL log_raw(" TESTING / DIAGNOSTICS")
      CALL log_raw(repeat("=", WSEC))
   END SUBROUTINE banner_begin

   SUBROUTINE banner_end(no_errors)
      INTEGER, INTENT(in) :: no_errors
      CHARACTER(len=256) :: tmp2
      CALL log_raw(repeat("-", WSEC))
      WRITE (tmp2, '(a,1x,'//FINT//')') "[INFO ] Total number of alerts:", no_errors
      CALL log_raw(tmp2)
      CALL log_raw(repeat("=", WSEC))
   END SUBROUTINE banner_end

   !=============================================================================
   ! Input checks (uniform messages + formatted numbers)
   !=============================================================================
   SUBROUTINE check_inputs(X, Y, Z, Vp, mu, F0, G0, FM, betaF, betaD, a1, tol, no_errors)
      REAL(wp), INTENT(in) :: X(:), Y(:), Z(:), Vp(:), mu(:)
      REAL(wp), INTENT(in) :: F0(:), G0(:), FM(:), betaF(:), betaD(:)
      REAL(wp), INTENT(in) :: a1, tol
      INTEGER, INTENT(inout) :: no_errors

      INTEGER :: n, nbad
      REAL(wp) :: pii, max_betaF_err, max_betaD_err

      n = size(X)
      pii = 3.14159265358979323846_WP

      IF (any([size(Y), size(Z), size(Vp), size(mu), size(F0), size(G0), size(FM), size(betaF), size(betaD)] /= n)) THEN
         no_errors = no_errors + 1
         CALL log_alert("Input arrays have inconsistent sizes")
         WRITE (TMP, '(a,1x,a,": X=",1x,'//FINT//',", Y=",1x,'//FINT//',", Z=",1x,'//FINT//', &
   &                   ", Vp=",1x,'//FINT//',", mu=",1x,'//FINT//',", F0=",1x,'//FINT//', &
   &                   ", G0=",1x,'//FINT//',", FM=",1x,'//FINT//',", betaF=",1x,'//FINT//', &
   &                   ", betaD=",1x,'//FINT//')') &
           "[INFO ]", pad_label("sizes"), size(X), size(Y), size(Z), size(Vp), size(mu), &
           size(F0), size(G0), size(FM), size(betaF), size(betaD)
         CALL log_raw(TMP)
         RETURN
      ELSE
         CALL log_ok("All input arrays have consistent size")
         WRITE (TMP, '(a,1x,a,":",1x,'//FINT//')') "[INFO ]", pad_label("n"), n
         CALL log_raw(TMP)
      END IF

      nbad = count(X <= 0.0_WP)
      IF (nbad > 0) THEN
         no_errors = no_errors + 1
         CALL log_alert("Condition failed: X must be > 0")
         CALL log_count("violations", nbad, n)
      ELSE
         CALL log_ok("X > 0 everywhere")
      END IF

      nbad = count(abs(Y) >= 1.0_WP)
      IF (nbad > 0) THEN
         no_errors = no_errors + 1
         CALL log_alert("Condition failed: abs(Y) must be < 1")
         CALL log_count("violations", nbad, n)
      ELSE
         CALL log_ok("abs(Y) < 1 everywhere")
      END IF

      nbad = count((Z <= -pii) .OR. (Z >= pii))
      IF (nbad > 0) THEN
         no_errors = no_errors + 1
         CALL log_alert("Condition failed: Z must be in (-pi, pi)")
         CALL log_count("violations", nbad, n)
      ELSE
         CALL log_ok("Z in (-pi, pi) everywhere")
      END IF

      nbad = count(mu <= 0.0_WP)
      IF (nbad > 0) THEN
         no_errors = no_errors + 1
         CALL log_alert("Condition failed: mu must be > 0")
         CALL log_count("violations", nbad, n)
      ELSE
         CALL log_ok("mu > 0 everywhere")
      END IF

      nbad = count(F0 <= 0.0_WP)
      IF (nbad > 0) THEN
         no_errors = no_errors + 1
         CALL log_alert("Condition failed: F0 must be > 0")
         CALL log_count("violations", nbad, n)
      ELSE
         CALL log_ok("F0 > 0 everywhere")
      END IF

      nbad = count(G0 <= 0.0_WP)
      IF (nbad > 0) THEN
         no_errors = no_errors + 1
         CALL log_alert("Condition failed: G0 must be > 0")
         CALL log_count("violations", nbad, n)
      ELSE
         CALL log_ok("G0 > 0 everywhere")
      END IF

      nbad = count(FM <= 0.0_WP)
      IF (nbad > 0) THEN
         no_errors = no_errors + 1
         CALL log_alert("Condition failed: FM must be > 0")
         CALL log_count("violations", nbad, n)
      ELSE
         CALL log_ok("FM > 0 everywhere")
      END IF

      IF (all((F0 > 0.0_WP) .AND. (FM > 0.0_WP))) THEN
         max_betaF_err = maxval(abs(FM*exp(betaF)/F0 - 1.0_WP))
         max_betaD_err = maxval(abs(FM*exp(betaD)/F0 - 1.0_WP))

         IF (max_betaF_err > tol) THEN
            no_errors = no_errors + 1
            CALL log_alert("Condition failed: FM*exp(betaF)/F0 must be 1")
         ELSE
            CALL log_ok("FM*exp(betaF)/F0 = 1 within tolerance")
         END IF
         CALL log_val("[INFO ]", "max betaF identity error", max_betaF_err)

         IF (max_betaD_err > tol) THEN
            no_errors = no_errors + 1
            CALL log_alert("Condition failed: FM*exp(betaD)/F0 must be 1")
         ELSE
            CALL log_ok("FM*exp(betaD)/F0 = 1 within tolerance")
         END IF
         CALL log_val("[INFO ]", "max betaD identity error", max_betaD_err)
         CALL log_val("[INFO ]", "tol", tol)
      END IF

      nbad = count(sqrt((X - 1.0_WP)**2 + Y**2) >= a1)
      IF (nbad > 0) THEN
         no_errors = no_errors + 1
         CALL log_alert("Condition failed: sqrt((X-1)^2 + Y^2) < a0")
         CALL log_count("violations", nbad, n)
         CALL log_val("[INFO ]", "a0", a1)
      ELSE
         CALL log_ok("sqrt((X-1)^2 + Y^2) < a0 everywhere")
         CALL log_val("[INFO ]", "a0", a1)
      END IF

   END SUBROUTINE check_inputs

   !=============================================================================
   ! Wave checks (rank-agnostic) WITHOUT L/Ls
   !=============================================================================
   SUBROUTINE check_waves_any(kxA, kyA, kzA, wA, phA, Vstar2, Vstar3, no_errors)
      REAL(wp), INTENT(in) :: kxA(..), kyA(..), kzA(..), wA(..), phA(..)
      REAL(wp), INTENT(in) :: Vstar2, Vstar3
      INTEGER, INTENT(inout) :: no_errors

      REAL(wp), PARAMETER :: tiny_wp = 1.0E-12_WP
      REAL(wp), PARAMETER :: eps_ref = 1.0E-30_WP

      REAL(wp) :: m1(5), m2(5), m4(5)
      REAL(wp) :: tol_mean(5), tol_m2(5), tol_m4(5)
      REAL(wp) :: m2_ref(5), m4_ref(5)
      REAL(wp) :: ratio, pii
      REAL(wp) :: t_m2_ky, t_m2_kz, t_m4_ky, t_m4_kz
      REAL(wp) :: fac

      CALL log_info("Wave diagnostics (moments / target checks)")

      pii = 3.14159265358979323846_WP

      tol_mean = 0.05_WP
      tol_m2 = 0.05_WP
      tol_m4 = 0.10_WP

      m2_ref(2) = 3.0_WP/(lambday*lambday) + k0i*k0i
      m2_ref(3) = (erfc(lambdaz/sqrt(2.0_WP)) &
                   + 3.0_WP*erfc(lambdaz*sqrt(2.0_WP)))/C3**2
      m2_ref(5) = (pii*pii)/3.0_WP

      ! Radial spectrum selection
      IF (x_corr == 1) THEN
         m2_ref(1) = 1.0_WP/(lambdax**2)
         m4_ref(1) = 3.0_WP/(lambdax**4)
      ELSEIF (x_corr == 2) THEN
         m2_ref(1) = (10.0_WP/atan(10.0_WP) - 1.0_WP)/(lambdax**2)
         m4_ref(1) = (1.0_WP + 10.0_WP*(10.0_WP**2 - 3)/3.0_WP/atan(10.0_WP))/(lambdax**4)
      END IF

      ! Temporal spectrum selection
      IF (t_corr == 1) THEN
         m2_ref(4) = 1.0_WP/(tauc**2)
         m4_ref(4) = 3.0_WP/(tauc**4)
      ELSEIF (t_corr == 2) THEN
         m2_ref(4) = (10.0_WP/atan(10.0_WP) - 1.0_WP)/(tauc**2)
         m4_ref(4) = (1.0_WP + 10.0_WP*(10.0_WP**2 - 3)/3.0_WP/atan(10.0_WP))/(tauc**4)
      END IF

      m4_ref(2) = (15.0_WP/(lambday**4) + 10.0_WP*k0i**2/(lambday**2) + k0i**4)
     m4_ref(3) = (erfc(lambdaz/sqrt(2.0_WP)) + 15.0_WP*erfc(lambdaz*sqrt(2.0_WP)) + 65.0_WP*erfc(lambdaz*3.0_WP/sqrt(2.0_WP)))/C3**4
      m4_ref(5) = (pii**4)/5.0_WP

      t_m2_ky = m2_ref(2)
      t_m2_kz = m2_ref(3)
      t_m4_ky = m4_ref(2)
      t_m4_kz = m4_ref(3)

      fac = (Ln*(rhoi/R0/(1.0_WP + r00))**2)

      m2_ref(4) = m2_ref(4) + (fac**2)*(Vstar2**2*t_m2_ky + Vstar3**2*t_m2_kz)

      m4_ref(4) = m4_ref(4) &
                  + (fac**4)*(Vstar2**4*t_m4_ky + Vstar3**4*t_m4_kz &
                              + 2.0_WP*Vstar2**2*Vstar3**2*t_m2_kz*t_m2_ky) &
                  + 2.0_WP/(tauc**2)*(fac**2)*(Vstar2**2*t_m2_ky + Vstar3**2*t_m2_kz)

      CALL moments_ar(kxA, m1(1), m2(1), m4(1))
      CALL moments_ar(kyA, m1(2), m2(2), m4(2))
      CALL moments_ar(kzA, m1(3), m2(3), m4(3))
      CALL moments_ar(wA, m1(4), m2(4), m4(4))
      CALL moments_ar(phA, m1(5), m2(5), m4(5))

      ratio = abs(m1(1))/max(sqrt(m2(1)), tiny_wp)
      CALL check_ratio("kx mean ~ 0", ratio, tol_mean(1), no_errors)

      ratio = abs(m1(2))/max(sqrt(m2(2)), tiny_wp)
      CALL check_ratio("ky mean ~ 0", ratio, tol_mean(2), no_errors)

      ratio = abs(m1(3))/max(sqrt(m2(3)), tiny_wp)
      CALL check_ratio("kz mean ~ 0", ratio, tol_mean(3), no_errors)

      ratio = abs(m1(4))/max(sqrt(m2(4)), tiny_wp)
      CALL check_ratio("w  mean ~ 0", ratio, tol_mean(4), no_errors)

      ratio = abs(m1(5))/max(sqrt(m2(5)), tiny_wp)
      CALL check_ratio("ph mean ~ 0", ratio, tol_mean(5), no_errors)

      CALL log_info("Second moments (targets)")
      CALL check_moment("kx^2 target", m2(1), m2_ref(1), tol_m2(1), no_errors, eps_ref)
      CALL check_moment("ky^2 target", m2(2), m2_ref(2), tol_m2(2), no_errors, eps_ref)
      CALL check_moment("kz^2 target", m2(3), m2_ref(3), tol_m2(3), no_errors, eps_ref)
      CALL check_moment("w^2  target", m2(4), m2_ref(4), tol_m2(4), no_errors, eps_ref)
      CALL check_moment("ph^2 target", m2(5), m2_ref(5), tol_m2(5), no_errors, eps_ref)

      CALL log_info("Fourth moments (targets)")
      CALL check_moment("kx^4 target", m4(1), m4_ref(1), tol_m4(1), no_errors, eps_ref)
      CALL check_moment("ky^4 target", m4(2), m4_ref(2), tol_m4(2), no_errors, eps_ref)
      CALL check_moment("kz^4 target", m4(3), m4_ref(3), tol_m4(3), no_errors, eps_ref)
      CALL check_moment("w^4  target", m4(4), m4_ref(4), tol_m4(4), no_errors, eps_ref)
      CALL check_moment("ph^4 target", m4(5), m4_ref(5), tol_m4(5), no_errors, eps_ref)

   END SUBROUTINE check_waves_any

   SUBROUTINE check_moment(label, val, ref, tol, no_errors, eps_ref)
      CHARACTER(len=*), INTENT(in) :: label
      REAL(wp), INTENT(in) :: val, ref, tol, eps_ref
      INTEGER, INTENT(inout) :: no_errors
      REAL(wp) :: relerr

      relerr = abs(val - ref)/max(abs(ref), eps_ref)

      IF (relerr > tol) THEN
         no_errors = no_errors + 1
         CALL log_moment("[ALERT]", label, val, ref, relerr, tol)
      ELSE
         CALL log_moment("[  OK ]", label, val, ref, relerr, tol)
      END IF
   END SUBROUTINE check_moment

   SUBROUTINE check_scalar_tolerance(label, err, tol, no_errors)
      CHARACTER(len=*), INTENT(in) :: label
      REAL(wp), INTENT(in) :: err, tol
      INTEGER, INTENT(inout) :: no_errors

      IF (ieee_is_nan(err) .OR. err > tol) THEN
         no_errors = no_errors + 1
      END IF
   END SUBROUTINE check_scalar_tolerance

   SUBROUTINE check_time_constancy(label, x, tol, eps_ref, no_errors)
      CHARACTER(len=*), INTENT(in) :: label
      REAL(wp), INTENT(in) :: x(:)
      REAL(wp), INTENT(in) :: tol, eps_ref
      INTEGER, INTENT(inout) :: no_errors

      REAL(wp) :: err
      REAL(wp) :: ref

      IF (size(x) < 1) THEN
         no_errors = no_errors + 1
         RETURN
      END IF

      IF (any(ieee_is_nan(x))) THEN
         no_errors = no_errors + 1
         RETURN
      END IF

      ref = x(1)
      err = maxval(abs(x - ref))/max(abs(ref), eps_ref)
      CALL check_scalar_tolerance(trim(label)//" time constancy", err, tol, no_errors)
   END SUBROUTINE check_time_constancy

   SUBROUTINE check_observable_identity(label, lhs, rhs, abs_tol, rel_tol, eps_ref, no_errors)
      CHARACTER(len=*), INTENT(in) :: label
      REAL(wp), INTENT(in) :: lhs(:, :), rhs(:, :)
      REAL(wp), INTENT(in) :: abs_tol, rel_tol, eps_ref
      INTEGER, INTENT(inout) :: no_errors

      REAL(wp) :: err, scale

      IF (size(lhs, 1) /= size(rhs, 1) .OR. size(lhs, 2) /= size(rhs, 2)) THEN
         no_errors = no_errors + 1
         RETURN
      END IF

      IF (any(ieee_is_nan(lhs)) .OR. any(ieee_is_nan(rhs))) THEN
         no_errors = no_errors + 1
         RETURN
      END IF

      err = maxval(abs(lhs - rhs))
      scale = max(maxval(abs(lhs)), maxval(abs(rhs)), eps_ref)

      IF (err > abs_tol + rel_tol*scale) THEN
         no_errors = no_errors + 1
      END IF

      CALL ignore_unused(REAL(len_trim(label), wp))
   END SUBROUTINE check_observable_identity

   SUBROUTINE log_moment(tag, label, val, ref, relerr, tol)
      CHARACTER(len=*), INTENT(in) :: tag, label
      REAL(wp), INTENT(in) :: val, ref, relerr, tol

      WRITE (TMP, '(a,1x,a,": val=",1x,'//FVAL//',", ref=",1x,'//FVAL//',", rel=",1x,'//FRAT//',", tol=",1x,'//FRAT//')') &
         tag, pad_label(label), val, ref, relerr, tol
      CALL log_raw(TMP)

      IF (trim(tag) == "[ALERT]") CALL alerts_add(trim(label))
   END SUBROUTINE log_moment

   !=============================================================================
   ! Moments for assumed-rank arrays (rank 0/1/2 supported)
   !=============================================================================
   SUBROUTINE moments_ar(x, m1, m2, m4)
      REAL(wp), INTENT(in) :: x(..)
      REAL(wp), INTENT(out) :: m1, m2, m4

      INTEGER :: n
      REAL(wp) :: rn

      SELECT rank (x)
      rank (0)
         n = 1
         rn = 1.0_WP
         m1 = x
         m2 = x*x
         m4 = (x*x)*(x*x)

      rank (1)
         n = size(x)
         rn = REAL(n, wp)
         m1 = sum(x)/rn
         m2 = sum(x**2)/rn
         m4 = sum(x**4)/rn

      rank (2)
         n = size(x)
         rn = REAL(n, wp)
         m1 = sum(x)/rn
         m2 = sum(x**2)/rn
         m4 = sum(x**4)/rn

      rank default
         CALL log_warn("moments_ar: unsupported rank encountered")
         m1 = huge(1.0_WP)
         m2 = huge(1.0_WP)
         m4 = huge(1.0_WP)
      END SELECT
   END SUBROUTINE moments_ar

   !=============================================================================
   ! Energy proxy: d(ln H)/dt finite difference stats
   !=============================================================================
   SUBROUTINE print_dlnH_stats(Hh, Np, Nt, T, epsH, no_errors)
      REAL(wp), INTENT(in) :: Hh(:, :)
      INTEGER, INTENT(in) :: Np, Nt
      REAL(wp), INTENT(in) :: T, epsH
      INTEGER, INTENT(inout) :: no_errors

      REAL(wp), ALLOCATABLE :: v(:)
      REAL(wp) :: mean_v, rms_v, std_v, maxabs_v, meanabs_v
      INTEGER :: imax(1)

      ALLOCATE (v(Np))

      v = (Hh(Nt, :) - Hh(1, :))/max(abs(Hh(1, :)), epsH)/T

      mean_v = sum(v)/REAL(Np, wp)
      rms_v = sqrt(sum(v**2)/REAL(Np, wp))
      std_v = sqrt(max(0.0_WP, rms_v**2 - mean_v**2))
      maxabs_v = maxval(abs(v))
      meanabs_v = sum(abs(v))/REAL(Np, wp)

      CALL log_val("[INFO ]", "mean (H(tf)/H(ti)-1)/(tf-ti)", mean_v)
      CALL log_val("[INFO ]", "rms  (H(tf)/H(ti)-1)/(tf-ti)", rms_v)
      CALL log_val("[INFO ]", "std  (H(tf)/H(ti)-1)/(tf-ti)", std_v)
      CALL log_val("[INFO ]", "max |(H(tf)/H(ti)-1)/(tf-ti)|", maxabs_v)
      CALL log_val("[INFO ]", "mean|(H(tf)/H(ti)-1)/(tf-ti)|", meanabs_v)

      imax = maxloc(abs(v))
      WRITE (TMP, '(a,1x,a,":",1x,'//FINT//')') "[INFO ]", pad_label("argmax |(H(tf)/H(ti)-1)/(tf-ti)| (particle index)"), imax(1)
      CALL log_raw(TMP)

      IF (maxabs_v > 1.0E-1_WP) THEN
         no_errors = no_errors + 1
         CALL log_alert("Energy proxy: large max |(H(tf)/H(ti)-1)/(tf-ti)| (threshold 1e-1)")
      END IF

      DEALLOCATE (v)
   END SUBROUTINE print_dlnH_stats

   !=============================================================================
   ! Energy proxy: d(ln H)/dt finite difference stats
   !=============================================================================
   SUBROUTINE print_dlnPc_stats(Hh, Np, Nt, T, epsH, no_errors)
      REAL(wp), INTENT(in) :: Hh(:, :)
      INTEGER, INTENT(in) :: Np, Nt
      REAL(wp), INTENT(in) :: T, epsH
      INTEGER, INTENT(inout) :: no_errors

      REAL(wp), ALLOCATABLE :: v(:)
      REAL(wp) :: mean_v, rms_v, std_v, maxabs_v, meanabs_v
      INTEGER :: imax(1)

      ALLOCATE (v(Np))

      v = (Hh(Nt, :) - Hh(1, :))/max(abs(Hh(1, :)), epsH)/T

      mean_v = sum(v)/REAL(Np, wp)
      rms_v = sqrt(sum(v**2)/REAL(Np, wp))
      std_v = sqrt(max(0.0_WP, rms_v**2 - mean_v**2))
      maxabs_v = maxval(abs(v))
      meanabs_v = sum(abs(v))/REAL(Np, wp)

      CALL log_val("[INFO ]", "mean (Pc(tf)/Pc(ti)-1)/(tf-ti)", mean_v)
      CALL log_val("[INFO ]", "rms  (Pc(tf)/Pc(ti)-1)/(tf-ti)", rms_v)
      CALL log_val("[INFO ]", "std  (Pc(tf)/Pc(ti)-1)/(tf-ti)", std_v)
      CALL log_val("[INFO ]", "max |(Pc(tf)/Pc(ti)-1)/(tf-ti)|", maxabs_v)
      CALL log_val("[INFO ]", "mean|(Pc(tf)/Pc(ti)-1)/(tf-ti)|", meanabs_v)

      imax = maxloc(abs(v))
      WRITE (TMP, '(a,1x,a,":",1x,'//FINT//')') "[INFO ]", pad_label("argmax |(Pc(tf)/Pc(ti)-1)/(tf-ti)| (particle index)"), imax(1)
      CALL log_raw(TMP)

      IF (maxabs_v > 1.0E-1_WP) THEN
         no_errors = no_errors + 1
         CALL log_alert("Energy proxy: large max |(Pc(tf)/Pc(ti)-1)/(tf-ti)| (threshold 1e-1)")
      END IF

      DEALLOCATE (v)
   END SUBROUTINE print_dlnPc_stats

   !=============================================================================
   ! Gaussianity summary over full 2D field (Nt x Np)
   !=============================================================================
   SUBROUTINE print_gaussianity_block(name, A, Np, Nt, no_errors, compare)
      CHARACTER(len=*), INTENT(in) :: name
      REAL(wp), INTENT(in) :: A(:, :), compare
      INTEGER, INTENT(in) :: Np, Nt
      INTEGER, INTENT(inout) :: no_errors

      REAL(wp) :: norm, m1, m2, m4, kurt3
      REAL(wp), PARAMETER :: tol_kurt3 = 0.15_WP

      norm = 1.0_WP/(REAL(Np, wp)*REAL(Nt, wp))
      m1 = sum(A)*norm
      m2 = sum(A**2)*norm
      m4 = sum(A**4)*norm
      kurt3 = m4/max(m2*m2, tiny(1.0_WP))/3.0_WP

      CALL log_gauss_line(name, m1, m2, kurt3, tol_kurt3, no_errors, compare)
   END SUBROUTINE print_gaussianity_block

   !=============================================================================
   ! Gaussianity summary for one time slice k (1..Nt): uses A(k,:)
   !=============================================================================
   SUBROUTINE print_gaussianity_slice(name, A, k, Np, no_errors, compare)
      CHARACTER(len=*), INTENT(in) :: name
      REAL(wp), INTENT(in) :: A(:, :), compare
      INTEGER, INTENT(in) :: k, Np
      INTEGER, INTENT(inout) :: no_errors

      REAL(wp) :: norm, m1, m2, m4, kurt3
      REAL(wp), PARAMETER :: tol_kurt3 = 0.15_WP

      norm = 1.0_WP/REAL(Np, wp)
      m1 = sum(A(k, :))*norm
      m2 = sum(A(k, :)**2)*norm
      m4 = sum(A(k, :)**4)*norm
      kurt3 = m4/max(m2*m2, tiny(1.0_WP))/3.0_WP

      CALL log_gauss_line(name, m1, m2, kurt3, tol_kurt3, no_errors, compare)
   END SUBROUTINE print_gaussianity_slice

   SUBROUTINE log_gauss_line(name, mean, ex2, kurt3, tol_kurt3, no_errors, compare)
      CHARACTER(len=*), INTENT(in) :: name
      REAL(wp), INTENT(in) :: mean, ex2, kurt3, tol_kurt3, compare
      INTEGER, INTENT(inout) :: no_errors

      CHARACTER(len=WTAG) :: tag
      LOGICAL :: bad

      bad = abs(kurt3 - 1.0_WP) > tol_kurt3
      IF (bad) THEN
         tag = "[ALERT]"
         no_errors = no_errors + 1
         CALL alerts_add("Gaussianity: "//trim(name)//" kurt/3 out of tol")
      ELSE
         tag = "[INFO ]"
      END IF

      WRITE (TMP, '(a,1x,a,": E[x]=",1x,'//FVAL//',", E[x^2]=",1x,'//FVAL//', &
  &              ", kurt/3=",1x,'//FVAL//',", tol=",1x,'//FVAL//',", compared to = ",1x,'//FVAL//')') &
           tag, pad_label(trim(name)), mean, ex2, kurt3, tol_kurt3, compare
      CALL log_raw(TMP)
   END SUBROUTINE log_gauss_line

   !=============================================================================
   ! Helper: check ratio against tol with uniform messages
   !=============================================================================
   SUBROUTINE check_ratio(label, ratio, tol, no_errors)
      CHARACTER(len=*), INTENT(in) :: label
      REAL(wp), INTENT(in) :: ratio, tol
      INTEGER, INTENT(inout) :: no_errors

      IF (ratio > tol) THEN
         no_errors = no_errors + 1
         CALL log_ratio("[ALERT]", trim(label), ratio, tol)
      ELSE
         CALL log_ratio("[  OK ]", trim(label), ratio, tol)
      END IF
   END SUBROUTINE check_ratio

   !=============================================================================
   ! Logging utilities (centralized formatting)
   !=============================================================================
   SUBROUTINE log_ok(msg)
      CHARACTER(len=*), INTENT(in) :: msg
      CALL log_line("[  OK ]", msg)
   END SUBROUTINE log_ok

   SUBROUTINE log_alert(msg)
      CHARACTER(len=*), INTENT(in) :: msg
      CALL alerts_add(trim(msg))
      CALL log_line("[ALERT]", msg)
   END SUBROUTINE log_alert

   SUBROUTINE log_warn(msg)
      CHARACTER(len=*), INTENT(in) :: msg
      CALL log_line("[WARN]", msg)
   END SUBROUTINE log_warn

   SUBROUTINE log_info(msg)
      CHARACTER(len=*), INTENT(in) :: msg
      CALL log_line("[INFO ]", msg)
   END SUBROUTINE log_info

   SUBROUTINE log_line(tag, msg)
      USE, INTRINSIC :: iso_fortran_env, ONLY: OUTPUT_UNIT
      CHARACTER(len=*), INTENT(in) :: tag, msg

      WRITE (OUTPUT_UNIT, '(a,1x,a)') tag, trim(msg)
      IF (LOG_TO_FILE) WRITE (LOGU, '(a,1x,a)') tag, trim(msg)
   END SUBROUTINE log_line

   SUBROUTINE log_val(tag, name, val)
      CHARACTER(len=*), INTENT(in) :: tag, name
      REAL(wp), INTENT(in) :: val
      WRITE (TMP, '(a,1x,a,":",1x,'//FVAL//')') tag, pad_label(name), val
      CALL log_raw(TMP)
   END SUBROUTINE log_val

   SUBROUTINE log_vals3(tag, name, v1, v2, v3)
      CHARACTER(len=*), INTENT(in) :: tag, name
      REAL(wp), INTENT(in) :: v1, v2, v3
      WRITE (TMP, '(a,1x,a,":",3(1x,'//FVAL//'))') tag, pad_label(name), v1, v2, v3
      CALL log_raw(TMP)
   END SUBROUTINE log_vals3

   SUBROUTINE log_ratio(tag, name, ratio, tol)
      CHARACTER(len=*), INTENT(in) :: tag, name
      REAL(wp), INTENT(in) :: ratio, tol
      WRITE (TMP, '(a,1x,a,": ratio=",1x,'//FRAT//',", tol=",1x,'//FRAT//')') tag, pad_label(name), ratio, tol
      CALL log_raw(TMP)
      IF (trim(tag) == "[ALERT]") CALL alerts_add(trim(name)//" ratio exceeded tol")
   END SUBROUTINE log_ratio

   SUBROUTINE log_count(label, nbad, n)
      CHARACTER(len=*), INTENT(in) :: label
      INTEGER, INTENT(in) :: nbad, n
      WRITE (TMP, '(a,1x,a,":",1x,'//FINT//',a,1x,'//FINT//')') "[INFO ]", pad_label(label), nbad, " /", n
      CALL log_raw(TMP)
   END SUBROUTINE log_count

   !=============================================================================
   ! Label padding
   !=============================================================================
   FUNCTION pad_label(s) RESULT(out)
      CHARACTER(len=*), INTENT(in) :: s
      CHARACTER(len=WLBL) :: out
      INTEGER :: ls

      out = " "
      ls = len_trim(s)

      IF (ls >= WLBL) THEN
         out = s(1:WLBL)
      ELSE
         out(1:ls) = s(1:ls)
      END IF
   END FUNCTION pad_label

   !=============================================================================
   ! Alert collection + recap
   !=============================================================================
   SUBROUTINE alerts_reset()
      INTEGER :: i
      n_alerts = 0
      DO i = 1, MAX_ALERTS
         alert_msgs(i) = " "
      END DO
   END SUBROUTINE alerts_reset

   SUBROUTINE alerts_add(msg)
      CHARACTER(len=*), INTENT(in) :: msg
      IF (n_alerts < MAX_ALERTS) THEN
         n_alerts = n_alerts + 1
         alert_msgs(n_alerts) = trim(msg)
      END IF
   END SUBROUTINE alerts_add

   SUBROUTINE alerts_recap()
      INTEGER :: i
      IF (n_alerts <= 0) RETURN

      CALL log_raw(" ")
      CALL log_info("Alerts recap (most recent run)")
      DO i = 1, n_alerts
         CALL log_line("[ALERT]", trim(alert_msgs(i)))
      END DO
   END SUBROUTINE alerts_recap

   !=============================================================================
   ! Small helper for unused dummy variable
   !=============================================================================
   SUBROUTINE ignore_unused(x)
      REAL(wp), INTENT(in) :: x
      IF (x /= x) THEN
         CALL log_warn("ignore_unused triggered (NaN)")
      END IF
   END SUBROUTINE ignore_unused

   SUBROUTINE log_open(filename, append)
      USE, INTRINSIC :: iso_fortran_env, ONLY: OUTPUT_UNIT, ERROR_UNIT
      CHARACTER(len=*), INTENT(in) :: filename
      LOGICAL, INTENT(in), OPTIONAL :: append
      LOGICAL :: do_append
      INTEGER :: ios

      do_append = .FALSE.
      IF (present(append)) do_append = append

      LOGFILE = trim(filename)

      IF (LOG_TO_FILE) CALL log_close()

      IF (do_append) THEN
         OPEN (newunit=LOGU, file=trim(LOGFILE), status="unknown", position="append", &
               action="write", iostat=ios)
      ELSE
         OPEN (newunit=LOGU, file=trim(LOGFILE), status="replace", action="write", iostat=ios)
      END IF

      IF (ios /= 0) THEN
         LOG_TO_FILE = .FALSE.
         WRITE (ERROR_UNIT, '(a,1x,a,1x,i0)') "[WARN]", "Could not open log file:", trim(LOGFILE), ios
         LOGU = -1
      ELSE
         LOG_TO_FILE = .TRUE.
         WRITE (OUTPUT_UNIT, '(a,1x,a)') "[INFO ]", "Logging to file: "//trim(LOGFILE)
         WRITE (LOGU, '(a,1x,a)') "[INFO ]", "Logging started"
      END IF
   END SUBROUTINE log_open

   SUBROUTINE log_close()
      USE, INTRINSIC :: iso_fortran_env, ONLY: ERROR_UNIT
      INTEGER :: ios

      IF (.NOT. LOG_TO_FILE) RETURN

      WRITE (LOGU, '(a)') " "
      WRITE (LOGU, '(a,1x,a)') "[INFO ]", "Logging ended"
      CLOSE (LOGU, iostat=ios)

      IF (ios /= 0) THEN
         WRITE (ERROR_UNIT, '(a,1x,a,1x,i0)') "[WARN]", "Could not close log file cleanly, iostat=", ios
      END IF

      LOG_TO_FILE = .FALSE.
      LOGU = -1
   END SUBROUTINE log_close

   SUBROUTINE log_raw(line)
      USE, INTRINSIC :: iso_fortran_env, ONLY: OUTPUT_UNIT
      CHARACTER(len=*), INTENT(in) :: line
      WRITE (OUTPUT_UNIT, '(a)') trim(line)
      IF (LOG_TO_FILE) WRITE (LOGU, '(a)') trim(line)
   END SUBROUTINE log_raw

END MODULE testing_T3ST_mod
