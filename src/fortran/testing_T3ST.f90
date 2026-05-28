module testing_T3ST_mod
  use constants
  use random_numbers              ! provides kx,ky,kz,w,ph,L and kxs,kys,kzs,ws,phs,Ls and USE_real
  use drift_kernels
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
  implicit none
  private

  public :: testing_T3ST
  public :: diagnose_saved_trajectories

  ! -------------------------
  ! Unified formatting knobs
  ! -------------------------
  integer, parameter :: WLBL = 52     ! label column width (bigger to avoid truncation)
  integer, parameter :: WTAG = 7      ! width for tags like "[ OK ]"
  integer, parameter :: WSEC = 70     ! banner width
  integer, parameter :: MAX_ALERTS = 256
  integer, parameter :: MSG_LEN = 220

  ! -------------------------
  ! Logging (tee: terminal + file)
  ! -------------------------
  integer,  save :: LOGU = -1
  logical,  save :: LOG_TO_FILE = .false.
  character(len=512), save :: LOGFILE = "testing_T3ST.log"

  ! Shared scratch buffer for formatted lines (available to all contained procs)
  character(len=1024), save :: TMP

  ! More legible scientific formats
  character(len=*), parameter :: FVAL = 'es10.3'  ! values
  character(len=*), parameter :: FRAT = 'es10.3'  ! ratios/tolerances
  character(len=*), parameter :: FINT = 'i0'      ! integers

  ! -------------------------
  ! Alert collection (so we can recap at the end)
  ! -------------------------
  integer :: n_alerts = 0
  character(len=MSG_LEN) :: alert_msgs(MAX_ALERTS)

contains
  subroutine diagnose_saved_trajectories(F0, FMtraj, betaFtraj, Htraj, Pctraj, &
                                         Obs_FF, Obs_FD, Obs_DF, Obs_DD, no_errors, fail_counts)
    real(rp), intent(in)  :: F0(:)
    real(rp), intent(in)  :: FMtraj(:, :), betaFtraj(:, :)
    real(rp), intent(in)  :: Htraj(:, :), Pctraj(:, :)
    real(rp), intent(in)  :: Obs_FF(:, :), Obs_FD(:, :), Obs_DF(:, :), Obs_DD(:, :)
    integer,  intent(out) :: no_errors
    integer,  intent(out), optional :: fail_counts(:)

    real(rp), parameter :: tol_identity = 1.0e-3_rp
    real(rp), parameter :: tol_obs_abs  = 1.0e-6_rp
    real(rp), parameter :: tol_obs_rel  = 1.0e-5_rp
    real(rp), parameter :: tol_H1       = 1.0e-5_rp
    real(rp), parameter :: tol_H2       = 2.0e-2_rp
    real(rp), parameter :: tol_Pc1      = 1.0e-4_rp
    real(rp), parameter :: eps_ref      = 1.0e-12_rp

    integer :: ntime, nsaved
    integer :: i, k
    integer :: before_errors
    real(rp) :: max_identity_err
    real(rp) :: value

    no_errors = 0
    if (present(fail_counts)) fail_counts = 0
    ntime = size(FMtraj, 1)
    nsaved = size(FMtraj, 2)

    if (ntime < 1 .or. nsaved < 1) then
      no_errors = no_errors + 1
      if (present(fail_counts)) fail_counts(1) = fail_counts(1) + 1
      return
    end if

    if (size(betaFtraj, 1) /= ntime .or. size(betaFtraj, 2) /= nsaved .or. size(F0) < nsaved) then
      no_errors = no_errors + 1
      if (present(fail_counts)) fail_counts(1) = fail_counts(1) + 1
      return
    end if

    max_identity_err = 0.0_rp
    do i = 1, nsaved
      if (F0(i) <= 0.0_rp) then
        no_errors = no_errors + 1
        if (present(fail_counts)) fail_counts(2) = fail_counts(2) + 1
        return
      end if

      do k = 1, ntime
        value = FMtraj(k, i) * exp(betaFtraj(k, i)) / F0(i)
        max_identity_err = max(max_identity_err, abs(value - 1.0_rp))
      end do
    end do

    before_errors = no_errors
    call check_scalar_tolerance("FM*exp(betaF)/F0 over saved trajectories", max_identity_err, tol_identity, no_errors)
    if (present(fail_counts)) fail_counts(3) = fail_counts(3) + no_errors - before_errors

    if (size(Htraj, 2) >= 1) then
      before_errors = no_errors
      call check_time_constancy("Htraj(:,1)", Htraj(:, 1), tol_H1, eps_ref, no_errors)
      if (present(fail_counts)) fail_counts(4) = fail_counts(4) + no_errors - before_errors
    else
      no_errors = no_errors + 1
      if (present(fail_counts)) fail_counts(4) = fail_counts(4) + 1
    end if

    if (size(Htraj, 2) >= 2) then
      before_errors = no_errors
      call check_time_constancy("Htraj(:,2)", Htraj(:, 2), tol_H2, eps_ref, no_errors)
      if (present(fail_counts)) fail_counts(5) = fail_counts(5) + no_errors - before_errors
    else
      no_errors = no_errors + 1
      if (present(fail_counts)) fail_counts(5) = fail_counts(5) + 1
    end if

    if (size(Pctraj, 2) >= 1) then
      before_errors = no_errors
      call check_time_constancy("Pctraj(:,1)", Pctraj(:, 1), tol_Pc1, eps_ref, no_errors)
      if (present(fail_counts)) fail_counts(6) = fail_counts(6) + no_errors - before_errors
    else
      no_errors = no_errors + 1
      if (present(fail_counts)) fail_counts(6) = fail_counts(6) + 1
    end if

    before_errors = no_errors
    call check_observable_identity("Obs_FF - Obs_FD = Obs_DF - Obs_DD", &
                                   Obs_FF - Obs_FD, Obs_DF - Obs_DD, &
                                   tol_obs_abs, tol_obs_rel, eps_ref, no_errors)
    if (present(fail_counts)) fail_counts(7) = fail_counts(7) + no_errors - before_errors

    before_errors = no_errors
    call check_observable_identity("Obs_FD - Obs_DD = Obs_FF - Obs_DF", &
                                   Obs_FD - Obs_DD, Obs_FF - Obs_DF, &
                                   tol_obs_abs, tol_obs_rel, eps_ref, no_errors)
    if (present(fail_counts)) fail_counts(8) = fail_counts(8) + no_errors - before_errors
  end subroutine diagnose_saved_trajectories

  !=============================================================================
  ! Main entry point
  !=============================================================================
  subroutine testing_T3ST(X, Y, Z, Vp, mu, F0, G0, FM, betaF, betaD, Vstar1, Vstar2, Vstar3, no_errors)
    real(rp), intent(in)  :: X(:), Y(:), Z(:), Vp(:), mu(:)
    real(rp), intent(in)  :: F0(:), G0(:), FM(:), betaF(:), betaD(:)
    real(rp), intent(in)  :: Vstar1, Vstar2, Vstar3
    integer,  intent(out) :: no_errors

    integer  :: Np_partial, Nt_partial, np_avail
    real(rp) :: dt, ratio32

    real(rp), allocatable :: Xh(:,:), Yh(:,:), Zh(:,:), Vph(:,:), muh(:,:)
    real(rp), allocatable :: Hh(:,:), Pch(:,:), c1h(:,:), c2h(:,:), c3h(:,:),  Vtx(:,:),  VFx(:,:)

    integer  :: nstep_eff
    real(rp) :: T, compare
    real(rp), parameter :: epsH = 1.0e-300_rp

    no_errors = 0
    ! -------------------------
    ! Controls (hard-coded here)
    ! -------------------------
    Np_partial = 10000
    Nt_partial = 200
    np_avail   = size(X)
    Np_partial = min(Np_partial, np_avail)
    dt         = real(tmax / Nt, rp)   ! base dt used by rk4_propagation below

    allocate(Xh (1:Nt_partial, 1:Np_partial), &
             Yh (1:Nt_partial, 1:Np_partial), &
             Zh (1:Nt_partial, 1:Np_partial), &
             Vph(1:Nt_partial, 1:Np_partial), &
             muh(1:Nt_partial, 1:Np_partial), &
             Hh (1:Nt_partial, 1:Np_partial), &
             Pch(1:Nt_partial, 1:Np_partial), &
             c1h(1:Nt_partial, 1:Np_partial), &
             c2h(1:Nt_partial, 1:Np_partial), &
             c3h(1:Nt_partial, 1:Np_partial),&
             Vtx(1:Nt_partial, 1:Np_partial), &
             VFx(1:Nt_partial, 1:Np_partial))

    call alerts_reset()
    call banner_begin()
    call log_open("testing_T3ST.log")   ! or any filename you want

    ! -------------------------
    ! Inputs: NaN checks
    ! -------------------------
    call section("Inputs: NaN checks")
    call check_nan("X",      X,      no_errors)
    call check_nan("Y",      Y,      no_errors)
    call check_nan("Z",      Z,      no_errors)
    call check_nan("Vp",     Vp,     no_errors)
    call check_nan("mu",     mu,     no_errors)
    call check_nan("F0",     F0,     no_errors)
    call check_nan("G0",     G0,     no_errors)
    call check_nan("FM",     FM,     no_errors)
    call check_nan("betaF",  betaF,  no_errors)
    call check_nan("betaD",  betaD,  no_errors)
    call check_nan("Vstar1", Vstar1, no_errors)
    call check_nan("Vstar2", Vstar2, no_errors)
    call check_nan("Vstar3", Vstar3, no_errors)

    ! -------------------------
    ! NaN checks: waves (choose set based on USE_real)
    ! (Skip L/Ls entirely)
    ! -------------------------
    call section("Waves: NaN checks")
    if (USE_real == ON) then
      call log_info("Wave-set selected: s (USE_real = ON)")
      call check_nan("kxs", kxs, no_errors)
      call check_nan("kys", kys, no_errors)
      call check_nan("kzs", kzs, no_errors)
      call check_nan("ws",  ws,  no_errors)
      call check_nan("phs", phs, no_errors)
    else
      call log_info("Wave-set selected: base (USE_real = OFF)")
      call check_nan("kx",  kx,  no_errors)
      call check_nan("ky",  ky,  no_errors)
      call check_nan("kz",  kz,  no_errors)
      call check_nan("w",   w,   no_errors)
      call check_nan("ph",  ph,  no_errors)
    end if

    ! -------------------------
    ! Vstar checks
    ! -------------------------
    call section("Vstar checks")
    call log_vals3("[INFO ]", "Vstar1 Vstar2 Vstar3", Vstar1, Vstar2, Vstar3)

    if (abs(Vstar1) > 1.0e-9_rp) then
      no_errors = no_errors + 1
      call log_alert("Vstar1 not ~ 0 (radial drift present)")
    else
      call log_ok("Vstar1 ~ 0 (radial drift OK)")
    end if

    if (abs(Vstar2) <= tiny(1.0_rp)) then
      no_errors = no_errors + 1
      call log_alert("Vstar2 ~ 0 (cannot evaluate Vstar3/Vstar2)")
    else
      ratio32 = Vstar3 / Vstar2
      call log_ratio("[INFO ]", "Vstar3/Vstar2", abs(ratio32), 0.1_rp)

      if (abs(ratio32) > 0.1_rp) then
        no_errors = no_errors + 1
        call log_alert("Large parallel drift (|Vstar3/Vstar2| > 0.1)")
      else
        call log_ok("Parallel drift OK (|Vstar3/Vstar2| <= 0.1)")
      end if
    end if

    ! -------------------------
    ! Input consistency checks
    ! -------------------------
    call section("Input consistency checks")
    call check_inputs(X, Y, Z, Vp, mu, F0, G0, FM, betaF, betaD, a0, 1.0e-12_rp, no_errors)

    ! -------------------------
    ! Wave moment/symmetry checks (same switch)
    ! -------------------------
    call section("Wave moment checks")
    if (USE_real == ON) then
      call check_waves_any(kxs, kys, kzs, ws, phs, Vstar2, Vstar3, no_errors)
    else
      call check_waves_any(kx,  ky,  kz,  w,  ph,  Vstar2, Vstar3, no_errors)
    end if

!!!!
!!!let's do something for eulerian correlation
!!!!

    ! =====================================================================
    ! Drift / RK4 propagation
    ! =====================================================================
    call section("RK4 propagation")
    call log_val("[INFO ]", "Np_partial", real(Np_partial, rp))
    call log_val("[INFO ]", "Nt_partial", real(Nt_partial, rp))
    call log_val("[INFO ]", "dt", dt)


    call rk4_propagation(X, Y, Z, Vp, mu, betaF, betaD, dt, Nt_partial, Np_partial, &
                         Xh, Yh, Zh, Vph, muh, Hh, Pch, c1h, c2h, c3h, Vtx, VFx)
                         
    nstep_eff = max(1, Nt_partial - 1)
    T         = dt * real(nstep_eff, rp)

    ! =====================================================================
    ! Energy proxy diagnostics (d ln H / dt)
    ! =====================================================================
    if (USE_freq /= OFF) then
      call log_line("[ALERT]", "Non-zero frequencies (USE_freq /= OFF): forget conservation")
    end if

    call section("Energy conservation proxy ")
    call print_dlnH_stats(Hh, Np_partial, Nt_partial, T, epsH, no_errors)
    call section("Toroidal(canonical) momentul conservation proxy")
    call print_dlnPc_stats(Pch, Np_partial, Nt_partial, T, epsH, no_errors)

    ! =====================================================================
    ! Gaussianity / moment checks: whole (time x particles)
    ! =====================================================================
    call section("Gaussianity summary: all times x particles (kurt/3 = 1 for Gaussian)")
    compare = 1.0_rp
    
    call print_gaussianity_block("phi",  c1h, Np_partial, Nt_partial, no_errors, compare)
        compare = 1.0_rp
    call print_gaussianity_block("phix", c2h, Np_partial, Nt_partial, no_errors, compare)
        compare = 1.0_rp
    call print_gaussianity_block("phiy", c3h, Np_partial, Nt_partial, no_errors, compare)

    call section("Gaussianity drivers (qualitative)")
    call log_info("Larmor & Frequencies -> NO effects on Gaussianity of fields")
    call log_info("Balloning breaks phix, phiz Gaussianity as I_0[4/Lbalon^2]/I_0[2/Lbalon^2]^2/2")
    call log_info("Balloning breaks amplitude as: I_0(1/lbalonz^2)*exp(-lbalonz^2)")
    call log_info("Tilting breaks phix/phiz Gaussianit when particles are z distributey; k_x^eff = kx+ky a --> ??")

    ! =====================================================================
    ! Gaussianity / moment checks: time slice (k = 1)
    ! =====================================================================
    call section("Gaussianity summary: time slice k=1 (kurt/3 = 1 for Gaussian)")
     compare = 1.0_rp
         if (USE_balloon.eq.ON) then
      compare = compare**exp(-1.0_rp/lbalonz**2)*(1 + 1.0_rp/lbalonz**2/4.0_rp + 1.0_rp/lbalonz**4/64.0_rp)
    endif
    call print_gaussianity_slice("phi",  c1h, 1, Np_partial, no_errors, compare)
           compare = 1.0_rp / lambdax**4 &
        + (pi**4 / 5.0_rp) * (3.0_rp/lambday**2 + k0i**2)**2 &
          * (2.0_rp*C2/(C1*C3))**4 * 9.0_rp * (r00/a0**2)**4 &
        + (2.0_rp*pi**2 / 3.0_rp) * (1.0_rp/lambdax**2) &
          * (3.0_rp/lambday**2 + k0i**2) &
          * (2.0_rp*C2/(C1*C3))**2 * 3.0_rp * (r00/a0**2)**2
          call print_gaussianity_slice("phix", c2h, 1, Np_partial, no_errors, compare)
 
    call print_gaussianity_slice("phiy", c3h, 1, Np_partial, no_errors, compare)

    if (USE_real.eq.OFF) then
      compare = ((R0/rhoi*Phi)**2)*sum(ky**2)/real(Np*Nc,rp)
    elseif (USE_real.eq.ON) then
      compare = ((R0/rhoi*Phi)**2)*sum(ky**2)/real(Np*Nc,rp)
    endif
    if (USE_balloon.eq.ON) then
      compare = compare**exp(-1.0_rp/lbalonz**2)*(1 + 1.0_rp/lbalonz**2/4.0_rp + 1.0_rp/lbalonz**4/64.0_rp)
    endif
    ! Larmor

    call log_raw(' --------------------------------------------------------------')
    call log_raw(' ------ A ballistic/QL estimation of turbulent transport ------')
    write(TMP,'(A,' // FVAL // ',A,' // FVAL // ',A,' // FVAL // ')') &
      '  E[Vtx] =', sum(Vtx(1,:))/Np_partial, ' E[Vtx^2] =', sum(Vtx(1,:)**2)/Np_partial,  ' compared to ', compare
    call log_raw(TMP)
    call log_raw(' --------------------------------------------------------------')
    call log_raw(' ------ A ballistic/QL estimation of neoclassical transport ------')
    write(TMP,'(A,' // FVAL // ',A,' // FVAL // ',A,' // FVAL // ')') &
      '  E[Vnx] =', sum(VFx(1,:)-Vtx(1,:))/Np_partial, ' E[Vnx^2] =', sum((VFx(1,:)-Vtx(1,:))**2)/Np_partial,  ' compared to (only if theta=pi/2)', compare
    call log_raw(TMP)

    ! -------------------------
    ! Final recap + banner end
    ! -------------------------
    call alerts_recap()
    call log_close()
    call banner_end(no_errors)

  end subroutine testing_T3ST


  !=============================================================================
  ! RK4 propagation (unchanged numerics; formatting/indent only)
  !=============================================================================
  subroutine rk4_propagation(X, Y, Z, Vp, mu, betaF, betaD, dt, Nt_partial, Np_partial, &
                             Xh, Yh, Zh, Vph, muh, Hh, Pch, c1h, c2h, c3h, Vtxall, VFxall)

    real(rp), intent(in) :: X(:), Y(:), Z(:), Vp(:), mu(:), betaF(:), betaD(:)
    integer,  intent(in) :: Np_partial, Nt_partial
    real(rp), intent(in) :: dt

    real(rp), intent(out) :: Xh(:,:), Yh(:,:), Zh(:,:), Vph(:,:), muh(:,:)
    real(rp), intent(out) :: Hh(:,:), Pch(:,:), c1h(:,:), c2h(:,:), c3h(:,:), Vtxall(:,:), VFxall(:,:)

    real(rp), pointer, contiguous :: Qx(:), Qy(:), Qz(:), Qw(:), Qph(:)

    real(rp) :: xi, yi, zi, mui, vpi, betaFi, betaDi
    real(rp) :: vx, vy, vz, vm, ap, vbF, vbD, kin
    real(rp) :: Wx, Wy, Wz, Wm, Wp, WbF, WbD
    real(rp) :: q1, q2, q3
    real(rp) :: Hi, FMi, B, Vtx, VFx, Pc
    real(rp) :: check_1, check_2, check_3, vs_1, vs_2, vs_3, time
    integer  :: ias, k

    do ias = 1, Np_partial

      ! Select per-particle or single-spectrum arrays
      if (USE_real == OFF) then
        Qx    => kx(:, ias)
        Qy    => ky(:, ias)
        Qz    => kz(:, ias)
        Qw    => w(:,  ias)
        Qph   => ph(:, ias)
      else
        Qx    => kxs
        Qy    => kys
        Qz    => kzs
        Qw    => ws
        Qph   => phs
      end if

      xi  = X(ias)
      yi  = Y(ias)
      zi  = Z(ias)
      mui = mu(ias)
      vpi = Vp(ias)
      betaFi = betaF(ias)
      betaDi = betaD(ias)

      do k = 1, Nt_partial

        time = dt * real(k - 1, rp)

call gyrocenter_drifts(xi, yi, zi, vpi, mui, q1, q2, q3, time,                       &
            vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx,     &
            check_1, check_2, check_3, vs_1, vs_2, vs_3, Qx, Qy, Qz, Qw, Qph)
            
        Xh(k, ias)  = xi
        Yh(k, ias)  = yi
        Zh(k, ias)  = zi
        Vph(k, ias) = vpi
        muh(k, ias) = mui
        Hh(k, ias)  = Hi
        Pch(k, ias) = Pc
        c1h(k, ias) = check_1
        c2h(k, ias) = check_2
        c3h(k, ias) = check_3
        Vtxall(k,ias) = Vtx
        VFxall(k,ias) = VFx

        Wx = vx / 6.0_rp
        Wy = vy / 6.0_rp
        Wz = vz / 6.0_rp
        Wm = vm / 6.0_rp
        WbF = vbF / 6.0_rp
        WbD = vbD / 6.0_rp
        Wp = ap / 6.0_rp

call gyrocenter_drifts(xi + vx*dt/2.0_rp, yi + vy*dt/2.0_rp, zi + vz*dt/2.0_rp,      &
            vpi + ap*dt/2.0_rp, mui + vm*dt/2.0_rp,                       &
            q1, q2, q3, time + dt/2.0_rp,                                  &
            vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx,     &
            check_1, check_2, check_3, vs_1, vs_2, vs_3, Qx, Qy, Qz, Qw, Qph)
                          
        Wx = Wx + vx / 3.0_rp
        Wy = Wy + vy / 3.0_rp
        Wz = Wz + vz / 3.0_rp
        Wm = Wm + vm / 3.0_rp
        WbF = WbF + vbF / 3.0_rp
        WbD = WbD + vbD / 3.0_rp
        Wp = Wp + ap / 3.0_rp


! RK2b
call gyrocenter_drifts(xi + vx*dt/2.0_rp, yi + vy*dt/2.0_rp, zi + vz*dt/2.0_rp,      &
            vpi + ap*dt/2.0_rp, mui + vm*dt/2.0_rp,                       &
            q1, q2, q3, time + dt/2.0_rp,                                  &
            vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx,     &
            check_1, check_2, check_3, vs_1, vs_2, vs_3, Qx, Qy, Qz, Qw, Qph)


        Wx = Wx + vx / 3.0_rp
        Wy = Wy + vy / 3.0_rp
        Wz = Wz + vz / 3.0_rp
        Wm = Wm + vm / 3.0_rp
        WbF = WbF + vbF / 3.0_rp
        WbD = WbD + vbD / 3.0_rp
        Wp = Wp + ap / 3.0_rp

! RK4
call gyrocenter_drifts(xi + vx*dt, yi + vy*dt, zi + vz*dt,                            &
            vpi + ap*dt, mui + vm*dt,                                      &
            q1, q2, q3, time + dt,                                         &
            vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx,     &
            check_1, check_2, check_3, vs_1, vs_2, vs_3, Qx, Qy, Qz, Qw, Qph)
            
            
        Wx = Wx + vx / 6.0_rp
        Wy = Wy + vy / 6.0_rp
        Wz = Wz + vz / 6.0_rp
        Wm = Wm + vm / 6.0_rp
        WbF = WbF + vbF / 6.0_rp
        WbD = WbD + vbD / 6.0_rp
        Wp = Wp + ap / 6.0_rp

        xi  = xi  + dt * Wx
        yi  = yi  + dt * Wy
        zi  = zi  + dt * Wz
        vpi = vpi + dt * Wp
        mui = mui + dt * Wm
        betaFi = betaFi + dt * WbF
        betaDi = betaDi + dt * WbD

      end do
    end do
  end subroutine rk4_propagation


  !=============================================================================
  ! Unified NaN checker (assumed-rank)
  !=============================================================================
  subroutine check_nan(name, x, no_errors)
    character(len=*), intent(in)    :: name
    real(rp),         intent(in)    :: x(..)
    integer,          intent(inout) :: no_errors

    logical :: has_nan

    has_nan = .false.
    select rank (x)
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
      call log_warn("check_nan: unsupported rank for " // trim(name))
      return
    end select

    if (has_nan) then
      no_errors = no_errors + 1
      call log_alert(trim(name) // " contains NaN")
    else
      call log_ok("No NaNs in " // trim(name))
    end if
  end subroutine check_nan

  !=============================================================================
  ! Section header
  !=============================================================================
  subroutine section(title)
    character(len=*), intent(in) :: title
    call log_raw(" ")
    call log_line("[INFO ]", repeat("-", 1) // " " // trim(title) // " " // repeat("-", max(0, WSEC - 9 - len_trim(title))))
  end subroutine section


  !=============================================================================
  ! Banners
  !=============================================================================
  subroutine banner_begin()
    call log_raw(" ")
    call log_raw(repeat("=", WSEC))
    call log_raw(" TESTING / DIAGNOSTICS")
    call log_raw(repeat("=", WSEC))
  end subroutine banner_begin

  subroutine banner_end(no_errors)
    integer, intent(in) :: no_errors
    character(len=256) :: tmp2
    call log_raw(repeat("-", WSEC))
    write(tmp2,'(a,1x,'//FINT//')') "[INFO ] Total number of alerts:", no_errors
    call log_raw(tmp2)
    call log_raw(repeat("=", WSEC))
  end subroutine banner_end


  !=============================================================================
  ! Input checks (uniform messages + formatted numbers)
  !=============================================================================
  subroutine check_inputs(X, Y, Z, Vp, mu, F0, G0, FM, betaF, betaD, a1, tol, no_errors)
    real(rp), intent(in)    :: X(:), Y(:), Z(:), Vp(:), mu(:)
    real(rp), intent(in)    :: F0(:), G0(:), FM(:), betaF(:), betaD(:)
    real(rp), intent(in)    :: a1, tol
    integer,  intent(inout) :: no_errors

    integer  :: n, nbad
    real(rp) :: pii, max_betaF_err, max_betaD_err

    n   = size(X)
    pii = 3.14159265358979323846_rp

    if (any([ size(Y), size(Z), size(Vp), size(mu), size(F0), size(G0), size(FM), size(betaF), size(betaD) ] /= n)) then
      no_errors = no_errors + 1
      call log_alert("Input arrays have inconsistent sizes")
      write(TMP,'(a,1x,a,": X=",1x,'//FINT//',", Y=",1x,'//FINT//',", Z=",1x,'//FINT//', &
                   ", Vp=",1x,'//FINT//',", mu=",1x,'//FINT//',", F0=",1x,'//FINT//', &
                   ", G0=",1x,'//FINT//',", FM=",1x,'//FINT//',", betaF=",1x,'//FINT//', &
                   ", betaD=",1x,'//FINT//')') &
        "[INFO ]", pad_label("sizes"), size(X), size(Y), size(Z), size(Vp), size(mu), &
        size(F0), size(G0), size(FM), size(betaF), size(betaD)
      call log_raw(TMP)
      return
    else
      call log_ok("All input arrays have consistent size")
      write(TMP,'(a,1x,a,":",1x,'//FINT//')') "[INFO ]", pad_label("n"), n
      call log_raw(TMP)
    end if

    nbad = count(X <= 0.0_rp)
    if (nbad > 0) then
      no_errors = no_errors + 1
      call log_alert("Condition failed: X must be > 0")
      call log_count("violations", nbad, n)
    else
      call log_ok("X > 0 everywhere")
    end if

    nbad = count(abs(Y) >= 1.0_rp)
    if (nbad > 0) then
      no_errors = no_errors + 1
      call log_alert("Condition failed: abs(Y) must be < 1")
      call log_count("violations", nbad, n)
    else
      call log_ok("abs(Y) < 1 everywhere")
    end if

    nbad = count((Z <= -pii) .or. (Z >= pii))
    if (nbad > 0) then
      no_errors = no_errors + 1
      call log_alert("Condition failed: Z must be in (-pi, pi)")
      call log_count("violations", nbad, n)
    else
      call log_ok("Z in (-pi, pi) everywhere")
    end if

    nbad = count(mu <= 0.0_rp)
    if (nbad > 0) then
      no_errors = no_errors + 1
      call log_alert("Condition failed: mu must be > 0")
      call log_count("violations", nbad, n)
    else
      call log_ok("mu > 0 everywhere")
    end if

    nbad = count(F0 <= 0.0_rp)
    if (nbad > 0) then
      no_errors = no_errors + 1
      call log_alert("Condition failed: F0 must be > 0")
      call log_count("violations", nbad, n)
    else
      call log_ok("F0 > 0 everywhere")
    end if

    nbad = count(G0 <= 0.0_rp)
    if (nbad > 0) then
      no_errors = no_errors + 1
      call log_alert("Condition failed: G0 must be > 0")
      call log_count("violations", nbad, n)
    else
      call log_ok("G0 > 0 everywhere")
    end if

    nbad = count(FM <= 0.0_rp)
    if (nbad > 0) then
      no_errors = no_errors + 1
      call log_alert("Condition failed: FM must be > 0")
      call log_count("violations", nbad, n)
    else
      call log_ok("FM > 0 everywhere")
    end if

    if (all((F0 > 0.0_rp) .and. (FM > 0.0_rp))) then
      max_betaF_err = maxval(abs(FM * exp(betaF) / F0 - 1.0_rp))
      max_betaD_err = maxval(abs(FM * exp(betaD) / F0 - 1.0_rp))

      if (max_betaF_err > tol) then
        no_errors = no_errors + 1
        call log_alert("Condition failed: FM*exp(betaF)/F0 must be 1")
      else
        call log_ok("FM*exp(betaF)/F0 = 1 within tolerance")
      end if
      call log_val("[INFO ]", "max betaF identity error", max_betaF_err)

      if (max_betaD_err > tol) then
        no_errors = no_errors + 1
        call log_alert("Condition failed: FM*exp(betaD)/F0 must be 1")
      else
        call log_ok("FM*exp(betaD)/F0 = 1 within tolerance")
      end if
      call log_val("[INFO ]", "max betaD identity error", max_betaD_err)
      call log_val("[INFO ]", "tol", tol)
    end if

    nbad = count(sqrt((X - 1.0_rp)**2 + Y**2) >= a1)
    if (nbad > 0) then
      no_errors = no_errors + 1
      call log_alert("Condition failed: sqrt((X-1)^2 + Y^2) < a0")
      call log_count("violations", nbad, n)
      call log_val("[INFO ]", "a0", a1)
    else
      call log_ok("sqrt((X-1)^2 + Y^2) < a0 everywhere")
      call log_val("[INFO ]", "a0", a1)
    end if

  end subroutine check_inputs


  !=============================================================================
  ! Wave checks (rank-agnostic) WITHOUT L/Ls
  !=============================================================================
  subroutine check_waves_any(kxA, kyA, kzA, wA, phA, Vstar2, Vstar3, no_errors)
    real(rp), intent(in)    :: kxA(..), kyA(..), kzA(..), wA(..), phA(..)
    real(rp), intent(in)    :: Vstar2, Vstar3
    integer,  intent(inout) :: no_errors

    real(rp), parameter :: tiny_rp = 1.0e-12_rp
    real(rp), parameter :: eps_ref = 1.0e-30_rp

    real(rp) :: m1(5), m2(5), m4(5)
    real(rp) :: tol_mean(5), tol_m2(5), tol_m4(5)
    real(rp) :: m2_ref(5), m4_ref(5)
    real(rp) :: ratio, pii
    real(rp) :: t_m2_ky, t_m2_kz, t_m4_ky, t_m4_kz
    real(rp) :: fac

    call log_info("Wave diagnostics (moments / target checks)")

    pii = 3.14159265358979323846_rp

    tol_mean = 0.05_rp
    tol_m2   = 0.05_rp
    tol_m4   = 0.10_rp


    m2_ref(2) = 3.0_rp / (lambday*lambday) + k0i*k0i
    m2_ref(3) = (erfc(lambdaz/sqrt(2.0_rp)) + 3.0_rp*erfc(lambdaz*sqrt(2.0_rp)))/C3**2!1.0_rp / (lambdaz*lambdaz) ! erfc(lambdaz/sqrt(2.0_rp))
    m2_ref(5) = (pii*pii) / 3.0_rp

        ! Radial spectrum selection
        IF (x_corr == 1) THEN
            m2_ref(1) = 1.0_rp / (lambdax**2)
            m4_ref(1) = 3.0_rp / (lambdax**4)
        ELSEIF (x_corr == 2) THEN
            m2_ref(1) = (10.0_rp / atan(10.0_rp)-1.0_rp)/ (lambdax**2)
            m4_ref(1) = (1.0_rp + 10.0_rp*(10.0_rp**2-3)/3.0_rp/atan(10.0_rp))/ (lambdax**4)
        END IF

        ! Temporal spectrum selection
        IF (t_corr == 1) THEN
            m2_ref(4) = 1.0_rp / (tauc**2)
            m4_ref(4) = 3.0_rp / (tauc**4)
        ELSEIF (t_corr == 2) THEN
            m2_ref(4) = (10.0_rp / atan(10.0_rp)-1.0_rp)/ (tauc**2)
            m4_ref(4) = (1.0_rp + 10.0_rp*(10.0_rp**2-3)/3.0_rp/atan(10.0_rp))/ (tauc**4)
        END IF

    m4_ref(2) = (15.0_rp/(lambday**4) + 10.0_rp*k0i**2/(lambday**2) + k0i**4)
    m4_ref(3) = (erfc(lambdaz/sqrt(2.0_rp)) + 15.0_rp*erfc(lambdaz*sqrt(2.0_rp)) + 65.0_rp*erfc(lambdaz*3.0_rp/sqrt(2.0_rp)))/C3**4
    m4_ref(5) = (pii**4) / 5.0_rp

    t_m2_ky = m2_ref(2)
    t_m2_kz = m2_ref(3)
    t_m4_ky = m4_ref(2)
    t_m4_kz = m4_ref(3)

    fac = (Ln*(rhoi/R0/(1.0_rp+r00))**2)

    m2_ref(4) = m2_ref(4) + (fac**2) * (Vstar2**2*t_m2_ky + Vstar3**2*t_m2_kz)

    m4_ref(4) = m4_ref(4)                                                        &
              + (fac**4) * (Vstar2**4*t_m4_ky + Vstar3**4*t_m4_kz                         &
              + 2.0_rp*Vstar2**2*Vstar3**2*t_m2_kz*t_m2_ky)                               &
              + 2.0_rp/(tauc**2) * (fac**2) * (Vstar2**2*t_m2_ky + Vstar3**2*t_m2_kz)

    call moments_ar(kxA, m1(1), m2(1), m4(1))
    call moments_ar(kyA, m1(2), m2(2), m4(2))
    call moments_ar(kzA, m1(3), m2(3), m4(3))
    call moments_ar(wA,  m1(4), m2(4), m4(4))
    call moments_ar(phA, m1(5), m2(5), m4(5))

    ratio = abs(m1(1)) / max(sqrt(m2(1)), tiny_rp)
    call check_ratio("kx mean ~ 0", ratio, tol_mean(1), no_errors)

    ratio = abs(m1(2)) / max(sqrt(m2(2)), tiny_rp)
    call check_ratio("ky mean ~ 0", ratio, tol_mean(2), no_errors)

    ratio = abs(m1(3)) / max(sqrt(m2(3)), tiny_rp)
    call check_ratio("kz mean ~ 0", ratio, tol_mean(3), no_errors)

    ratio = abs(m1(4)) / max(sqrt(m2(4)), tiny_rp)
    call check_ratio("w  mean ~ 0", ratio, tol_mean(4), no_errors)

    ratio = abs(m1(5)) / max(sqrt(m2(5)), tiny_rp)
    call check_ratio("ph mean ~ 0", ratio, tol_mean(5), no_errors)

    call log_info("Second moments (targets)")
    call check_moment("kx^2 target", m2(1), m2_ref(1), tol_m2(1), no_errors, eps_ref)
    call check_moment("ky^2 target", m2(2), m2_ref(2), tol_m2(2), no_errors, eps_ref)
    call check_moment("kz^2 target", m2(3), m2_ref(3), tol_m2(3), no_errors, eps_ref)
    call check_moment("w^2  target", m2(4), m2_ref(4), tol_m2(4), no_errors, eps_ref)
    call check_moment("ph^2 target", m2(5), m2_ref(5), tol_m2(5), no_errors, eps_ref)

    call log_info("Fourth moments (targets)")
    call check_moment("kx^4 target", m4(1), m4_ref(1), tol_m4(1), no_errors, eps_ref)
    call check_moment("ky^4 target", m4(2), m4_ref(2), tol_m4(2), no_errors, eps_ref)
    call check_moment("kz^4 target", m4(3), m4_ref(3), tol_m4(3), no_errors, eps_ref)
    call check_moment("w^4  target", m4(4), m4_ref(4), tol_m4(4), no_errors, eps_ref)
    call check_moment("ph^4 target", m4(5), m4_ref(5), tol_m4(5), no_errors, eps_ref)

  end subroutine check_waves_any


  subroutine check_moment(label, val, ref, tol, no_errors, eps_ref)
    character(len=*), intent(in)    :: label
    real(rp),         intent(in)    :: val, ref, tol, eps_ref
    integer,          intent(inout) :: no_errors
    real(rp) :: relerr

    relerr = abs(val - ref) / max(abs(ref), eps_ref)

    if (relerr > tol) then
      no_errors = no_errors + 1
      call log_moment("[ALERT]", label, val, ref, relerr, tol)
    else
      call log_moment("[  OK ]", label, val, ref, relerr, tol)
    end if
  end subroutine check_moment

  subroutine check_scalar_tolerance(label, err, tol, no_errors)
    character(len=*), intent(in)    :: label
    real(rp),         intent(in)    :: err, tol
    integer,          intent(inout) :: no_errors

    if (ieee_is_nan(err) .or. err > tol) then
      no_errors = no_errors + 1
    end if
  end subroutine check_scalar_tolerance

  subroutine check_time_constancy(label, x, tol, eps_ref, no_errors)
    character(len=*), intent(in)    :: label
    real(rp),         intent(in)    :: x(:)
    real(rp),         intent(in)    :: tol, eps_ref
    integer,          intent(inout) :: no_errors

    real(rp) :: err
    real(rp) :: ref

    if (size(x) < 1) then
      no_errors = no_errors + 1
      return
    end if

    if (any(ieee_is_nan(x))) then
      no_errors = no_errors + 1
      return
    end if

    ref = x(1)
    err = maxval(abs(x - ref)) / max(abs(ref), eps_ref)
    call check_scalar_tolerance(trim(label) // " time constancy", err, tol, no_errors)
  end subroutine check_time_constancy

  subroutine check_observable_identity(label, lhs, rhs, abs_tol, rel_tol, eps_ref, no_errors)
    character(len=*), intent(in)    :: label
    real(rp),         intent(in)    :: lhs(:, :), rhs(:, :)
    real(rp),         intent(in)    :: abs_tol, rel_tol, eps_ref
    integer,          intent(inout) :: no_errors

    real(rp) :: err, scale

    if (size(lhs, 1) /= size(rhs, 1) .or. size(lhs, 2) /= size(rhs, 2)) then
      no_errors = no_errors + 1
      return
    end if

    if (any(ieee_is_nan(lhs)) .or. any(ieee_is_nan(rhs))) then
      no_errors = no_errors + 1
      return
    end if

    err = maxval(abs(lhs - rhs))
    scale = max(maxval(abs(lhs)), maxval(abs(rhs)), eps_ref)

    if (err > abs_tol + rel_tol * scale) then
      no_errors = no_errors + 1
    end if

    call ignore_unused(real(len_trim(label), rp))
  end subroutine check_observable_identity


  subroutine log_moment(tag, label, val, ref, relerr, tol)
    character(len=*), intent(in) :: tag, label
    real(rp),         intent(in) :: val, ref, relerr, tol

    write(TMP,'(a,1x,a,": val=",1x,'//FVAL//',", ref=",1x,'//FVAL//',", rel=",1x,'//FRAT//',", tol=",1x,'//FRAT//')') &
      tag, pad_label(label), val, ref, relerr, tol
    call log_raw(TMP)

    if (trim(tag) == "[ALERT]") call alerts_add(trim(label))
  end subroutine log_moment


  !=============================================================================
  ! Moments for assumed-rank arrays (rank 0/1/2 supported)
  !=============================================================================
  subroutine moments_ar(x, m1, m2, m4)
    real(rp), intent(in)  :: x(..)
    real(rp), intent(out) :: m1, m2, m4

    integer  :: n
    real(rp) :: rn

    select rank (x)
    rank (0)
      n  = 1
      rn = 1.0_rp
      m1 = x
      m2 = x*x
      m4 = (x*x)*(x*x)

    rank (1)
      n  = size(x)
      rn = real(n, rp)
      m1 = sum(x)    / rn
      m2 = sum(x**2) / rn
      m4 = sum(x**4) / rn

    rank (2)
      n  = size(x)
      rn = real(n, rp)
      m1 = sum(x)    / rn
      m2 = sum(x**2) / rn
      m4 = sum(x**4) / rn

    rank default
      call log_warn("moments_ar: unsupported rank encountered")
      m1 = huge(1.0_rp)
      m2 = huge(1.0_rp)
      m4 = huge(1.0_rp)
    end select
  end subroutine moments_ar


  !=============================================================================
  ! Energy proxy: d(ln H)/dt finite difference stats
  !=============================================================================
  subroutine print_dlnH_stats(Hh, Np, Nt, T, epsH, no_errors)
    real(rp), intent(in)    :: Hh(:,:)
    integer,  intent(in)    :: Np, Nt
    real(rp), intent(in)    :: T, epsH
    integer,  intent(inout) :: no_errors

    real(rp), allocatable :: v(:)
    real(rp) :: mean_v, rms_v, std_v, maxabs_v, meanabs_v
    integer  :: imax(1)

    allocate(v(Np))

    v = (Hh(Nt, :) - Hh(1, :)) / max(abs(Hh(1, :)), epsH) / T

    mean_v    = sum(v) / real(Np, rp)
    rms_v     = sqrt(sum(v**2) / real(Np, rp))
    std_v     = sqrt(max(0.0_rp, rms_v**2 - mean_v**2))
    maxabs_v  = maxval(abs(v))
    meanabs_v = sum(abs(v)) / real(Np, rp)

    call log_val("[INFO ]", "mean (H(tf)/H(ti)-1)/(tf-ti)",    mean_v)
    call log_val("[INFO ]", "rms  (H(tf)/H(ti)-1)/(tf-ti)",    rms_v)
    call log_val("[INFO ]", "std  (H(tf)/H(ti)-1)/(tf-ti)",    std_v)
    call log_val("[INFO ]", "max |(H(tf)/H(ti)-1)/(tf-ti)|",   maxabs_v)
    call log_val("[INFO ]", "mean|(H(tf)/H(ti)-1)/(tf-ti)|",   meanabs_v)

    imax = maxloc(abs(v))
    write(TMP,'(a,1x,a,":",1x,'//FINT//')') "[INFO ]", pad_label("argmax |(H(tf)/H(ti)-1)/(tf-ti)| (particle index)"), imax(1)
    call log_raw(TMP)

    if (maxabs_v > 1.0e-1_rp) then
      no_errors = no_errors + 1
      call log_alert("Energy proxy: large max |(H(tf)/H(ti)-1)/(tf-ti)| (threshold 1e-1)")
    end if

    deallocate(v)
  end subroutine print_dlnH_stats


  !=============================================================================
  ! Energy proxy: d(ln H)/dt finite difference stats
  !=============================================================================
  subroutine print_dlnPc_stats(Hh, Np, Nt, T, epsH, no_errors)
    real(rp), intent(in)    :: Hh(:,:)
    integer,  intent(in)    :: Np, Nt
    real(rp), intent(in)    :: T, epsH
    integer,  intent(inout) :: no_errors

    real(rp), allocatable :: v(:)
    real(rp) :: mean_v, rms_v, std_v, maxabs_v, meanabs_v
    integer  :: imax(1)

    allocate(v(Np))

    v = (Hh(Nt, :) - Hh(1, :)) / max(abs(Hh(1, :)), epsH) / T

    mean_v    = sum(v) / real(Np, rp)
    rms_v     = sqrt(sum(v**2) / real(Np, rp))
    std_v     = sqrt(max(0.0_rp, rms_v**2 - mean_v**2))
    maxabs_v  = maxval(abs(v))
    meanabs_v = sum(abs(v)) / real(Np, rp)

    call log_val("[INFO ]", "mean (Pc(tf)/Pc(ti)-1)/(tf-ti)" ,    mean_v)
    call log_val("[INFO ]", "rms  (Pc(tf)/Pc(ti)-1)/(tf-ti)",    rms_v)
    call log_val("[INFO ]", "std  (Pc(tf)/Pc(ti)-1)/(tf-ti)",    std_v)
    call log_val("[INFO ]", "max |(Pc(tf)/Pc(ti)-1)/(tf-ti)|",  maxabs_v)
    call log_val("[INFO ]", "mean|(Pc(tf)/Pc(ti)-1)/(tf-ti)|",   meanabs_v)

    imax = maxloc(abs(v))
    write(TMP,'(a,1x,a,":",1x,'//FINT//')') "[INFO ]", pad_label("argmax |(Pc(tf)/Pc(ti)-1)/(tf-ti)| (particle index)"), imax(1)
    call log_raw(TMP)

    if (maxabs_v > 1.0e-1_rp) then
      no_errors = no_errors + 1
      call log_alert("Energy proxy: large max |(Pc(tf)/Pc(ti)-1)/(tf-ti)| (threshold 1e-1)")
    end if

    deallocate(v)
  end subroutine print_dlnPc_stats


  !=============================================================================
  ! Gaussianity summary over full 2D field (Nt x Np)
  !=============================================================================
  subroutine print_gaussianity_block(name, A, Np, Nt, no_errors, compare)
    character(len=*), intent(in)    :: name
    real(rp),         intent(in)    :: A(:,:), compare
    integer,          intent(in)    :: Np, Nt
    integer,          intent(inout) :: no_errors

    real(rp) :: norm, m1, m2, m4, kurt3
    real(rp), parameter :: tol_kurt3 = 0.15_rp

    norm  = 1.0_rp / (real(Np, rp) * real(Nt, rp))
    m1    = sum(A)      * norm
    m2    = sum(A**2)   * norm
    m4    = sum(A**4)   * norm
    kurt3 = m4 / max(m2*m2, tiny(1.0_rp)) / 3.0_rp

    call log_gauss_line(name, m1, m2, kurt3, tol_kurt3, no_errors, compare)
  end subroutine print_gaussianity_block


  !=============================================================================
  ! Gaussianity summary for one time slice k (1..Nt): uses A(k,:)
  !=============================================================================
  subroutine print_gaussianity_slice(name, A, k, Np, no_errors, compare)
    character(len=*), intent(in)    :: name
    real(rp),         intent(in)    :: A(:,:), compare
    integer,          intent(in)    :: k, Np
    integer,          intent(inout) :: no_errors

    real(rp) :: norm, m1, m2, m4, kurt3
    real(rp), parameter :: tol_kurt3 = 0.15_rp

    norm  = 1.0_rp / real(Np, rp)
    m1    = sum(A(k, :))      * norm
    m2    = sum(A(k, :)**2)   * norm
    m4    = sum(A(k, :)**4)   * norm
    kurt3 = m4 / max(m2*m2, tiny(1.0_rp)) / 3.0_rp

    call log_gauss_line(name, m1, m2, kurt3, tol_kurt3, no_errors, compare)
  end subroutine print_gaussianity_slice


  subroutine log_gauss_line(name, mean, ex2, kurt3, tol_kurt3, no_errors, compare)
    character(len=*), intent(in)    :: name
    real(rp),         intent(in)    :: mean, ex2, kurt3, tol_kurt3, compare
    integer,          intent(inout) :: no_errors

    character(len=WTAG) :: tag
    logical :: bad

    bad = abs(kurt3 - 1.0_rp) > tol_kurt3
    if (bad) then
      tag = "[ALERT]"
      no_errors = no_errors + 1
      call alerts_add("Gaussianity: " // trim(name) // " kurt/3 out of tol")
    else
      tag = "[INFO ]"
    end if

    write(TMP,'(a,1x,a,": E[x]=",1x,'//FVAL//',", E[x^2]=",1x,'//FVAL//', &
              ", kurt/3=",1x,'//FVAL//',", tol=",1x,'//FVAL//',", compared to = ",1x,'//FVAL//')') &
         tag, pad_label(trim(name)), mean, ex2, kurt3, tol_kurt3, compare
    call log_raw(TMP)
  end subroutine log_gauss_line


  !=============================================================================
  ! Helper: check ratio against tol with uniform messages
  !=============================================================================
  subroutine check_ratio(label, ratio, tol, no_errors)
    character(len=*), intent(in)    :: label
    real(rp),         intent(in)    :: ratio, tol
    integer,          intent(inout) :: no_errors

    if (ratio > tol) then
      no_errors = no_errors + 1
      call log_ratio("[ALERT]", trim(label), ratio, tol)
    else
      call log_ratio("[  OK ]", trim(label), ratio, tol)
    end if
  end subroutine check_ratio


  !=============================================================================
  ! Logging utilities (centralized formatting)
  !=============================================================================
  subroutine log_ok(msg)
    character(len=*), intent(in) :: msg
    call log_line("[  OK ]", msg)
  end subroutine log_ok

  subroutine log_alert(msg)
    character(len=*), intent(in) :: msg
    call alerts_add(trim(msg))
    call log_line("[ALERT]", msg)
  end subroutine log_alert

  subroutine log_warn(msg)
    character(len=*), intent(in) :: msg
    call log_line("[WARN]", msg)
  end subroutine log_warn

  subroutine log_info(msg)
    character(len=*), intent(in) :: msg
    call log_line("[INFO ]", msg)
  end subroutine log_info

  subroutine log_line(tag, msg)
    use, intrinsic :: iso_fortran_env, only : output_unit
    character(len=*), intent(in) :: tag, msg

    write(output_unit,'(a,1x,a)') tag, trim(msg)
    if (LOG_TO_FILE) write(LOGU,'(a,1x,a)') tag, trim(msg)
  end subroutine log_line

  subroutine log_val(tag, name, val)
    character(len=*), intent(in) :: tag, name
    real(rp),         intent(in) :: val
    write(TMP,'(a,1x,a,":",1x,'//FVAL//')') tag, pad_label(name), val
    call log_raw(TMP)
  end subroutine log_val

  subroutine log_vals3(tag, name, v1, v2, v3)
    character(len=*), intent(in) :: tag, name
    real(rp),         intent(in) :: v1, v2, v3
    write(TMP,'(a,1x,a,":",3(1x,'//FVAL//'))') tag, pad_label(name), v1, v2, v3
    call log_raw(TMP)
  end subroutine log_vals3

  subroutine log_ratio(tag, name, ratio, tol)
    character(len=*), intent(in) :: tag, name
    real(rp),         intent(in) :: ratio, tol
    write(TMP,'(a,1x,a,": ratio=",1x,'//FRAT//',", tol=",1x,'//FRAT//')') tag, pad_label(name), ratio, tol
    call log_raw(TMP)
    if (trim(tag) == "[ALERT]") call alerts_add(trim(name) // " ratio exceeded tol")
  end subroutine log_ratio

  subroutine log_count(label, nbad, n)
    character(len=*), intent(in) :: label
    integer,          intent(in) :: nbad, n
    write(TMP,'(a,1x,a,":",1x,'//FINT//',a,1x,'//FINT//')') "[INFO ]", pad_label(label), nbad, " /", n
    call log_raw(TMP)
  end subroutine log_count


  !=============================================================================
  ! Label padding
  !=============================================================================
  function pad_label(s) result(out)
    character(len=*), intent(in) :: s
    character(len=WLBL)          :: out
    integer :: ls

    out = " "
    ls  = len_trim(s)

    if (ls >= WLBL) then
      out = s(1:WLBL)
    else
      out(1:ls) = s(1:ls)
    end if
  end function pad_label


  !=============================================================================
  ! Alert collection + recap
  !=============================================================================
  subroutine alerts_reset()
    integer :: i
    n_alerts = 0
    do i = 1, MAX_ALERTS
      alert_msgs(i) = " "
    end do
  end subroutine alerts_reset

  subroutine alerts_add(msg)
    character(len=*), intent(in) :: msg
    if (n_alerts < MAX_ALERTS) then
      n_alerts = n_alerts + 1
      alert_msgs(n_alerts) = trim(msg)
    end if
  end subroutine alerts_add

  subroutine alerts_recap()
    integer :: i
    if (n_alerts <= 0) return

    call log_raw(" ")
    call log_info("Alerts recap (most recent run)")
    do i = 1, n_alerts
      call log_line("[ALERT]", trim(alert_msgs(i)))
    end do
  end subroutine alerts_recap


  !=============================================================================
  ! Small helper for unused dummy variable
  !=============================================================================
  subroutine ignore_unused(x)
    real(rp), intent(in) :: x
    if (x /= x) then
      call log_warn("ignore_unused triggered (NaN)")
    end if
  end subroutine ignore_unused


  subroutine log_open(filename, append)
    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
    character(len=*), intent(in)           :: filename
    logical,          intent(in), optional :: append
    logical :: do_append
    integer :: ios

    do_append = .false.
    if (present(append)) do_append = append

    LOGFILE = trim(filename)

    if (LOG_TO_FILE) call log_close()

    if (do_append) then
      open(newunit=LOGU, file=trim(LOGFILE), status="unknown", position="append", &
           action="write", iostat=ios)
    else
      open(newunit=LOGU, file=trim(LOGFILE), status="replace", action="write", iostat=ios)
    end if

    if (ios /= 0) then
      LOG_TO_FILE = .false.
      write(error_unit,'(a,1x,a,1x,i0)') "[WARN]", "Could not open log file:", trim(LOGFILE), ios
      LOGU = -1
    else
      LOG_TO_FILE = .true.
      write(output_unit,'(a,1x,a)') "[INFO ]", "Logging to file: " // trim(LOGFILE)
      write(LOGU,'(a,1x,a)') "[INFO ]", "Logging started"
    end if
  end subroutine log_open


  subroutine log_close()
    use, intrinsic :: iso_fortran_env, only : error_unit
    integer :: ios

    if (.not. LOG_TO_FILE) return

    write(LOGU,'(a)') " "
    write(LOGU,'(a,1x,a)') "[INFO ]", "Logging ended"
    close(LOGU, iostat=ios)

    if (ios /= 0) then
      write(error_unit,'(a,1x,a,1x,i0)') "[WARN]", "Could not close log file cleanly, iostat=", ios
    end if

    LOG_TO_FILE = .false.
    LOGU = -1
  end subroutine log_close


  subroutine log_raw(line)
    use, intrinsic :: iso_fortran_env, only : output_unit
    character(len=*), intent(in) :: line
    write(output_unit,'(a)') trim(line)
    if (LOG_TO_FILE) write(LOGU,'(a)') trim(line)
  end subroutine log_raw

end module testing_T3ST_mod
