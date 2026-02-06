module testing_T3ST_mod
  use constants
  use random_numbers              ! provides kx,ky,kz,w,ph,L and kxs,kys,kzs,ws,phs,Ls and USE_real
  use drift_kernels
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
  implicit none
  private

  public :: testing_T3ST

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
  !=============================================================================
  ! Main entry point
  !=============================================================================
  subroutine testing_T3ST(X, Y, Z, Vp, mu, pb, Vstar1, Vstar2, Vstar3, no_errors)
    real(dp), intent(in)  :: X(:), Y(:), Z(:), Vp(:), mu(:), pb(:)
    real(dp), intent(in)  :: Vstar1, Vstar2, Vstar3
    integer,  intent(out) :: no_errors

    integer  :: Np_partial, Nt_partial, np_avail
    real(dp) :: dt, ratio32

    real(dp), allocatable :: Xh(:,:), Yh(:,:), Zh(:,:), Vph(:,:), muh(:,:)
    real(dp), allocatable :: Hh(:,:), Pch(:,:), c1h(:,:), c2h(:,:), c3h(:,:),  Vtx(:,:),  Vty(:,:)
    real(dp), allocatable :: Xaux(:), Yaux(:), Zaux(:), Vpaux(:), muaux(:)

    integer  :: nstep_eff
    real(dp) :: T, compare
    real(dp), parameter :: epsH = 1.0e-300_dp

    no_errors = 0
    ! -------------------------
    ! Controls (hard-coded here)
    ! -------------------------
    Np_partial = 16000
    Nt_partial = 50

    np_avail   = size(X)
    Np_partial = min(Np_partial, np_avail)
    dt         = real(tmax / Nt, dp)   ! base dt used by rk4_propagation below

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
             Vty(1:Nt_partial, 1:Np_partial))

    allocate(Xaux(Np), Yaux(Np), Zaux(Np), Vpaux(Np), muaux(Np))
    Xaux=X; Yaux=Y; Zaux=Z; Vpaux=Vp; muaux=mu;

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
    call check_nan("pb",     pb,     no_errors)
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
    call log_vals3("[INFO]", "Vstar1 Vstar2 Vstar3", Vstar1, Vstar2, Vstar3)

    if (abs(Vstar1) > 1.0e-9_dp) then
      no_errors = no_errors + 1
      call log_alert("Vstar1 not ~ 0 (radial drift present)")
    else
      call log_ok("Vstar1 ~ 0 (radial drift OK)")
    end if

    if (abs(Vstar2) <= tiny(1.0_dp)) then
      no_errors = no_errors + 1
      call log_alert("Vstar2 ~ 0 (cannot evaluate Vstar3/Vstar2)")
    else
      ratio32 = Vstar3 / Vstar2
      call log_ratio("[INFO]", "Vstar3/Vstar2", abs(ratio32), 0.1_dp)

      if (abs(ratio32) > 0.1_dp) then
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
    call check_inputs(X, Y, Z, Vp, mu, pb, a0, 1.0e-12_dp, no_errors)

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
    call log_val("[INFO]", "Np_partial", real(Np_partial, dp))
    call log_val("[INFO]", "Nt_partial", real(Nt_partial, dp))
    call log_val("[INFO]", "dt", dt)

    call rk4_propagation(X, Y, Z, Vp, mu, pb, dt, Nt_partial, Np_partial, &
                         Xh, Yh, Zh, Vph, muh, Hh, Pch, c1h, c2h, c3h, Vtx, Vty)

    nstep_eff = max(1, Nt_partial - 1)
    T         = dt * real(nstep_eff, dp)

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
    compare = 1.0_dp
    call print_gaussianity_block("phi",  c1h, Np_partial, Nt_partial, no_errors, compare)
    call print_gaussianity_block("phix", c2h, Np_partial, Nt_partial, no_errors, compare)
    call print_gaussianity_block("phiy", c3h, Np_partial, Nt_partial, no_errors, compare)

    call section("Gaussianity drivers (qualitative)")
    call log_info("Larmor & Frequencies -> NO effects on Gaussianity of fields")
    call log_info("Balloning breaks phix, phiz Gaussianity as I_0[4/Lbalon^2]/I_0[2/Lbalon^2]^2/2")
    call log_info("Balloning breaks amplitude as: I_0(1/lbalonz^2)*exp(-lbalonz^2)")
    call log_info("Tilting breaks phix/phiz Gaussianity; k_x^eff = kx+ky a --> ??")

    ! =====================================================================
    ! Gaussianity / moment checks: time slice (k = 1)
    ! =====================================================================
    call section("Gaussianity summary: time slice k=1 (kurt/3 = 1 for Gaussian)")
    call print_gaussianity_slice("phi",  c1h, 1, Np_partial, no_errors, compare)
    call print_gaussianity_slice("phix", c2h, 1, Np_partial, no_errors, compare)
    call print_gaussianity_slice("phiy", c3h, 1, Np_partial, no_errors, compare)

    call rk4_propagation(Xaux, Yaux, Zaux, Vpaux, muaux, pb, dt, 1, Np_partial, &
                         Xh, Yh, Zh, Vph, muh, Hh, Pch, c1h, c2h, c3h, Vtx, Vty)

    if (USE_real.eq.OFF) then
      compare = ((R0/rhoi*Phi)**2)*sum(ky**2)/real(Np*Nc,dp)
    elseif (USE_real.eq.ON) then
      compare = ((R0/rhoi*Phi)**2)*sum(ky**2)/real(Np*Nc,dp)
    endif
    if (USE_balloon.eq.ON) then
      compare = compare**exp(-1.0_dp/lbalonz**2)*(1 + 1.0_dp/lbalonz**2/4.0_dp + 1.0_dp/lbalonz**4/64.0_dp)
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
      '  E[Vnx] =', sum(Vty(1,:)-Vtx(1,:))/Np_partial, ' E[Vnx^2] =', sum((Vty(1,:)-Vtx(1,:))**2)/Np_partial,  ' compared to (only if theta=pi/2)', compare
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
  subroutine rk4_propagation(X, Y, Z, Vp, mu, pb, dt, Nt_partial, Np_partial, &
                             Xh, Yh, Zh, Vph, muh, Hh, Pch, c1h, c2h, c3h, Vtxall, Vtyall)

    real(dp), intent(in) :: X(:), Y(:), Z(:), Vp(:), mu(:), pb(:)
    integer,  intent(in) :: Np_partial, Nt_partial
    real(dp), intent(in) :: dt

    real(dp), intent(out) :: Xh(:,:), Yh(:,:), Zh(:,:), Vph(:,:), muh(:,:)
    real(dp), intent(out) :: Hh(:,:), Pch(:,:), c1h(:,:), c2h(:,:), c3h(:,:), Vtxall(:,:), Vtyall(:,:)

    real(dp), pointer, contiguous :: Qx(:), Qy(:), Qz(:), Qw(:), Qph(:), QL(:), Q0wrp(:)

    real(dp) :: xi, yi, zi, mui, vpi, pbi
    real(dp) :: vx, vy, vz, vm, ap
    real(dp) :: Wx, Wy, Wz, Wm, Wp
    real(dp) :: q1, q2, q3
    real(dp) :: Hi, B, Vtx, Vty, Pc
    real(dp) :: check_1, check_2, check_3, time
    integer  :: ias, k

    do ias = 1, Np_partial

      ! Select per-particle or single-spectrum arrays
      if (USE_real == OFF) then
        Qx    => kx(:, ias)
        Qy    => ky(:, ias)
        Qz    => kz(:, ias)
        Qw    => w(:,  ias)
        Qph   => ph(:, ias)
        QL    => L(:,  ias)
        Q0wrp => q00wrap(:, ias)
      else
        Qx    => kxs
        Qy    => kys
        Qz    => kzs
        Qw    => ws
        Qph   => phs
        QL    => Ls(:, ias)
        Q0wrp => q00wraps
      end if

      xi  = X(ias)
      yi  = Y(ias)
      zi  = Z(ias)
      mui = mu(ias)
      vpi = Vp(ias)
      pbi = pb(ias)
      call ignore_unused(pbi)

      do k = 1, Nt_partial

        time = dt * real(k - 1, dp)

        call Drift2(dt, xi, yi, zi, vpi, mui, q1, q2, q3, time,                 &
                    vx, vy, vz, ap, vm, Hi, Pc, B, Vtx, Vty,                    &
                    check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph, QL, Q0wrp)

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
        Vtyall(k,ias) = Vty

        Wx = vx / 6.0_dp
        Wy = vy / 6.0_dp
        Wz = vz / 6.0_dp
        Wm = vm / 6.0_dp
        Wp = ap / 6.0_dp

        call Drift2(dt, xi + vx*dt/2.0_dp, yi + vy*dt/2.0_dp, zi + vz*dt/2.0_dp, &
                    vpi + ap*dt/2.0_dp, mui + vm*dt/2.0_dp, q1, q2, q3,         &
                    time + dt/2.0_dp,                                           &
                    vx, vy, vz, ap, vm, Hi, Pc, B, Vtx, Vty,                    &
                    check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph, QL, Q0wrp)

        Wx = Wx + vx / 3.0_dp
        Wy = Wy + vy / 3.0_dp
        Wz = Wz + vz / 3.0_dp
        Wm = Wm + vm / 3.0_dp
        Wp = Wp + ap / 3.0_dp

        call Drift2(dt, xi + vx*dt/2.0_dp, yi + vy*dt/2.0_dp, zi + vz*dt/2.0_dp, &
                    vpi + ap*dt/2.0_dp, mui + vm*dt/2.0_dp, q1, q2, q3,         &
                    time + dt/2.0_dp,                                           &
                    vx, vy, vz, ap, vm, Hi, Pc, B, Vtx, Vty,                    &
                    check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph, QL, Q0wrp)

        Wx = Wx + vx / 3.0_dp
        Wy = Wy + vy / 3.0_dp
        Wz = Wz + vz / 3.0_dp
        Wm = Wm + vm / 3.0_dp
        Wp = Wp + ap / 3.0_dp

        call Drift2(dt, xi + vx*dt, yi + vy*dt, zi + vz*dt,                      &
                    vpi + ap*dt, mui + vm*dt, q1, q2, q3, time + dt,            &
                    vx, vy, vz, ap, vm, Hi, Pc, B, Vtx, Vty,                    &
                    check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph, QL, Q0wrp)

        Wx = Wx + vx / 6.0_dp
        Wy = Wy + vy / 6.0_dp
        Wz = Wz + vz / 6.0_dp
        Wm = Wm + vm / 6.0_dp
        Wp = Wp + ap / 6.0_dp

        xi  = xi  + dt * Wx
        yi  = yi  + dt * Wy
        zi  = zi  + dt * Wz
        vpi = vpi + dt * Wp
        mui = mui + dt * Wm

      end do
    end do

  end subroutine rk4_propagation


  !=============================================================================
  ! Unified NaN checker (assumed-rank)
  !=============================================================================
  subroutine check_nan(name, x, no_errors)
    character(len=*), intent(in)    :: name
    real(dp),         intent(in)    :: x(..)
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
    call log_line("[INFO]", repeat("-", 1) // " " // trim(title) // " " // repeat("-", max(0, WSEC - 9 - len_trim(title))))
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
    write(tmp2,'(a,1x,'//FINT//')') "[INFO] Total number of alerts:", no_errors
    call log_raw(tmp2)
    call log_raw(repeat("=", WSEC))
  end subroutine banner_end


  !=============================================================================
  ! Input checks (uniform messages + formatted numbers)
  !=============================================================================
  subroutine check_inputs(X, Y, Z, Vp, mu, pb, a1, tol, no_errors)
    real(dp), intent(in)    :: X(:), Y(:), Z(:), Vp(:), mu(:), pb(:)
    real(dp), intent(in)    :: a1, tol
    integer,  intent(inout) :: no_errors

    integer  :: n, nbad
    real(dp) :: pii, mean_pb

    n   = size(X)
    pii = 3.14159265358979323846_dp

    if (any([ size(Y), size(Z), size(Vp), size(mu), size(pb) ] /= n)) then
      no_errors = no_errors + 1
      call log_alert("Input arrays have inconsistent sizes")
      write(TMP,'(a,1x,a,": X=",1x,'//FINT//',", Y=",1x,'//FINT//',", Z=",1x,'//FINT//', &
                   ", Vp=",1x,'//FINT//',", mu=",1x,'//FINT//',", pb=",1x,'//FINT//')') &
        "[INFO]", pad_label("sizes"), size(X), size(Y), size(Z), size(Vp), size(mu), size(pb)
      call log_raw(TMP)
      return
    else
      call log_ok("All input arrays have consistent size")
      write(TMP,'(a,1x,a,":",1x,'//FINT//')') "[INFO]", pad_label("n"), n
      call log_raw(TMP)
    end if

    nbad = count(X <= 0.0_dp)
    if (nbad > 0) then
      no_errors = no_errors + 1
      call log_alert("Condition failed: X must be > 0")
      call log_count("violations", nbad, n)
    else
      call log_ok("X > 0 everywhere")
    end if

    nbad = count(abs(Y) >= 1.0_dp)
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

    nbad = count(mu <= 0.0_dp)
    if (nbad > 0) then
      no_errors = no_errors + 1
      call log_alert("Condition failed: mu must be > 0")
      call log_count("violations", nbad, n)
    else
      call log_ok("mu > 0 everywhere")
    end if

    nbad = count(pb <= 0.0_dp)
    if (nbad > 0) then
      no_errors = no_errors + 1
      call log_alert("Condition failed: pb must be > 0")
      call log_count("violations", nbad, n)
    else
      call log_ok("pb > 0 everywhere")
    end if

    mean_pb = sum(pb) / real(n, dp)
    if (abs(mean_pb - 1.0_dp) > tol) then
      no_errors = no_errors + 1
      call log_alert("Condition failed: mean(pb) must be 1")
      call log_val("[INFO]", "mean(pb)", mean_pb)
      call log_val("[INFO]", "tol",      tol)
    else
      call log_ok("mean(pb) = 1 within tolerance")
      call log_val("[INFO]", "mean(pb)", mean_pb)
      call log_val("[INFO]", "tol",      tol)
    end if

    nbad = count(sqrt((X - 1.0_dp)**2 + Y**2) >= a1)
    if (nbad > 0) then
      no_errors = no_errors + 1
      call log_alert("Condition failed: sqrt((X-1)^2 + Y^2) < a0")
      call log_count("violations", nbad, n)
      call log_val("[INFO]", "a0", a1)
    else
      call log_ok("sqrt((X-1)^2 + Y^2) < a0 everywhere")
      call log_val("[INFO]", "a0", a1)
    end if

  end subroutine check_inputs


  !=============================================================================
  ! Wave checks (rank-agnostic) WITHOUT L/Ls
  !=============================================================================
  subroutine check_waves_any(kxA, kyA, kzA, wA, phA, Vstar2, Vstar3, no_errors)
    real(dp), intent(in)    :: kxA(..), kyA(..), kzA(..), wA(..), phA(..)
    real(dp), intent(in)    :: Vstar2, Vstar3
    integer,  intent(inout) :: no_errors

    real(dp), parameter :: tiny_dp = 1.0e-12_dp
    real(dp), parameter :: eps_ref = 1.0e-30_dp

    real(dp) :: m1(5), m2(5), m4(5)
    real(dp) :: tol_mean(5), tol_m2(5), tol_m4(5)
    real(dp) :: m2_ref(5), m4_ref(5)
    real(dp) :: ratio, pii
    real(dp) :: t_m2_ky, t_m2_kz, t_m4_ky, t_m4_kz
    real(dp) :: fac

    call log_info("Wave diagnostics (moments / target checks)")

    pii = 3.14159265358979323846_dp

    tol_mean = 0.05_dp
    tol_m2   = 0.05_dp
    tol_m4   = 0.10_dp

    m2_ref(1) = 1.0_dp / (lambdax*lambdax)
    m2_ref(2) = 3.0_dp / (lambday*lambday) + k0i*k0i
    m2_ref(3) = 1.0_dp / (lambdaz*lambdaz)
    m2_ref(5) = (pii*pii) / 3.0_dp

    m4_ref(1) = 3.0_dp * (m2_ref(1)**2)
    m4_ref(2) = (15.0_dp/(lambday**4) + 10.0_dp*k0i**2/(lambday**2) + k0i**4)
    m4_ref(3) = 3.0_dp * (m2_ref(3)**2)
    m4_ref(5) = (pii**4) / 5.0_dp

    t_m2_ky = m2_ref(2)
    t_m2_kz = m2_ref(3)
    t_m4_ky = m4_ref(2)
    t_m4_kz = m4_ref(3)

    fac = (Ln*(rhoi/R0)**2)

    m2_ref(4) = 1.0_dp/(tauc*tauc) + (fac**2) * (Vstar2**2*t_m2_ky + Vstar3**2*t_m2_kz)

    m4_ref(4) = 3.0_dp/(tauc**4)                                                         &
              + (fac**4) * (Vstar2**4*t_m4_ky + Vstar3**4*t_m4_kz                         &
              + 2.0_dp*Vstar2**2*Vstar3**2*t_m2_kz*t_m2_ky)                               &
              + 2.0_dp/(tauc**2) * (fac**2) * (Vstar2**2*t_m2_ky + Vstar3**2*t_m2_kz)

    call moments_ar(kxA, m1(1), m2(1), m4(1))
    call moments_ar(kyA, m1(2), m2(2), m4(2))
    call moments_ar(kzA, m1(3), m2(3), m4(3))
    call moments_ar(wA,  m1(4), m2(4), m4(4))
    call moments_ar(phA, m1(5), m2(5), m4(5))

    ratio = abs(m1(1)) / max(sqrt(m2(1)), tiny_dp)
    call check_ratio("kx mean ~ 0", ratio, tol_mean(1), no_errors)

    ratio = abs(m1(2)) / max(sqrt(m2(2)), tiny_dp)
    call check_ratio("ky mean ~ 0", ratio, tol_mean(2), no_errors)

    ratio = abs(m1(3)) / max(sqrt(m2(3)), tiny_dp)
    call check_ratio("kz mean ~ 0", ratio, tol_mean(3), no_errors)

    ratio = abs(m1(4)) / max(sqrt(m2(4)), tiny_dp)
    call check_ratio("w  mean ~ 0", ratio, tol_mean(4), no_errors)

    ratio = abs(m1(5)) / max(sqrt(m2(5)), tiny_dp)
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
    real(dp),         intent(in)    :: val, ref, tol, eps_ref
    integer,          intent(inout) :: no_errors
    real(dp) :: relerr

    relerr = abs(val - ref) / max(abs(ref), eps_ref)

    if (relerr > tol) then
      no_errors = no_errors + 1
      call log_moment("[ALERT]", label, val, ref, relerr, tol)
    else
      call log_moment("[ OK ]", label, val, ref, relerr, tol)
    end if
  end subroutine check_moment


  subroutine log_moment(tag, label, val, ref, relerr, tol)
    character(len=*), intent(in) :: tag, label
    real(dp),         intent(in) :: val, ref, relerr, tol

    write(TMP,'(a,1x,a,": val=",1x,'//FVAL//',", ref=",1x,'//FVAL//',", rel=",1x,'//FRAT//',", tol=",1x,'//FRAT//')') &
      tag, pad_label(label), val, ref, relerr, tol
    call log_raw(TMP)

    if (trim(tag) == "[ALERT]") call alerts_add(trim(label))
  end subroutine log_moment


  !=============================================================================
  ! Moments for assumed-rank arrays (rank 0/1/2 supported)
  !=============================================================================
  subroutine moments_ar(x, m1, m2, m4)
    real(dp), intent(in)  :: x(..)
    real(dp), intent(out) :: m1, m2, m4

    integer  :: n
    real(dp) :: rn

    select rank (x)
    rank (0)
      n  = 1
      rn = 1.0_dp
      m1 = x
      m2 = x*x
      m4 = (x*x)*(x*x)

    rank (1)
      n  = size(x)
      rn = real(n, dp)
      m1 = sum(x)    / rn
      m2 = sum(x**2) / rn
      m4 = sum(x**4) / rn

    rank (2)
      n  = size(x)
      rn = real(n, dp)
      m1 = sum(x)    / rn
      m2 = sum(x**2) / rn
      m4 = sum(x**4) / rn

    rank default
      call log_warn("moments_ar: unsupported rank encountered")
      m1 = huge(1.0_dp)
      m2 = huge(1.0_dp)
      m4 = huge(1.0_dp)
    end select
  end subroutine moments_ar


  !=============================================================================
  ! Energy proxy: d(ln H)/dt finite difference stats
  !=============================================================================
  subroutine print_dlnH_stats(Hh, Np, Nt, T, epsH, no_errors)
    real(dp), intent(in)    :: Hh(:,:)
    integer,  intent(in)    :: Np, Nt
    real(dp), intent(in)    :: T, epsH
    integer,  intent(inout) :: no_errors

    real(dp), allocatable :: v(:)
    real(dp) :: mean_v, rms_v, std_v, maxabs_v, meanabs_v
    integer  :: imax(1)

    allocate(v(Np))

    v = (Hh(Nt, :) - Hh(1, :)) / max(abs(Hh(1, :)), epsH) / T

    mean_v    = sum(v) / real(Np, dp)
    rms_v     = sqrt(sum(v**2) / real(Np, dp))
    std_v     = sqrt(max(0.0_dp, rms_v**2 - mean_v**2))
    maxabs_v  = maxval(abs(v))
    meanabs_v = sum(abs(v)) / real(Np, dp)

    call log_val("[INFO]", "mean (H(tf)/H(ti)-1)/(tf-ti)",    mean_v)
    call log_val("[INFO]", "rms  (H(tf)/H(ti)-1)/(tf-ti)",    rms_v)
    call log_val("[INFO]", "std  (H(tf)/H(ti)-1)/(tf-ti)",    std_v)
    call log_val("[INFO]", "max |(H(tf)/H(ti)-1)/(tf-ti)|",   maxabs_v)
    call log_val("[INFO]", "mean|(H(tf)/H(ti)-1)/(tf-ti)|",   meanabs_v)

    imax = maxloc(abs(v))
    write(TMP,'(a,1x,a,":",1x,'//FINT//')') "[INFO]", pad_label("argmax |(H(tf)/H(ti)-1)/(tf-ti)| (particle index)"), imax(1)
    call log_raw(TMP)

    if (maxabs_v > 1.0e-1_dp) then
      no_errors = no_errors + 1
      call log_alert("Energy proxy: large max |(H(tf)/H(ti)-1)/(tf-ti)| (threshold 1e-1)")
    end if

    deallocate(v)
  end subroutine print_dlnH_stats


  !=============================================================================
  ! Energy proxy: d(ln H)/dt finite difference stats
  !=============================================================================
  subroutine print_dlnPc_stats(Hh, Np, Nt, T, epsH, no_errors)
    real(dp), intent(in)    :: Hh(:,:)
    integer,  intent(in)    :: Np, Nt
    real(dp), intent(in)    :: T, epsH
    integer,  intent(inout) :: no_errors

    real(dp), allocatable :: v(:)
    real(dp) :: mean_v, rms_v, std_v, maxabs_v, meanabs_v
    integer  :: imax(1)

    allocate(v(Np))

    v = (Hh(Nt, :) - Hh(1, :)) / max(abs(Hh(1, :)), epsH) / T

    mean_v    = sum(v) / real(Np, dp)
    rms_v     = sqrt(sum(v**2) / real(Np, dp))
    std_v     = sqrt(max(0.0_dp, rms_v**2 - mean_v**2))
    maxabs_v  = maxval(abs(v))
    meanabs_v = sum(abs(v)) / real(Np, dp)

    call log_val("[INFO]", "mean (Pc(tf)/Pc(ti)-1)/(tf-ti)" ,    mean_v)
    call log_val("[INFO]", "rms  (Pc(tf)/Pc(ti)-1)/(tf-ti)",    rms_v)
    call log_val("[INFO]", "std  (Pc(tf)/Pc(ti)-1)/(tf-ti)",    std_v)
    call log_val("[INFO]", "max |(Pc(tf)/Pc(ti)-1)/(tf-ti)|",  maxabs_v)
    call log_val("[INFO]", "mean|(Pc(tf)/Pc(ti)-1)/(tf-ti)|",   meanabs_v)

    imax = maxloc(abs(v))
    write(TMP,'(a,1x,a,":",1x,'//FINT//')') "[INFO]", pad_label("argmax |(Pc(tf)/Pc(ti)-1)/(tf-ti)| (particle index)"), imax(1)
    call log_raw(TMP)

    if (maxabs_v > 1.0e-1_dp) then
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
    real(dp),         intent(in)    :: A(:,:), compare
    integer,          intent(in)    :: Np, Nt
    integer,          intent(inout) :: no_errors

    real(dp) :: norm, m1, m2, m4, kurt3
    real(dp), parameter :: tol_kurt3 = 0.15_dp

    norm  = 1.0_dp / (real(Np, dp) * real(Nt, dp))
    m1    = sum(A)      * norm
    m2    = sum(A**2)   * norm
    m4    = sum(A**4)   * norm
    kurt3 = m4 / max(m2*m2, tiny(1.0_dp)) / 3.0_dp

    call log_gauss_line(name, m1, m2, kurt3, tol_kurt3, no_errors, compare)
  end subroutine print_gaussianity_block


  !=============================================================================
  ! Gaussianity summary for one time slice k (1..Nt): uses A(k,:)
  !=============================================================================
  subroutine print_gaussianity_slice(name, A, k, Np, no_errors, compare)
    character(len=*), intent(in)    :: name
    real(dp),         intent(in)    :: A(:,:), compare
    integer,          intent(in)    :: k, Np
    integer,          intent(inout) :: no_errors

    real(dp) :: norm, m1, m2, m4, kurt3
    real(dp), parameter :: tol_kurt3 = 0.15_dp

    norm  = 1.0_dp / real(Np, dp)
    m1    = sum(A(k, :))      * norm
    m2    = sum(A(k, :)**2)   * norm
    m4    = sum(A(k, :)**4)   * norm
    kurt3 = m4 / max(m2*m2, tiny(1.0_dp)) / 3.0_dp

    call log_gauss_line(name, m1, m2, kurt3, tol_kurt3, no_errors, compare)
  end subroutine print_gaussianity_slice


  subroutine log_gauss_line(name, mean, ex2, kurt3, tol_kurt3, no_errors, compare)
    character(len=*), intent(in)    :: name
    real(dp),         intent(in)    :: mean, ex2, kurt3, tol_kurt3, compare
    integer,          intent(inout) :: no_errors

    character(len=WTAG) :: tag
    logical :: bad

    bad = abs(kurt3 - 1.0_dp) > tol_kurt3
    if (bad) then
      tag = "[ALERT]"
      no_errors = no_errors + 1
      call alerts_add("Gaussianity: " // trim(name) // " kurt/3 out of tol")
    else
      tag = "[INFO]"
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
    real(dp),         intent(in)    :: ratio, tol
    integer,          intent(inout) :: no_errors

    if (ratio > tol) then
      no_errors = no_errors + 1
      call log_ratio("[ALERT]", trim(label), ratio, tol)
    else
      call log_ratio("[ OK ]", trim(label), ratio, tol)
    end if
  end subroutine check_ratio


  !=============================================================================
  ! Logging utilities (centralized formatting)
  !=============================================================================
  subroutine log_ok(msg)
    character(len=*), intent(in) :: msg
    call log_line("[ OK ]", msg)
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
    call log_line("[INFO]", msg)
  end subroutine log_info

  subroutine log_line(tag, msg)
    use, intrinsic :: iso_fortran_env, only : output_unit
    character(len=*), intent(in) :: tag, msg

    write(output_unit,'(a,1x,a)') tag, trim(msg)
    if (LOG_TO_FILE) write(LOGU,'(a,1x,a)') tag, trim(msg)
  end subroutine log_line

  subroutine log_val(tag, name, val)
    character(len=*), intent(in) :: tag, name
    real(dp),         intent(in) :: val
    write(TMP,'(a,1x,a,":",1x,'//FVAL//')') tag, pad_label(name), val
    call log_raw(TMP)
  end subroutine log_val

  subroutine log_vals3(tag, name, v1, v2, v3)
    character(len=*), intent(in) :: tag, name
    real(dp),         intent(in) :: v1, v2, v3
    write(TMP,'(a,1x,a,":",3(1x,'//FVAL//'))') tag, pad_label(name), v1, v2, v3
    call log_raw(TMP)
  end subroutine log_vals3

  subroutine log_ratio(tag, name, ratio, tol)
    character(len=*), intent(in) :: tag, name
    real(dp),         intent(in) :: ratio, tol
    write(TMP,'(a,1x,a,": ratio=",1x,'//FRAT//',", tol=",1x,'//FRAT//')') tag, pad_label(name), ratio, tol
    call log_raw(TMP)
    if (trim(tag) == "[ALERT]") call alerts_add(trim(name) // " ratio exceeded tol")
  end subroutine log_ratio

  subroutine log_count(label, nbad, n)
    character(len=*), intent(in) :: label
    integer,          intent(in) :: nbad, n
    write(TMP,'(a,1x,a,":",1x,'//FINT//',a,1x,'//FINT//')') "[INFO]", pad_label(label), nbad, " /", n
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
    real(dp), intent(in) :: x
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
      write(output_unit,'(a,1x,a)') "[INFO]", "Logging to file: " // trim(LOGFILE)
      write(LOGU,'(a,1x,a)') "[INFO]", "Logging started"
    end if
  end subroutine log_open


  subroutine log_close()
    use, intrinsic :: iso_fortran_env, only : error_unit
    integer :: ios

    if (.not. LOG_TO_FILE) return

    write(LOGU,'(a)') " "
    write(LOGU,'(a,1x,a)') "[INFO]", "Logging ended"
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

