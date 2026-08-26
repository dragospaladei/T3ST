!========================================================================================================
! T3ST :: Turbulent Transport in Tokamaks from Statistical Trajectories
!
! Purpose :: test-particle (DNS) simulations of ion trajectories in turbulent tokamak-like environments
!         :: Tokamak-like MHD equilibrium (realistic or modeled)
!         :: Turbulence via ensemble of random fields with given statistical properties (ITG / TEM)
!         :: Linear analytical dispersion relation; no magnetic perturbations
!         :: For each random field, compute a trajectory from given initial conditions
!         :: Equations of motion consistent with GK gyro-centers
!         :: Projected onto cylindrical coordinates (X,Y,Z)
!         :: Turbulence in field-aligned (q1,q2,q3) ~ (psi, zeta-nu, nu)
!
! Note    :: In the turbulence description, (x,y,z) are field-aligned coordinates.
!         :: Elsewhere, (X,Y,Z) = (R,Z,varphi).
!
! Record of revisions:
!    Date       | Programmer          | Description
!  -----------  | ------------------- | -----------------------------------------
!  --/--/2021   | D. I. Palade        | Initial version
!  --/--/2023   | L. M. Pomarjanschi  | Corrections and development
!  01/01/2024   | D. I. Palade        | Remake
!  01/02/2025   | D. I. Palade        | Further smaller corrections
!  01/09/2025   | D. I. Palade        | Parallelization over Np; performance boost
!  --/--/2026   | D. I. Palade        | Continuous improvements (new day, new code, old me)
!========================================================================================================

PROGRAM T3ST

   USE sims
   USE output_files
   USE constants
   USE random_numbers
   USE omp_lib
   USE drift_kernels
   USE testing_T3ST_mod
   USE rng_mod
   USE, INTRINSIC :: IEEE_ARITHMETIC

   IMPLICIT NONE

!========================================================================================================
! VARIABLE DECLARATION
!========================================================================================================

   !-----------------------------------------------------------------------------------------------------
   ! Folders and simulations
   !-----------------------------------------------------------------------------------------------------
   CHARACTER(LEN=:), ALLOCATABLE :: address
   CHARACTER(LEN=:), ALLOCATABLE :: folder
   CHARACTER(LEN=35) :: runs, nsim

   INTEGER :: run
   INTEGER :: simi, simf

   !-----------------------------------------------------------------------------------------------------
   ! Timing
   !-----------------------------------------------------------------------------------------------------
   INTEGER :: count0, count3
   INTEGER :: count_rate
   CHARACTER(LEN=18) :: date_time

   !-----------------------------------------------------------------------------------------------------
   ! Random seeds
   !-----------------------------------------------------------------------------------------------------
   INTEGER, ALLOCATABLE :: seed(:)
   INTEGER :: ne, clock

   !-----------------------------------------------------------------------------------------------------
   ! Trajectory variables
   !-----------------------------------------------------------------------------------------------------
   REAL(KIND=wp), ALLOCATABLE :: X(:), Y(:), Z(:)
   REAL(KIND=wp), ALLOCATABLE :: Vp(:), mu(:), Einit(:), Ham(:), Pcc(:)
   REAL(KIND=wp), ALLOCATABLE :: ck1(:), ck2(:), ck3(:)
   REAL(KIND=wp), ALLOCATABLE :: F0(:), G0(:), FM(:), betaF(:), betaD(:)
   REAL(KIND=wp), ALLOCATABLE :: alpha_1(:), alpha_2(:), alpha_3(:)
   REAL(KIND=wp), ALLOCATABLE :: q1al(:), q2al(:), q3al(:)
   REAL(KIND=wp), ALLOCATABLE :: Vtxal(:), mask2al(:)

   REAL(KIND=wp), ALLOCATABLE :: Xtraj(:, :), Ytraj(:, :), Ztraj(:, :)
   REAL(KIND=wp), ALLOCATABLE :: Vptraj(:, :), mutraj(:, :)
   REAL(KIND=wp), ALLOCATABLE :: betaFtraj(:, :), betaDtraj(:, :), FMtraj(:, :)
   REAL(KIND=wp), ALLOCATABLE :: Htraj(:, :), Pctraj(:, :)
   REAL(KIND=wp), ALLOCATABLE :: q1traj(:, :), q2traj(:, :), q3traj(:, :)
   REAL(KIND=wp), ALLOCATABLE :: ck1traj(:, :), ck2traj(:, :), ck3traj(:, :)

   REAL(KIND=wp) :: normB, Vstar1, Vstar2, Vstar3
   REAL(KIND=wp), DIMENSION(3, 3) :: gees

   !-----------------------------------------------------------------------------------------------------
   ! Transport quantities
   !-----------------------------------------------------------------------------------------------------
   REAL(KIND=wp), ALLOCATABLE :: Obs_GF(:, :), Obs_FF(:, :), Obs_FD(:, :), Obs_DF(:, :), Obs_DD(:, :)
   REAL(KIND=wp), ALLOCATABLE :: Obs_SS(:, :)
   REAL(KIND=wp), ALLOCATABLE :: Vcorff(:, :), VcorTT(:, :), VcorTN(:, :), VcorNT(:, :)

   !-----------------------------------------------------------------------------------------------------
   ! Numerical and auxiliary variables
   !-----------------------------------------------------------------------------------------------------
   INTEGER :: i, k, re, lene, k2, i8
   INTEGER :: filled, ias, eror_flag
   INTEGER :: diag_loops_checked, diag_loops_failed, diag_checks_failed
   INTEGER, PARAMETER :: n_diag_fail_types = 8
   INTEGER :: diag_fail_counts(n_diag_fail_types), diag_fail_counts_loop(n_diag_fail_types)
   INTEGER :: type_choice
   INTEGER :: arg_count, arg_type_choice, arg_sim_number, arg_status

   INTEGER, PARAMETER :: barWidth = 25

   REAL :: progress = 25.0

   REAL(KIND=wp), ALLOCATABLE :: t(:)
   REAL(KIND=wp), ALLOCATABLE :: Vtx_init(:)
   REAL(KIND=wp), ALLOCATABLE :: kin_init(:)
   REAL(KIND=wp), ALLOCATABLE :: Lagr_corr(:)

   REAL(KIND=wp) :: dt, dt_half
   REAL(KIND=wp) :: q1min, q1max
   INTEGER, PARAMETER :: nq = 55
   REAL(wp) :: quan(nq)
   REAL(wp), POINTER, CONTIGUOUS :: Qx(:), Qy(:), Qz(:), Qw(:), Qph(:)

   CHARACTER(LEN=80) :: bar
   CHARACTER(LEN=80) :: arg

!========================================================================================================
! ADDRESS DECLARATION FOR EXPORTING DATA
!========================================================================================================

   arg_count = COMMAND_ARGUMENT_COUNT()

   IF (arg_count == 0) THEN
      CALL load_simulation_matrix(type_choice)
   ELSEIF (arg_count == 2) THEN
      CALL GET_COMMAND_ARGUMENT(1, arg)
      READ (arg, *, IOSTAT=arg_status) arg_type_choice
      IF (arg_status /= 0) THEN
         ERROR STOP 'First argument must be 1 (Sim) or 2 (DB).'
      END IF

      CALL GET_COMMAND_ARGUMENT(2, arg)
      READ (arg, *, IOSTAT=arg_status) arg_sim_number
      IF (arg_status /= 0) THEN
         ERROR STOP 'Second argument must be the Sim/DB file number.'
      END IF

      CALL load_simulation_matrix(type_choice, arg_type_choice, arg_sim_number)
   ELSE
      ERROR STOP 'Usage: T3ST [file_type file_number], where file_type is 1=Sim or 2=DB.'
   END IF

   IF (type_choice == 1) THEN
      address = project_root//'/data/Sim_'//noofsimstr//'/'
   ELSE
      address = project_root//'/data/DB_'//noofsimstr//'/'
   END IF

   OPEN (333, FILE=address//'/logs.dat')

!========================================================================================================
! RANDOM SEED
!========================================================================================================

   CALL RANDOM_SEED(SIZE=ne)

   ALLOCATE (seed(ne))

   CALL SYSTEM_CLOCK(COUNT=clock)

   seed = clock + 37*(/(i - 1, i=1, ne)/)

   CALL RANDOM_SEED(PUT=seed)

   DEALLOCATE (seed)

!========================================================================================================
! START OF THE RUN LOOP
!========================================================================================================

   simi = 1
   simf = rows1

   DO run = simi, simf

      CALL DATE_AND_TIME(date=date_time(1:8), time=date_time(9:18))
      CALL SYSTEM_CLOCK(count0, count_rate)
      CALL load_run_parameters(run)

!========================================================================================================
! FOLDER DECLARATION
!========================================================================================================

      WRITE (runs, '(I5.5)') run
      folder = TRIM(address)//'Run_'//TRIM(runs)//'/'

      CALL init_output_files()
      CALL open_output_files(folder)

      WRITE (runs, *) run
      WRITE (nsim, *) (simf - simi + 1)

      lene = LEN(folder)

      CALL print_run_header(runs, nsim, date_time, folder, lene, simf, simi, 333)

!========================================================================================================
! VARIABLE ALLOCATION
!========================================================================================================

      ALLOCATE (X(Np), Xtraj(Nt + 1, ntraj), Obs_GF(10, Nt + 1), Obs_FF(10, Nt + 1), &
                Obs_FD(10, Nt + 1), Obs_DF(10, Nt + 1), Obs_DD(10, Nt + 1), &
                Obs_SS(4, Nt + 1), Vcorff(Nt + 1, Nt + 1), t(Nt + 1))

      ALLOCATE (Y, Z, Vp, mu, Einit, Ham, F0, G0, FM, betaF, betaD, q1al, q2al, q3al, Vtxal, mask2al, Pcc, &
                ck1, ck2, ck3, alpha_1, alpha_2, alpha_3, MOLD=X)

      ALLOCATE (Ytraj, Ztraj, Vptraj, mutraj, betaFtraj, betaDtraj, FMtraj, Htraj, q1traj, q2traj, &
                q3traj, Pctraj, ck1traj, ck2traj, ck3traj, MOLD=Xtraj)

      ALLOCATE (VcorTT, VcorTN, VcorNT, MOLD=Vcorff)

      ALLOCATE (Vtx_init(Np), kin_init(Np), Lagr_corr(Nt + 1))

!========================================================================================================
! EXPORT PARAMETERS OF THIS RUN
!========================================================================================================

      WRITE (funit(f_parameters), *) '# '//TRIM(sim_comment)

      WRITE (funit(f_parameters), *) [CHARACTER(LEN=15) :: &
                                      "t0", "tc", "tt", "tmax", "Nt", "Nreal", "Nc", "Nloop", "Np", "ntraj", "Nci", "Nce", &
                                      "Ti", "Te", "Zeff", "Aeff", "ndens", "vth", "rhoi", "wi", "Ln", "Li", "Le", &
                                 "magnetic_model", "B0", "R0", "a0", "safety", "shearmag", "s1", "s2", "s3", "q00", "r00", "psi0", &
                                      "elong", "device", "shot", "shotslice", "NgridR", "NgridZ", "Omgt0", "Omgtprim", &
                                      "turb_model", "x_corr", "y_corr", "t_corr", "Phi", "beta", "turbprof", "Ai", "Ae", &
                                      "lambdax", "lambday", "lambdaz", "lbalonz", "tauc", "k0i", "k0e", "gamma_ZF", "gamma_E", &
                                      "X0", "Y0", "Z0", "Ts", "Es", "pitch", "As", "Zs", "taucc", &
                                      "position_type", "annulus_width", "pitch_type", "energy_type", "Lns", "Lts", &
                                      "USE_larmor", "USE_coll", "USE_turb", "USE_magnturb", "USE_freq", "USE_polar", &
                                      "USE_PC", "USE_real", "USE_corr", "USE_balloon", "USE_tilt", "USE_testing", &
                                      "C1", "C2", "C3"]

      WRITE (funit(f_parameters), *) [ &
         t0, tc, tt, tmax, REAL(Nt, wp), REAL(Nreal, wp), REAL(Nc, wp), REAL(Nloop, wp), &
         REAL(Np, wp), REAL(ntraj, wp), REAL(Nci, wp), REAL(Nce, wp), Ti, Te, Zeff, Aeff, &
         ndens, vth, rhoi, wi, Ln, Li, Le, REAL(magnetic_model, wp), B0, R0, a0, safety, shear_mag, s1, s2, &
         s3, q00, r00, psi0, elong, REAL(device, wp), REAL(shot, wp), &
         REAL(shotslice, wp), REAL(NgridR, wp), REAL(NgridZ, wp), Omgt0, Omgtprim, &
         REAL(turb_model, wp), REAL(x_corr, wp), REAL(y_corr, wp), REAL(t_corr, wp), &
         Phi, Beta, turbprof, Ai, Ae, lambdax, lambday, lambdaz, lbalonz, tauc, k0i, k0e, &
         gamma_ZF, gamma_E, X0, Y0, Z0, Ts, Es, pitch, As, Zs, taucc, &
         REAL(position_type, wp), annulus_width, REAL(pitch_type, wp), REAL(energy_type, wp), &
         REAL(Lns, wp), REAL(Lts, wp), REAL(USE_larmor, wp), &
         REAL(USE_coll, wp), REAL(USE_turb, wp), REAL(USE_magnturb, wp), REAL(USE_freq, wp), &
         REAL(USE_polar, wp), REAL(USE_PC, wp), REAL(USE_real, wp), REAL(USE_corr, wp), &
         REAL(USE_balloon, wp), REAL(USE_tilt, wp), REAL(USE_testing, wp), C1, C2, C3]

!========================================================================================================
! TIME GRID   & minimal diagnostic counters
!========================================================================================================

      dt = (tmax - t0)/REAL(Nt)
      dt_half = dt/2.0_WP
      t = [(i8*dt + t0, i8=0, Nt)]
      Vtx_init = 0.0_WP
      Lagr_corr = 0.0_WP
      q1min = C1*MAX(0.0_WP, rhot0 - 0.1_WP*annulus_width)
      q1max = C1*MIN(1.0_WP, rhot0 + 0.1_WP*annulus_width)

      diag_loops_checked = 0
      diag_loops_failed = 0
      diag_checks_failed = 0
      diag_fail_counts = 0

!========================================================================================================
! REALIZATIONS LOOP
!========================================================================================================

      DO re = 1, Nreal
         DO k2 = 1, Nloop

            CALL dispersion(normB, Vstar1, Vstar2, Vstar3, gees)
            CALL wavenum_unified(normB, Vstar1, Vstar2, Vstar3, gees, USE_real)
            CALL initialize_particles(X, Y, Z, mu, Vp, Einit, F0, G0, FM, betaF, betaD, alpha_1, alpha_2, alpha_3)
            CALL Larmor(USE_larmor, mu)

            IF (USE_testing == ON) THEN
               CALL testing_T3ST(X, Y, Z, Vp, mu, F0, G0, FM, betaF, betaD, Vstar1, Vstar2, Vstar3, &
                                 eror_flag)
               STOP
            END IF

            CALL write_initial_state(X, Y, Z, Vp, mu, Einit, betaF, betaD, FM, F0, G0)

            !========================================================================================================
! TIME PROPAGATION
!========================================================================================================

            !$omp parallel default(shared) private(k, i)

            DO k = 1, Nt + 1

               !$omp single
               quan = 0.0_WP
               !$omp end single

               !$omp do schedule(static) private(Qx, Qy, Qz, Qw, Qph) reduction(+:quan)

               DO i = 1, Np
                  BLOCK

                     REAL(wp) :: xi, yi, zi, mui, betaFi, betaDi, vpi
                     REAL(wp) :: alpha_1i, alpha_2i, alpha_3i
                     REAL(wp) :: mask, mask2, maskR, maskZ
                     REAL(wp) :: vx, vy, vz, vm, ap, vbF, vbD
                     REAL(wp) :: Wx, Wy, Wz, Wm, Wvp, WbF, WbD
                     REAL(wp) :: Walpha_1, Walpha_2, Walpha_3
                     REAL(wp) :: q1, q2, q3, FMi, F0i, G0i
                     REAL(wp) :: pbGF, pbFF, pbFD, pbDF, pbDD, pbSS
                     REAL(wp) :: Hi, B, Vtx, VFx, sm1, sm2, Pc, kin
                     REAL(wp) :: check_1, check_2, check_3, vs_1, vs_2, vs_3

                     ! Select per-particle or single-spectrum arrays.
                     IF (USE_real == OFF) THEN
                        Qx => kx(:, i)
                        Qy => ky(:, i)
                        Qz => kz(:, i)
                        Qw => w(:, i)
                        Qph => ph(:, i)
                     ELSE IF (USE_real == ON) THEN
                        Qx => kxs
                        Qy => kys
                        Qz => kzs
                        Qw => ws
                        Qph => phs
                     END IF

                     ! Load to scalars.
                     xi = X(i)
                     yi = Y(i)
                     zi = Z(i)
                     mui = mu(i)
                     vpi = Vp(i)
                     F0i = F0(i)
                     G0i = G0(i)
                     betaFi = betaF(i)
                     betaDi = betaD(i)
                     alpha_1i = alpha_1(i)
                     alpha_2i = alpha_2(i)
                     alpha_3i = alpha_3(i)

                     ! Half-collisional Euler-Maruyama step.
                     ! This keeps the Stratonovich interpretation in play.
                     IF ((USE_coll == ON) .AND. (t(k) > tc)) THEN

                        sm1 = rng_uniform()
                        sm2 = rng_uniform()

                        CALL gyrocenter_drifts(xi, yi, zi, vpi, mui, q1, q2, q3, t(k), &
                                               vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx, &
                                               check_1, check_2, check_3, vs_1, vs_2, vs_3, Qx, Qy, Qz, Qw, Qph, Sol_data)

                        CALL collision_kicks(sm1, sm2, dt_half, xi, yi, zi, vpi, mui, &
                                             B, vx, vy, vz, vm, ap)

                        xi = xi + dt_half*vx
                        yi = yi + dt_half*vy
                        zi = zi + dt_half*vz
                        vpi = vpi + dt_half*ap
                        mui = mui + dt_half*vm
                        betaFi = betaFi + dt_half*vbF
                        betaDi = betaDi + dt_half*vbD
                        alpha_1i = alpha_1i + abs(dt_half)*vs_1
                        alpha_2i = alpha_2i + abs(dt_half)*vs_2
                        alpha_3i = alpha_3i + abs(dt_half)*vs_3

                     END IF

                     ! RK4 stage 1.
                     CALL gyrocenter_drifts(xi, yi, zi, vpi, mui, q1, q2, q3, t(k), &
                                            vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx, &
                                            check_1, check_2, check_3, vs_1, vs_2, vs_3, Qx, Qy, Qz, Qw, Qph, Sol_data)

                     IF (k == 1) THEN
                        Vtx_init(i) = Vtx
                        kin_init(i) = kin
                     END IF

                     Wx = vx/6.0_WP
                     Wy = vy/6.0_WP
                     Wz = vz/6.0_WP
                     Wm = vm/6.0_WP
                     Wvp = ap/6.0_WP
                     WbF = vbF/6.0_WP
                     WbD = vbD/6.0_WP
                     Walpha_1 = vs_1/6.0_WP
                     Walpha_2 = vs_2/6.0_WP
                     Walpha_3 = vs_3/6.0_WP

                     ! Mask for particles outside the domain.
                     mask = 1.0_WP

                     IF (magnetic_model == 1) THEN
                        maskR = (xi - minR)*(maxR - xi)
                        maskZ = (yi - minZ)*(maxZ - yi)
                        mask = MERGE(1.0_WP, 0.0_WP, (maskR >= 0.0) .AND. (maskZ >= 0.0))
                     END IF
                     mask2 = MERGE(1.0_WP, 0.0_WP, (q1 > q1min) .AND. (q1 < q1max))

                     ! Observables.
                     IF (.NOT. ISNAN(q1)) THEN
                        pbGF = F0i/G0i/Np
                        pbFF = FMi*exp(betaFi)/G0i/Np
                        pbFD = FMi*exp(betaDi)/G0i/Np
                        pbDF = FMi*(exp(betaFi) - 1.0_WP)/G0i/Np
                        pbDD = FMi*(exp(betaDi) - 1.0_WP)/G0i/Np
                        pbSS = FMi/G0i/Np

                        ! Obs columns: Vtx1, Q1, Vp1, H1, Heat1, Vtx2, Q2, Vp2, H2, Heat2.
                        quan(1:10) = quan(1:10) + pbGF*mask &
                                     *[Vtx, q1, vpi, Hi, Vtx*kin, Vtx**2, q1**2, vpi**2, Hi**2, Vtx*kin]
                        quan(11:20) = quan(11:20) + pbFF*mask &
                                      *[Vtx, q1, vpi, Hi, Vtx*kin, Vtx**2, q1**2, vpi**2, Hi**2, Vtx*kin]
                        quan(21:30) = quan(21:30) + pbFD*mask &
                                      *[Vtx, q1, vpi, Hi, Vtx*kin, Vtx**2, q1**2, vpi**2, Hi**2, Vtx*kin]
                        quan(31:40) = quan(31:40) + pbDF*mask &
                                      *[Vtx, q1, vpi, Hi, Vtx*kin, Vtx**2, q1**2, vpi**2, Hi**2, Vtx*kin]
                        quan(41:50) = quan(41:50) + pbDD*mask &
                                      *[Vtx, q1, vpi, Hi, Vtx*kin, Vtx**2, q1**2, vpi**2, Hi**2, Vtx*kin]
                        quan(51) = quan(51) + pbSS*mask*Vtx_init(i)*Vtx         ! green function
                        ! Signs follow the convection-diffusion convention G/n = V - D grad(log n).
                        quan(52) = quan(52) + pbSS*mask*Vtx_init(i)*alpha_1i
                        quan(53) = quan(53) - pbSS*mask*Vtx_init(i)*alpha_2i
                        quan(54) = quan(54) - pbSS*mask*Vtx_init(i)*alpha_3i
                        quan(55) = quan(55) - pbSS*mask*Vtx_init(i)*alpha_3i*kin_init(i)

                     END IF

!             IF (i == 12) THEN
!               write(*,*) 'rk', exp((q1 - C1*rhot0)*a0/C1/R0*Lns), q1 - C1*rhot0, kin, FMi, F0i, exp((q1 - C1*rhot0)*a0/C1/R0*Lns)*exp(-kin/ Ts) / Ts**1.5_wp
!             write(*,*) 'rk', exp((q1 - C1*rhot0)*a0/C1/R0*Lns), FMi!, F0i, exp((q1 - C1*rhot0)*a0/C1/R0*Lns)*exp(-kin/ Ts) / Ts**1.5_wp
!                  WRITE(*, '("particle | ",I0,"| FM*Exp/F0 | ",G0.6,"| deltaFM/F0 | ",G0.6,"| kin/Ts | ",G0.6)') &
!            i, 100.0*(FMi*exp(betaFi)/F0i-1.0), 100.0*(FMi/F0i-1.0_wp), kin-Hi
!           END IF
                     ! RK4 stage 2.
                     CALL gyrocenter_drifts(xi + vx*dt/2.0_WP, &
                                            yi + vy*dt/2.0_WP, &
                                            zi + vz*dt/2.0_WP, &
                                            vpi + ap*dt/2.0_WP, &
                                            mui + vm*dt/2.0_WP, &
                                            q1, q2, q3, t(k) + dt/2.0_WP, &
                                            vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx, &
                                            check_1, check_2, check_3, vs_1, vs_2, vs_3, Qx, Qy, Qz, Qw, Qph, Sol_data)

                     Wx = Wx + vx/3.0_WP
                     Wy = Wy + vy/3.0_WP
                     Wz = Wz + vz/3.0_WP
                     Wm = Wm + vm/3.0_WP
                     Wvp = Wvp + ap/3.0_WP
                     WbF = WbF + vbF/3.0_WP
                     WbD = WbD + vbD/3.0_WP
                     Walpha_1 = Walpha_1 + vs_1/3.0_WP
                     Walpha_2 = Walpha_2 + vs_2/3.0_WP
                     Walpha_3 = Walpha_3 + vs_3/3.0_WP

                     ! RK4 stage 3.
                     CALL gyrocenter_drifts(xi + vx*dt/2.0_WP, &
                                            yi + vy*dt/2.0_WP, &
                                            zi + vz*dt/2.0_WP, &
                                            vpi + ap*dt/2.0_WP, &
                                            mui + vm*dt/2.0_WP, &
                                            q1, q2, q3, t(k) + dt/2.0_WP, &
                                            vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx, &
                                            check_1, check_2, check_3, vs_1, vs_2, vs_3, Qx, Qy, Qz, Qw, Qph, Sol_data)

                     Wx = Wx + vx/3.0_WP
                     Wy = Wy + vy/3.0_WP
                     Wz = Wz + vz/3.0_WP
                     Wm = Wm + vm/3.0_WP
                     Wvp = Wvp + ap/3.0_WP
                     WbF = WbF + vbF/3.0_WP
                     WbD = WbD + vbD/3.0_WP
                     Walpha_1 = Walpha_1 + vs_1/3.0_WP
                     Walpha_2 = Walpha_2 + vs_2/3.0_WP
                     Walpha_3 = Walpha_3 + vs_3/3.0_WP

                     ! RK4 stage 4.
                     CALL gyrocenter_drifts(xi + vx*dt, yi + vy*dt, zi + vz*dt, &
                                            vpi + ap*dt, mui + vm*dt, &
                                            q1, q2, q3, t(k) + dt, &
                                            vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx, &
                                            check_1, check_2, check_3, vs_1, vs_2, vs_3, Qx, Qy, Qz, Qw, Qph, Sol_data)

                     Wx = Wx + vx/6.0_WP
                     Wy = Wy + vy/6.0_WP
                     Wz = Wz + vz/6.0_WP
                     Wm = Wm + vm/6.0_WP
                     Wvp = Wvp + ap/6.0_WP
                     WbF = WbF + vbF/6.0_WP
                     WbD = WbD + vbD/6.0_WP
                     Walpha_1 = Walpha_1 + vs_1/6.0_WP
                     Walpha_2 = Walpha_2 + vs_2/6.0_WP
                     Walpha_3 = Walpha_3 + vs_3/6.0_WP

                     ! Commit back.
                     X(i) = xi + dt*Wx
                     Y(i) = yi + dt*Wy
                     Z(i) = zi + dt*Wz
                     Vp(i) = vpi + dt*Wvp
                     mu(i) = mui + dt*Wm
                     betaF(i) = betaFi + dt*WbF
                     betaD(i) = betaDi + dt*WbD
                     alpha_1(i) = alpha_1i + abs(dt)*Walpha_1
                     alpha_2(i) = alpha_2i + abs(dt)*Walpha_2
                     alpha_3(i) = alpha_3i + abs(dt)*Walpha_3

                     Ham(i) = Hi
                     FM(i) = FMi
                     Pcc(i) = Pc
                     q1al(i) = q1
                     q2al(i) = q2
                     q3al(i) = q3
                     Vtxal(i) = Vtx
                     mask2al(i) = MERGE(1.0_WP, 0.0_WP, (q1 > q1min) .AND. (q1 < q1max))
                     ck1(i) = check_1
                     ck2(i) = check_2
                     ck3(i) = check_3

                     ! Half-collisional Euler-Maruyama step.
                     IF ((USE_coll == ON) .AND. (t(k) > tc)) THEN

                        xi = X(i)
                        yi = Y(i)
                        zi = Z(i)
                        mui = mu(i)
                        vpi = Vp(i)
                        betaFi = betaF(i)
                        betaDi = betaD(i)
                        alpha_1i = alpha_1(i)
                        alpha_2i = alpha_2(i)
                        alpha_3i = alpha_3(i)

                        sm1 = rng_uniform()
                        sm2 = rng_uniform()

                        CALL gyrocenter_drifts(xi, yi, zi, vpi, mui, q1, q2, q3, t(k), &
                                               vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx, &
                                               check_1, check_2, check_3, vs_1, vs_2, vs_3, Qx, Qy, Qz, Qw, Qph, Sol_data)

                        CALL collision_kicks(sm1, sm2, dt_half, xi, yi, zi, vpi, mui, &
                                             B, vx, vy, vz, vm, ap)

                        X(i) = X(i) + dt_half*vx
                        Y(i) = Y(i) + dt_half*vy
                        Z(i) = Z(i) + dt_half*vz
                        Vp(i) = Vp(i) + dt_half*ap
                        mu(i) = mu(i) + dt_half*vm
                        betaF(i) = betaF(i) + dt_half*vbF
                        betaD(i) = betaD(i) + dt_half*vbD
                        ! These are quadratures, so use |dt| even for backward trajectories.
                        alpha_1(i) = alpha_1i + abs(dt_half)*vs_1
                        alpha_2(i) = alpha_2i + abs(dt_half)*vs_2
                        alpha_3(i) = alpha_3i + abs(dt_half)*vs_3

                     END IF

                  END BLOCK
               END DO

               !$omp end do

               !$omp single

               Obs_GF(:, k) = quan(1:10)    ! green function
               Obs_FF(:, k) = quan(11:20)  ! fullF+fullV
               Obs_FD(:, k) = quan(21:30)
               Obs_DF(:, k) = quan(31:40)  ! fullF+fullV
               Obs_DD(:, k) = quan(41:50)
               Obs_SS(:, k) = quan(52:55)
               Lagr_corr(k) = quan(51)

               Xtraj(k, :) = X(1:ntraj)
               Ytraj(k, :) = Y(1:ntraj)
               Ztraj(k, :) = Z(1:ntraj)
               Vptraj(k, :) = Vp(1:ntraj)
               mutraj(k, :) = mu(1:ntraj)
               betaFtraj(k, :) = betaF(1:ntraj)
               betaDtraj(k, :) = betaD(1:ntraj)
               FMtraj(k, :) = FM(1:ntraj)
               Htraj(k, :) = Ham(1:ntraj)
               Pctraj(k, :) = Pcc(1:ntraj)

               ck1traj(k, :) = ck1(1:ntraj)
               ck2traj(k, :) = ck2(1:ntraj)
               ck3traj(k, :) = ck3(1:ntraj)

               q1traj(k, :) = q1al(1:ntraj)
               q2traj(k, :) = q2al(1:ntraj)
               q3traj(k, :) = q3al(1:ntraj)

               !$omp end single

               !$omp barrier

            END DO

            !$omp end parallel

            CALL diagnose_saved_trajectories(F0, FMtraj, betaFtraj, Htraj, Pctraj, &
                                             Obs_FF, Obs_FD, Obs_DF, Obs_DD, eror_flag, diag_fail_counts_loop)

            diag_loops_checked = diag_loops_checked + 1
            diag_checks_failed = diag_checks_failed + eror_flag
            diag_fail_counts = diag_fail_counts + diag_fail_counts_loop

            IF (eror_flag > 0) THEN
               diag_loops_failed = diag_loops_failed + 1
            END IF

!========================================================================================================
! EXPORT DATA
!========================================================================================================
            CALL write_trajectories(Xtraj, Ytraj, Ztraj, Vptraj, mutraj, betaFtraj, betaDtraj, FMtraj, Htraj, &
                                    q1traj, q2traj, q3traj, Pctraj, ck1traj, ck2traj, ck3traj)

            CALL write_transport(Obs_GF, Obs_FF, Obs_FD, Obs_DF, Obs_DD, Obs_SS)
            CALL write_final_state(X, Y, Z, Vp, mu, betaF, betaD, FM, Ham, q1al, q2al, Vtxal, mask2al, &
                                   F0, G0, Vtx_init, alpha_1, alpha_2, alpha_3)
            CALL write_lagrangian(Lagr_corr)

            IF (USE_corr == 1) THEN
               CALL write_correlations(Vcorff, VcorTT, VcorTN, VcorNT)
            END IF

            IF (magnetic_model == 1) THEN

               CALL write_mask(REAL( &
                               COUNT(((X - minR)*(maxR - X) >= 0.0_WP) .AND. &
                                     ((Y - minZ)*(maxZ - Y) >= 0.0_WP) .AND. &
                                     .NOT. IEEE_IS_NAN(X) .AND. .NOT. IEEE_IS_NAN(Y)), wp))

            ELSE IF (ANY(magnetic_model == [2, 3])) THEN

               CALL write_mask(REAL( &
                               COUNT(((X - (1.0_WP - a0))*((1.0_WP + a0) - X) >= 0.0_WP) .AND. &
                                     ((Y - (-a0))*(a0 - Y) >= 0.0_WP) .AND. &
                                     .NOT. IEEE_IS_NAN(X) .AND. .NOT. IEEE_IS_NAN(Y)), wp))

            END IF
!========================================================================================================
! PROGRESS BAR
!========================================================================================================

            progress = REAL((re - 1)*Nloop + k2)/REAL(Nloop*Nreal)
            filled = INT(progress*barWidth)

            bar = REPEAT(' ', LEN(bar))

            DO ias = 1, barWidth
               IF (ias <= filled) THEN
                  bar(ias:ias) = '#'
               ELSE
                  bar(ias:ias) = '-'
               END IF
            END DO

            WRITE (*, '(1X,I6,"/",I6," realizations & ",I6,"/",I6," loops [",A,"] ",F6.1,"%",A)', &
                   ADVANCE='yes') &
               re, Nreal, k2, Nloop, bar(1:barWidth), progress*100.0, ACHAR(13)

            CALL flush (6)

         END DO

         PRINT *, ''

      END DO
!========================================================================================================
! MINIMAL DIAGNOSIS RESULTS
!========================================================================================================
      IF (diag_loops_failed == 0) THEN
         WRITE (*, *) 'Minimal diagnosis check: PASSED'
      ELSE
         WRITE (*, *) 'Minimal diagnosis check: FAILED'
         WRITE (*, *) '  checked loops: ', diag_loops_checked
         WRITE (*, *) '  failed loops : ', diag_loops_failed
         WRITE (*, *) '  failed checks: ', diag_checks_failed
         WRITE (*, *) '  failure breakdown:'
         WRITE (*, *) '    shape/size checks              : ', diag_fail_counts(1)
         WRITE (*, *) '    F0 <= 0                        : ', diag_fail_counts(2)
         WRITE (*, *) '    FM*exp(betaF)/F0 identity      : ', diag_fail_counts(3)
         WRITE (*, *) '    Htraj(:,1) time constancy      : ', diag_fail_counts(4)
         WRITE (*, *) '    Htraj(:,2) time constancy      : ', diag_fail_counts(5)
         WRITE (*, *) '    Pctraj(:,1) time constancy     : ', diag_fail_counts(6)
         WRITE (*, *) '    Obs_FF-Obs_FD vs Obs_DF-Obs_DD : ', diag_fail_counts(7)
         WRITE (*, *) '    Obs_FD-Obs_DD vs Obs_FF-Obs_DF : ', diag_fail_counts(8)
      END IF

      WRITE (333, *) 'Saved-trajectory diagnostics: ', MERGE('GOOD', 'BAD ', diag_loops_failed == 0)
      WRITE (333, *) '  checked loops: ', diag_loops_checked
      WRITE (333, *) '  failed loops : ', diag_loops_failed
      WRITE (333, *) '  failed checks: ', diag_checks_failed
      WRITE (333, *) '  failure breakdown:'
      WRITE (333, *) '    shape/size checks              : ', diag_fail_counts(1)
      WRITE (333, *) '    F0 <= 0                        : ', diag_fail_counts(2)
      WRITE (333, *) '    FM*exp(betaF)/F0 identity      : ', diag_fail_counts(3)
      WRITE (333, *) '    Htraj(:,1) time constancy      : ', diag_fail_counts(4)
      WRITE (333, *) '    Htraj(:,2) time constancy      : ', diag_fail_counts(5)
      WRITE (333, *) '    Pctraj(:,1) time constancy     : ', diag_fail_counts(6)
      WRITE (333, *) '    Obs_FF-Obs_FD vs Obs_DF-Obs_DD : ', diag_fail_counts(7)
      WRITE (333, *) '    Obs_FD-Obs_DD vs Obs_FF-Obs_DF : ', diag_fail_counts(8)
      WRITE (333, *) '------------------------------------------------------------------------------'

!========================================================================================================
! CLOSE FILES / TIMING / DEALLOCATIONS
!========================================================================================================

      CALL close_output_files()

      CALL SYSTEM_CLOCK(count3, count_rate)

      WRITE (*, *) ' '
      WRITE (*, *) 'Total computation time  = ', REAL(count3 - count0)/REAL(count_rate)
      WRITE (*, *) '------------------------------------------------------------------------------'

      WRITE (333, *) 'Total computation time  = ', REAL(count3 - count0)/REAL(count_rate)
      WRITE (333, *) '------------------------------------------------------------------------------'

      DEALLOCATE (X, Y, Z, Vp, mu, Einit, Ham, F0, G0, FM, betaF, betaD, alpha_1, alpha_2, alpha_3, &
                  Obs_GF, Obs_FF, Obs_FD, Obs_DF, Obs_DD, Obs_SS, &
                  Xtraj, Ytraj, Ztraj, Vptraj, mutraj, betaFtraj, betaDtraj, FMtraj, Htraj, &
                  q1traj, q2traj, q3traj, Vcorff, VcorTT, VcorTN, VcorNT, &
                  t, q1al, q2al, q3al, Vtxal, mask2al, Pcc, ck1, ck2, ck3, Pctraj, &
                  ck1traj, ck2traj, ck3traj, Vtx_init, kin_init, Lagr_corr)

      IF (magnetic_model == 1) THEN
         DEALLOCATE (Efit_data)
      END IF

   END DO

   CLOSE (333)

END PROGRAM T3ST
