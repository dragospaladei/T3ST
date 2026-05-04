!========================================================================================================
! FAST2 / T3ST :: Turbulent Transport in Tokamaks from Statistical Trajectories
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
!  --/--/2026   | D. I. Palade        | Continuous improvements
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
   CHARACTER(LEN=10000)          :: CWD
   CHARACTER(LEN=:), ALLOCATABLE :: address
   CHARACTER(LEN=:), ALLOCATABLE :: folder
   CHARACTER(LEN=35)             :: runs, nsim

   INTEGER :: lenCWD
   INTEGER :: run
   INTEGER :: simi, simf

   !-----------------------------------------------------------------------------------------------------
   ! Timing
   !-----------------------------------------------------------------------------------------------------
   INTEGER           :: count0, count1, count2, count3
   INTEGER           :: count_rate, count_max
   CHARACTER(LEN=18) :: date_time

   !-----------------------------------------------------------------------------------------------------
   ! Random seeds
   !-----------------------------------------------------------------------------------------------------
   INTEGER, ALLOCATABLE :: seed(:)
   INTEGER              :: ne, clock, ja

   !-----------------------------------------------------------------------------------------------------
   ! Trajectory variables
   !-----------------------------------------------------------------------------------------------------
   REAL(KIND=rp), ALLOCATABLE :: X(:), Y(:), Z(:)
   REAL(KIND=rp), ALLOCATABLE :: Vp(:), mu(:), Ham(:), Pb(:), Pcc(:)
   REAL(KIND=rp), ALLOCATABLE :: ck1(:), ck2(:), ck3(:)
   REAL(KIND=rp), ALLOCATABLE :: Weight(:), FM0(:)
   REAL(KIND=rp), ALLOCATABLE :: q1al(:), q2al(:), q3al(:)

   REAL(KIND=rp), ALLOCATABLE :: q1tot(:, :), q2tot(:, :), q3tot(:, :)

   REAL(KIND=rp), ALLOCATABLE :: Xtraj(:, :), Ytraj(:, :), Ztraj(:, :)
   REAL(KIND=rp), ALLOCATABLE :: Vptraj(:, :), mutraj(:, :), Wetraj(:, :)
   REAL(KIND=rp), ALLOCATABLE :: Htraj(:, :), Pctraj(:, :)
   REAL(KIND=rp), ALLOCATABLE :: q1traj(:, :), q2traj(:, :), q3traj(:, :)
   REAL(KIND=rp), ALLOCATABLE :: ck1traj(:, :), ck2traj(:, :), ck3traj(:, :)

   REAL(KIND=rp)                  :: normB, Vstar1, Vstar2, Vstar3
   REAL(KIND=rp), DIMENSION(3, 3) :: gees

   !-----------------------------------------------------------------------------------------------------
   ! Transport quantities
   !-----------------------------------------------------------------------------------------------------
   REAL(KIND=rp), ALLOCATABLE :: Obs_full(:, :), Obs_delta(:, :)
   REAL(KIND=rp), ALLOCATABLE :: Vcorff(:, :), VcorTT(:, :), VcorTN(:, :), VcorNT(:, :)

   !-----------------------------------------------------------------------------------------------------
   ! Numerical and auxiliary variables
   !-----------------------------------------------------------------------------------------------------
   INTEGER :: i, j, k, m, re, lene, k2, i8
   INTEGER :: filled, ias, eror_flag
   INTEGER :: type_choice

   INTEGER, PARAMETER :: barWidth = 25

   REAL :: progress = 25.0

   REAL(KIND=rp), ALLOCATABLE :: t(:)
   REAL(KIND=rp), ALLOCATABLE :: Lagr_ref(:)
   REAL(KIND=rp), ALLOCATABLE :: Lagr_corr(:)

   REAL(KIND=rp) :: dt1, dt2, dt, dt_half
!   REAL(KIND=rp) :: quan00, quan11, quan12, quan13, quan14, quanVtxCorr
!   REAL(KIND=rp) :: quan21, quan22, quan23, quan24
   INTEGER, PARAMETER :: nq = 21
   REAL(rp) :: quan(nq)
   REAL(KIND=rp) :: gnorm_local

   REAL(rp), POINTER, CONTIGUOUS :: Qx(:), Qy(:), Qz(:), Qw(:), Qph(:), QL(:)

   CHARACTER(LEN=80) :: bar

   LOGICAL :: has_nan, ok

!========================================================================================================
! ADDRESS DECLARATION FOR EXPORTING DATA
!========================================================================================================

   CALL sim_values(type_choice)

   CALL GetCWD(CWD)

   address = TRIM(ADJUSTL(CWD))
   lenCWD  = LEN(address)

   ! Must have 'source' and 'data' folders at same level.
   IF (type_choice == 1) THEN
      address = address(1:lenCWD - 7)//'/data/Sim_'//noofsimstr//'/'
   ELSE
      address = address(1:lenCWD - 7)//'/data/DB_'//noofsimstr//'/'
   END IF

   OPEN(333, FILE=address//'/logs.dat')

!========================================================================================================
! RANDOM SEED
!========================================================================================================

   CALL RANDOM_SEED(SIZE=ne)

   ALLOCATE(seed(ne))

   CALL SYSTEM_CLOCK(COUNT=clock)

   seed = clock + 37 * (/(i - 1, i = 1, ne)/)

   CALL RANDOM_SEED(PUT=seed)

   DEALLOCATE(seed)

!========================================================================================================
! START OF THE RUN LOOP
!========================================================================================================

   simi = 1
   simf = rows1

   DO run = simi, simf

      CALL DATE_AND_TIME(date=date_time(1:8), time=date_time(9:18))
      CALL SYSTEM_CLOCK(count0, count_rate, count_max)
      CALL parameters(run)

!========================================================================================================
! FOLDER DECLARATION
!========================================================================================================

WRITE(runs, '(I5.5)') run
folder = TRIM(address)//'Run_'//TRIM(runs)//'/'


      CALL init_output_files()
      CALL open_output_files(folder)

      WRITE(runs, *) run
      WRITE(nsim, *) (simf - simi + 1)

      lene = LEN(folder)

      CALL printing(runs, nsim, date_time, folder, lene, simf, simi, 333)

!========================================================================================================
! VARIABLE ALLOCATION
!========================================================================================================

      ALLOCATE( X(Np), Xtraj(Nt + 1, ntraj), Obs_full(10, Nt + 1), Obs_delta(10, Nt + 1), &
                Vcorff(Nt + 1, Nt + 1), t(Nt + 1) )

      ALLOCATE( Y, Z, Vp, mu, Ham, Pb, Weight, FM0, q1al, q2al, q3al, Pcc, &
                ck1, ck2, ck3, MOLD=X )

      ALLOCATE( Ytraj, Ztraj, Vptraj, mutraj, Wetraj, Htraj, q1traj, q2traj, &
                q3traj, Pctraj, ck1traj, ck2traj, ck3traj, MOLD=Xtraj )

      ALLOCATE( VcorTT, VcorTN, VcorNT, MOLD=Vcorff )

      ALLOCATE( Lagr_ref(Np), Lagr_corr(Nt + 1) )

!========================================================================================================
! EXPORT PARAMETERS OF THIS RUN
!========================================================================================================

      WRITE(funit(f_parameters), *) '# '//TRIM(sim_comment)

      WRITE(funit(f_parameters), *) [ CHARACTER(LEN=15) :: &
         "t0", "tc", "tt", "tmax", "Nt", "Nreal", "Nc", "Nloop", "Np", "ntraj", "Nci", "Nce", &
         "Ti", "Te", "Zeff", "Aeff", "ndens", "vth", "rhoi", "wi", "Ln", "Li", "Le", &
         "magnetic_model", "B0", "R0", "a0", "s1", "s2", "s3", "q00", "r00", "psi0", &
         "amp", "elong", "device", "shot", "shotslice", "NgridR", "NgridZ", "Omgt0", "Omgtprim", &
         "turb_model", "x_corr", "y_corr", "t_corr", "Phi", "turbprof", "Ai", "Ae", &
         "lambdax", "lambday", "lambdaz", "lbalonz", "tauc", "k0i", "k0e", "gamma_ZF", "gamma_E", &
         "X0", "Y0", "Z0", "Ts", "Es", "pitch", "As", "Zs", "taucc", &
         "position_type", "pitch_type", "energy_type", "Lns", "Lts", &
          "USE_larmor", "USE_coll", "USE_turb", "USE_magnturb", "USE_freq", "USE_polar", &
         "USE_PC", "USE_real", "USE_corr", "USE_balloon", "USE_tilt", "USE_testing", &
         "C1", "C2", "C3" ]

      WRITE(funit(f_parameters), *) [ &
         t0, tc, tt, tmax, REAL(Nt, rp), REAL(Nreal, rp), REAL(Nc, rp), REAL(Nloop, rp), &
         REAL(Np, rp), REAL(ntraj, rp), REAL(Nci, rp), REAL(Nce, rp), Ti, Te, Zeff, Aeff, &
         ndens, vth, rhoi, wi, Ln, Li, Le, REAL(magnetic_model, rp), B0, R0, a0, s1, s2, &
         s3, q00, r00, psi0, amp, elong, REAL(device, rp), REAL(shot, rp), &
         REAL(shotslice, rp), REAL(NgridR, rp), REAL(NgridZ, rp), Omgt0, Omgtprim, &
         REAL(turb_model, rp), REAL(x_corr, rp), REAL(y_corr, rp), REAL(t_corr, rp), &
         Phi, turbprof, Ai, Ae, lambdax, lambday, lambdaz, lbalonz, tauc, k0i, k0e, &
         gamma_ZF, gamma_E, X0, Y0, Z0, Ts, Es, pitch, As, Zs, taucc, &
         REAL(position_type, rp), REAL(pitch_type, rp), REAL(energy_type, rp), &
         REAL(Lns, rp), REAL(Lts, rp),  REAL(USE_larmor, rp), &
         REAL(USE_coll, rp), REAL(USE_turb, rp), REAL(USE_magnturb, rp), REAL(USE_freq, rp), &
         REAL(USE_polar, rp), REAL(USE_PC, rp), REAL(USE_real, rp), REAL(USE_corr, rp), &
         REAL(USE_balloon, rp), REAL(USE_tilt, rp), REAL(USE_testing, rp), C1, C2, C3 ]

!========================================================================================================
! TIME GRID
!========================================================================================================

      dt       = (tmax - t0) / REAL(Nt)
      dt_half  = dt / 2.0_rp
      t        = [(i8 * dt + t0, i8 = 0, Nt)]
      Lagr_ref = 0.0_rp
      Lagr_corr = 0.0_rp

!========================================================================================================
! REALIZATIONS LOOP
!========================================================================================================

      DO re = 1, Nreal
         DO k2 = 1, Nloop

            CALL dispersion(normB, Vstar1, Vstar2, Vstar3, gees)
            CALL wavenum_unified(normB, Vstar1, Vstar2, Vstar3, gees, USE_real, gnorm_local)
            CALL initial_conditions(X, Y, Z, mu, Vp, Pb, Weight, FM0)
            CALL Larmor(USE_larmor, mu)

            IF (USE_testing == ON) THEN
               CALL testing_T3ST(X, Y, Z, Vp, mu, Weight, Pb, Vstar1, Vstar2, Vstar3, &
                                 eror_flag, gnorm_local)
               STOP
            END IF

            CALL write_initial_state(X, Y, Z, Vp, mu, Weight, Pb)

!========================================================================================================
! TIME PROPAGATION
!========================================================================================================

            !$omp parallel default(shared) private(k, i)

            DO k = 1, Nt + 1

               !$omp single
                quan = 0.0_rp
               !$omp end single

               !$omp do schedule(static) private(Qx, Qy, Qz, Qw, Qph, QL) reduction(+:quan)

               DO i = 1, Np
                  BLOCK

                     REAL(rp) :: xi, yi, zi, mui, weighti, vpi, pbi
                     REAL(rp) :: mask, maskR, maskZ
                     REAL(rp) :: vx, vy, vz, vm, vw, vq, ap
                     REAL(rp) :: Wx, Wy, Wz, Wm, Wp, Ww
                     REAL(rp) :: q1, q2, q3, Vrt, Vrr, temp, pbW, FMi, FM0i
                     REAL(rp) :: Hi, B, Vtx, Vty, sm1, sm2, Pc
                     REAL(rp) :: check_1, check_2, check_3

                     ! Select per-particle or single-spectrum arrays.
                     IF (USE_real == OFF) THEN
                        Qx  => kx(:, i)
                        Qy  => ky(:, i)
                        Qz  => kz(:, i)
                        Qw  => w(:, i)
                        Qph => ph(:, i)
                        QL  => L(:, i)
                     ELSE IF (USE_real == ON) THEN
                        Qx  => kxs
                        Qy  => kys
                        Qz  => kzs
                        Qw  => ws
                        Qph => phs
                        QL  => Ls(:, i)
                     END IF

                     ! Load to scalars.
                     xi      = X(i)
                     yi      = Y(i)
                     zi      = Z(i)
                     mui     = mu(i)
                     vpi     = Vp(i)
                     pbi     = Pb(i)
                     weighti = Weight(i)
                     FM0i    = FM0(i)

                     ! Half-collisional Euler-Maruyama step.
                     ! This keeps the Stratonovich interpretation in play.
                     IF ((USE_coll == ON) .AND. (t(k) > tc)) THEN

                        sm1 = rng_uniform()
                        sm2 = rng_uniform()

                        CALL Drift2(dt, xi, yi, zi, vpi, mui, weighti, q1, q2, q3, t(k), &
                                    vx, vy, vz, ap, vm, vw, vq, Hi, Pc, B, Vtx, Vty, &
                                    check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph, QL, &
                                    gnorm_local)

                        CALL collisions_1_MP(sm1, sm2, dt_half, xi, yi, zi, vpi, mui, &
                                             weighti, B, vx, vy, vz, vm, vw, ap)

                        xi      = xi      + dt_half * vx
                        yi      = yi      + dt_half * vy
                        zi      = zi      + dt_half * vz
                        vpi     = vpi     + dt_half * ap
                        mui     = mui     + dt_half * vm
                        weighti = weighti + dt_half * vw

                     END IF

                     ! RK4 stage 1.
                     CALL Drift2(dt, xi, yi, zi, vpi, mui, weighti, q1, q2, q3, t(k), &
                                 vx, vy, vz, ap, vm, vw, vq, Hi, Pc, B, Vtx, Vty, &
                                 check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph, QL, &
                                 gnorm_local)

                     IF (k == 1) Lagr_ref(i) = Vtx

                     Wx = vx / 6.0_rp
                     Wy = vy / 6.0_rp
                     Wz = vz / 6.0_rp
                     Wm = vm / 6.0_rp
                     Wp = ap / 6.0_rp
                     Ww = vw / 6.0_rp

                     ! Mask for particles outside the domain.
                     mask = 1.0_rp

                     IF (ANY(magnetic_model == [1, 2])) THEN
                        maskR = (xi - minR) * (maxR - xi)
                        maskZ = (yi - minZ) * (maxZ - yi)
                        mask  = MERGE(1.0_rp, 0.0_rp, (maskR >= 0.0) .AND. (maskZ >= 0.0))
                     END IF

                     ! Observables.
                     IF (.NOT. ISNAN(q1)) THEN
                        FMi = EXP(-Hi / Ts) / Ts**(1.5_rp)
                        pbW = weighti * FMi / FM0i * Xi !  Xi is here for the same reason as in the initial Pb: the initial sampling (uniform theta) is not consistent with the jacobian of the toroidal coordinates and flux surface homogeneity
                        
                        quan(1:10)   = quan(1:10)   + pbi/Np*mask * [vq, q1, vpi, Hi, Vtx, vq**2, q1**2, vpi**2, Hi**2, Vtx**2]
                        quan(11:20)  = quan(11:20)  + pbW/Np*mask * [vq, q1, vpi, Hi, Vtx, vq**2, q1**2, vpi**2, Hi**2, Vtx**2]
                        quan(21)     = quan(21)     + pbi/Np*mask * (Lagr_ref(i) * Vtx)

                     END IF

                     ! RK4 stage 2.
                     CALL Drift2(dt, xi + vx * dt / 2.0_rp, &
                                     yi + vy * dt / 2.0_rp, &
                                     zi + vz * dt / 2.0_rp, &
                                     vpi + ap * dt / 2.0_rp, &
                                     mui + vm * dt / 2.0_rp, &
                                     weighti + vw * dt / 2.0_rp, &
                                     q1, q2, q3, t(k) + dt / 2.0_rp, &
                                     vx, vy, vz, ap, vm, vw, vq, Hi, Pc, B, Vtx, Vty, &
                                     check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph, QL, &
                                     gnorm_local)

                     Wx = Wx + vx / 3.0_rp
                     Wy = Wy + vy / 3.0_rp
                     Wz = Wz + vz / 3.0_rp
                     Wm = Wm + vm / 3.0_rp
                     Wp = Wp + ap / 3.0_rp
                     Ww = Ww + vw / 3.0_rp

                     ! RK4 stage 3.
                     CALL Drift2(dt, xi + vx * dt / 2.0_rp, &
                                     yi + vy * dt / 2.0_rp, &
                                     zi + vz * dt / 2.0_rp, &
                                     vpi + ap * dt / 2.0_rp, &
                                     mui + vm * dt / 2.0_rp, &
                                     weighti + vw * dt / 2.0_rp, &
                                     q1, q2, q3, t(k) + dt / 2.0_rp, &
                                     vx, vy, vz, ap, vm, vw, vq, Hi, Pc, B, Vtx, Vty, &
                                     check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph, QL, &
                                     gnorm_local)

                     Wx = Wx + vx / 3.0_rp
                     Wy = Wy + vy / 3.0_rp
                     Wz = Wz + vz / 3.0_rp
                     Wm = Wm + vm / 3.0_rp
                     Wp = Wp + ap / 3.0_rp
                     Ww = Ww + vw / 3.0_rp

                     ! RK4 stage 4.
                     CALL Drift2(dt, xi + vx * dt, yi + vy * dt, zi + vz * dt, &
                                 vpi + ap * dt, mui + vm * dt, weighti + vw * dt, &
                                 q1, q2, q3, t(k) + dt, &
                                 vx, vy, vz, ap, vm, vw, vq, Hi, Pc, B, Vtx, Vty, &
                                 check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph, QL, &
                                 gnorm_local)

                     Wx = Wx + vx / 6.0_rp
                     Wy = Wy + vy / 6.0_rp
                     Wz = Wz + vz / 6.0_rp
                     Wm = Wm + vm / 6.0_rp
                     Wp = Wp + ap / 6.0_rp
                     Ww = Ww + vw / 6.0_rp

                     ! Commit back.
                     X(i)      = xi      + dt * Wx
                     Y(i)      = yi      + dt * Wy
                     Z(i)      = zi      + dt * Wz
                     Vp(i)     = vpi     + dt * Wp
                     mu(i)     = mui     + dt * Wm
                     Weight(i) = weighti + dt * Ww

                     Ham(i)  = Hi
                     Pcc(i)  = Pc
                     q1al(i) = q1
                     q2al(i) = q2
                     q3al(i) = q3
                     ck1(i)  = check_1
                     ck2(i)  = check_2
                     ck3(i)  = check_3

                     ! Half-collisional Euler-Maruyama step.
                     IF ((USE_coll == ON) .AND. (t(k) > tc)) THEN

                        xi      = X(i)
                        yi      = Y(i)
                        zi      = Z(i)
                        mui     = mu(i)
                        vpi     = Vp(i)
                        pbi     = Pb(i)
                        weighti = Weight(i)

                        sm1 = rng_uniform()
                        sm2 = rng_uniform()

                        CALL Drift2(dt, xi, yi, zi, vpi, mui, weighti, q1, q2, q3, t(k), &
                                    vx, vy, vz, ap, vm, vw, vq, Hi, Pc, B, Vtx, Vty, &
                                    check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph, QL, &
                                    gnorm_local)

                        CALL collisions_1_MP(sm1, sm2, dt_half, xi, yi, zi, vpi, mui, &
                                             weighti, B, vx, vy, vz, vm, vw, ap)

                        X(i)      = X(i)      + dt_half * vx
                        Y(i)      = Y(i)      + dt_half * vy
                        Z(i)      = Z(i)      + dt_half * vz
                        Vp(i)     = Vp(i)     + dt_half * ap
                        mu(i)     = mu(i)     + dt_half * vm
                        Weight(i) = Weight(i) + dt_half * vw

                     END IF

                  END BLOCK
               END DO

               !$omp end do

               !$omp single

               Obs_full(:, k)  = quan(1:10) 
               Obs_delta(:, k) = quan(11:20)
               Lagr_corr(k) = quan(21) / Np

               Xtraj(k, :)  = X(1:ntraj)
               Ytraj(k, :)  = Y(1:ntraj)
               Ztraj(k, :)  = Z(1:ntraj)
               Vptraj(k, :) = Vp(1:ntraj)
               mutraj(k, :) = mu(1:ntraj)
               Wetraj(k, :) = Weight(1:ntraj)
               Htraj(k, :)  = Ham(1:ntraj)
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

!========================================================================================================
! EXPORT DATA
!========================================================================================================
            CALL write_trajectories(Xtraj, Ytraj, Ztraj, Vptraj, mutraj, Wetraj, Htraj, &
                                    q1traj, q2traj, q3traj, Pctraj, ck1traj, ck2traj, ck3traj)

            CALL write_transport(Obs_full, Obs_delta)
            CALL write_final_state(X, Y, Z, Vp, mu, Weight, Pb, Ham, q1al, q2al)
            CALL write_lagrangian_and_pb(Lagr_corr, Pb)

            IF (USE_corr == 1) THEN
               CALL write_correlations(Vcorff, VcorTT, VcorTN, VcorNT)
            END IF

            IF (ANY(magnetic_model == [1, 2])) THEN

               CALL write_mask(REAL( &
                  COUNT( ((X - minR) * (maxR - X) >= 0.0_rp) .AND. &
                         ((Y - minZ) * (maxZ - Y) >= 0.0_rp) .AND. &
                         .NOT. IEEE_IS_NAN(X) .AND. .NOT. IEEE_IS_NAN(Y) ), rp))

            ELSE IF (ANY(magnetic_model == [3, 4])) THEN

               CALL write_mask(REAL( &
                  COUNT( ((X - (1.0_rp - a0)) * ((1.0_rp + a0) - X) >= 0.0_rp) .AND. &
                         ((Y - (-a0)) * (a0 - Y) >= 0.0_rp) .AND. &
                         .NOT. IEEE_IS_NAN(X) .AND. .NOT. IEEE_IS_NAN(Y) ), rp))

            END IF
!========================================================================================================
! PROGRESS BAR
!========================================================================================================

            progress = REAL(k2 * re) / REAL(Nloop * Nreal)
            filled   = INT(progress * barWidth)

            bar = REPEAT(' ', LEN(bar))

            DO ias = 1, barWidth
               IF (ias <= filled) THEN
                  bar(ias:ias) = '#'
               ELSE
                  bar(ias:ias) = '-'
               END IF
            END DO

            WRITE(*, '(1X,I6,"/",I6," realizations & ",I6,"/",I6," loops [",A,"] ",F6.1,"%",A)', &
                  ADVANCE='yes') &
                  re, Nreal, k2, Nloop, bar(1:barWidth), progress * 100.0, ACHAR(13)

            CALL flush(6)

         END DO

         PRINT *, ''

      END DO

!========================================================================================================
! CLOSE FILES / TIMING / DEALLOCATIONS
!========================================================================================================

      CALL close_output_files()

      CALL SYSTEM_CLOCK(count3, count_rate, count_max)

      WRITE(*, *) ' '
      WRITE(*, *) 'Total computation time  = ', REAL(count3 - count0) / REAL(count_rate)
      WRITE(*, *) '------------------------------------------------------------------------------'

      WRITE(333, *) 'Total computation time  = ', REAL(count3 - count0) / REAL(count_rate)
      WRITE(333, *) '------------------------------------------------------------------------------'

      DEALLOCATE( X, Y, Z, Vp, mu, Ham, Pb, Weight, FM0, Obs_delta, Obs_full, &
                  Xtraj, Ytraj, Ztraj, Vptraj, mutraj, Wetraj, Htraj, &
                  q1traj, q2traj, q3traj, Vcorff, VcorTT, VcorTN, VcorNT, &
                  t, q1al, q2al, q3al, Pcc, ck1, ck2, ck3, Pctraj, &
                  ck1traj, ck2traj, ck3traj, Lagr_ref, Lagr_corr )

      IF ((magnetic_model == 1) .OR. (magnetic_model == 2)) THEN
         DEALLOCATE(Efit_data)
      END IF

   END DO

   CLOSE(333)

END PROGRAM T3ST
