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
   CHARACTER(LEN=10000)          :: CWD
   CHARACTER(LEN=:), ALLOCATABLE :: address
   CHARACTER(LEN=:), ALLOCATABLE :: folder
   CHARACTER(LEN=35)             :: runs, nsim

   INTEGER :: run
   INTEGER :: simi, simf

   !-----------------------------------------------------------------------------------------------------
   ! Timing
   !-----------------------------------------------------------------------------------------------------
   INTEGER           :: count0, count3
   INTEGER           :: count_rate, count_max
   CHARACTER(LEN=18) :: date_time

   !-----------------------------------------------------------------------------------------------------
   ! Random seeds
   !-----------------------------------------------------------------------------------------------------
   INTEGER, ALLOCATABLE :: seed(:)
   INTEGER              :: ne, clock

   !-----------------------------------------------------------------------------------------------------
   ! Trajectory variables
   !-----------------------------------------------------------------------------------------------------
   REAL(KIND=rp), ALLOCATABLE :: X(:), Y(:), Z(:)
   REAL(KIND=rp), ALLOCATABLE :: Vp(:), mu(:), Ham(:), Pcc(:)
   REAL(KIND=rp), ALLOCATABLE :: ck1(:), ck2(:), ck3(:)
   REAL(KIND=rp), ALLOCATABLE :: F0(:), G0(:), FM(:), betaF(:), betaD(:)
   REAL(KIND=rp), ALLOCATABLE :: q1al(:), q2al(:), q3al(:)

   REAL(KIND=rp), ALLOCATABLE :: Xtraj(:, :), Ytraj(:, :), Ztraj(:, :)
   REAL(KIND=rp), ALLOCATABLE :: Vptraj(:, :), mutraj(:, :)
   REAL(KIND=rp), ALLOCATABLE :: betaFtraj(:, :), betaDtraj(:, :), FMtraj(:, :)
   REAL(KIND=rp), ALLOCATABLE :: Htraj(:, :), Pctraj(:, :)
   REAL(KIND=rp), ALLOCATABLE :: q1traj(:, :), q2traj(:, :), q3traj(:, :)
   REAL(KIND=rp), ALLOCATABLE :: ck1traj(:, :), ck2traj(:, :), ck3traj(:, :)

   REAL(KIND=rp)                  :: normB, Vstar1, Vstar2, Vstar3
   REAL(KIND=rp), DIMENSION(3, 3) :: gees

   !-----------------------------------------------------------------------------------------------------
   ! Transport quantities
   !-----------------------------------------------------------------------------------------------------
   REAL(KIND=rp), ALLOCATABLE :: Obs_GF(:, :), Obs_FF(:, :), Obs_FD(:, :), Obs_DF(:, :), Obs_DD(:, :)
   REAL(KIND=rp), ALLOCATABLE :: Vcorff(:, :), VcorTT(:, :), VcorTN(:, :), VcorNT(:, :)

   !-----------------------------------------------------------------------------------------------------
   ! Numerical and auxiliary variables
   !-----------------------------------------------------------------------------------------------------
   INTEGER :: i, k, re, lene, k2, i8
   INTEGER :: filled, ias, eror_flag
   INTEGER :: type_choice
   INTEGER :: arg_count, arg_type_choice, arg_sim_number, arg_status

   INTEGER, PARAMETER :: barWidth = 25

   REAL :: progress = 25.0

   REAL(KIND=rp), ALLOCATABLE :: t(:)
   REAL(KIND=rp), ALLOCATABLE :: Lagr_ref(:)
   REAL(KIND=rp), ALLOCATABLE :: Lagr_corr(:)

   REAL(KIND=rp) :: dt, dt_half
   INTEGER, PARAMETER :: nq = 51
   REAL(rp) :: quan(nq)
   REAL(rp), POINTER, CONTIGUOUS :: Qx(:), Qy(:), Qz(:), Qw(:), Qph(:)

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
      READ(arg, *, IOSTAT=arg_status) arg_type_choice
      IF (arg_status /= 0) THEN
         ERROR STOP 'First argument must be 1 (Sim) or 2 (DB).'
      END IF

      CALL GET_COMMAND_ARGUMENT(2, arg)
      READ(arg, *, IOSTAT=arg_status) arg_sim_number
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
      CALL load_run_parameters(run)

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

      CALL print_run_header(runs, nsim, date_time, folder, lene, simf, simi, 333)

!========================================================================================================
! VARIABLE ALLOCATION
!========================================================================================================

      ALLOCATE( X(Np), Xtraj(Nt + 1, ntraj), Obs_GF(10, Nt + 1), Obs_FF(10, Nt + 1), &
                Obs_FD(10, Nt + 1), Obs_DF(10, Nt + 1), Obs_DD(10, Nt + 1), &
                Vcorff(Nt + 1, Nt + 1), t(Nt + 1) )

      ALLOCATE( Y, Z, Vp, mu, Ham, F0, G0, FM, betaF, betaD, q1al, q2al, q3al, Pcc, &
                ck1, ck2, ck3, MOLD=X )

      ALLOCATE( Ytraj, Ztraj, Vptraj, mutraj, betaFtraj, betaDtraj, FMtraj, Htraj, q1traj, q2traj, &
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
         "position_type", "annulus_width", "pitch_type", "energy_type", "Lns", "Lts", &
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
         REAL(position_type, rp), annulus_width, REAL(pitch_type, rp), REAL(energy_type, rp), &
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
            CALL wavenum_unified(normB, Vstar1, Vstar2, Vstar3, gees, USE_real)
            CALL initialize_particles(X, Y, Z, mu, Vp, F0, G0, FM, betaF, betaD)
            CALL Larmor(USE_larmor, mu)

            IF (USE_testing == ON) THEN
               CALL testing_T3ST(X, Y, Z, Vp, mu, F0, G0, FM, betaF, betaD, Vstar1, Vstar2, Vstar3, &
                                 eror_flag)
               STOP
            END IF

            CALL write_initial_state(X, Y, Z, Vp, mu, betaF, betaD, FM)

!========================================================================================================
! TIME PROPAGATION
!========================================================================================================

            !$omp parallel default(shared) private(k, i)

            DO k = 1, Nt + 1

               !$omp single
                quan = 0.0_rp
               !$omp end single

               !$omp do schedule(static) private(Qx, Qy, Qz, Qw, Qph) reduction(+:quan)

               DO i = 1, Np
                  BLOCK

                     REAL(rp) :: xi, yi, zi, mui, betaFi, betaDi, vpi
                     REAL(rp) :: mask, maskR, maskZ
                     REAL(rp) :: vx, vy, vz, vm, ap, vbF, vbD
                     REAL(rp) :: Wx, Wy, Wz, Wm, Wp, WbF, WbD
                     REAL(rp) :: q1, q2, q3, Vrt, Vrr, FMi, F0i, G0i
                     REAL(rp) :: pbGF, pbFF, pbFD, pbDF, pbDD
                     REAL(rp) :: Hi, B, Vtx, VFx, sm1, sm2, Pc, kin
                     REAL(rp) :: check_1, check_2, check_3

                     ! Select per-particle or single-spectrum arrays.
                     IF (USE_real == OFF) THEN
                        Qx  => kx(:, i)
                        Qy  => ky(:, i)
                        Qz  => kz(:, i)
                        Qw  => w(:, i)
                        Qph => ph(:, i)
                     ELSE IF (USE_real == ON) THEN
                        Qx  => kxs
                        Qy  => kys
                        Qz  => kzs
                        Qw  => ws
                        Qph => phs
                     END IF

                     ! Load to scalars.
                     xi      = X(i)
                     yi      = Y(i)
                     zi      = Z(i)
                     mui     = mu(i)
                     vpi     = Vp(i)
                     F0i     = F0(i)
                     G0i     = G0(i)
                     betaFi  = betaF(i)
                     betaDi  = betaD(i)

                     ! Half-collisional Euler-Maruyama step.
                     ! This keeps the Stratonovich interpretation in play.
                     IF ((USE_coll == ON) .AND. (t(k) > tc)) THEN

                        sm1 = rng_uniform()
                        sm2 = rng_uniform()

                        CALL gyrocenter_drifts(xi, yi, zi, vpi, mui, q1, q2, q3, t(k), &
                                    vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx, &
                                    check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph)

                        CALL collision_kicks(sm1, sm2, dt_half, xi, yi, zi, vpi, mui, &
                                            B, vx, vy, vz, vm, ap)

                        xi      = xi      + dt_half * vx
                        yi      = yi      + dt_half * vy
                        zi      = zi      + dt_half * vz
                        vpi     = vpi     + dt_half * ap
                        mui     = mui     + dt_half * vm
                        betaFi  = betaFi  + dt_half * vbF
                        betaDi  = betaDi  + dt_half * vbD

                     END IF

                     ! RK4 stage 1.
                     CALL gyrocenter_drifts(xi, yi, zi, vpi, mui, q1, q2, q3, t(k), &
                                 vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx, &
                                 check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph)

                     IF (k == 1) Lagr_ref(i) = Vtx

                     Wx = vx / 6.0_rp
                     Wy = vy / 6.0_rp
                     Wz = vz / 6.0_rp
                     Wm = vm / 6.0_rp
                     Wp = ap / 6.0_rp
                     WbF = vbF / 6.0_rp
                     WbD = vbD / 6.0_rp

                     ! Mask for particles outside the domain.
                     mask = 1.0_rp

                     IF (ANY(magnetic_model == [1, 2])) THEN
                        maskR = (xi - minR) * (maxR - xi)
                        maskZ = (yi - minZ) * (maxZ - yi)
                        mask  = MERGE(1.0_rp, 0.0_rp, (maskR >= 0.0) .AND. (maskZ >= 0.0))
                     END IF

                     ! Observables.
                     IF (.NOT. ISNAN(q1)) THEN
                        pbGF = 1.0_rp/Np
                        pbFF = FMi*exp(betaFi)/G0i/Np
                        pbFD = FMi*exp(betaDi)/G0i/Np
                        pbDF = FMi*(exp(betaFi)-1.0_rp)/G0i/Np
                        pbDD = FMi*(exp(betaDi)-1.0_rp)/G0i/Np
                        
                        ! Obs columns: Vtx1, Q1, Vp1, H1, Heat1, Vtx2, Q2, Vp2, H2, Heat2.
                        quan(1:10)   = quan(1:10)   + pbGF*mask * [Vtx, q1, vpi, Hi, Vtx*kin, Vtx**2, q1**2, vpi**2, Hi**2, Vtx*kin] !green function
                        quan(11:20)  = quan(11:20)  + pbFF*mask * [Vtx, q1, vpi, Hi, Vtx*kin, Vtx**2, q1**2, vpi**2, Hi**2, Vtx*kin] !fullF + full V
                        quan(21:30)  = quan(21:30)  + pbFD*mask * [Vtx, q1, vpi, Hi, Vtx*kin, Vtx**2, q1**2, vpi**2, Hi**2, Vtx*kin] !fullF + delta V
                        quan(31:40)  = quan(31:40)  + pbDF*mask * [Vtx, q1, vpi, Hi, Vtx*kin, Vtx**2, q1**2, vpi**2, Hi**2, Vtx*kin] !deltaF + full V
                        quan(41:50)  = quan(41:50)  + pbDD*mask * [Vtx, q1, vpi, Hi, Vtx*kin, Vtx**2, q1**2, vpi**2, Hi**2, Vtx*kin] !deltaF + delta V
                        quan(51)     = quan(51)     + pbGF*mask * (Lagr_ref(i) * Vtx)    ! green function

                     END IF
!                    IF (i == 12) THEN
 !                     WRITE(*, '("particle | ",I0,"| FM*Exp/F0 | ",G0.6,"| FM/F0-1.0 | ",G0.6,"| beta | ",G0.6,"| dbeta/dt | ",G0.6)') &
  !                       i, 100.0*(FMi*exp(betaFi)/F0i-1.0), 100.0*(FMi/F0i-1.0_rp), betaFi, vbF/Vtx
   !               END IF
                ! RK4 stage 2.
                     CALL gyrocenter_drifts(xi + vx * dt / 2.0_rp, &
                                     yi + vy * dt / 2.0_rp, &
                                     zi + vz * dt / 2.0_rp, &
                                     vpi + ap * dt / 2.0_rp, &
                                     mui + vm * dt / 2.0_rp, &
                                     q1, q2, q3, t(k) + dt / 2.0_rp, &
                                     vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx, &
                                     check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph)

                     Wx = Wx + vx / 3.0_rp
                     Wy = Wy + vy / 3.0_rp
                     Wz = Wz + vz / 3.0_rp
                     Wm = Wm + vm / 3.0_rp
                     Wp = Wp + ap / 3.0_rp
                     WbF = WbF + vbF / 3.0_rp
                     WbD = WbD + vbD / 3.0_rp

                     ! RK4 stage 3.
                     CALL gyrocenter_drifts(xi + vx * dt / 2.0_rp, &
                                     yi + vy * dt / 2.0_rp, &
                                     zi + vz * dt / 2.0_rp, &
                                     vpi + ap * dt / 2.0_rp, &
                                     mui + vm * dt / 2.0_rp, &
                                     q1, q2, q3, t(k) + dt / 2.0_rp, &
                                     vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx, &
                                     check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph)

                     Wx = Wx + vx / 3.0_rp
                     Wy = Wy + vy / 3.0_rp
                     Wz = Wz + vz / 3.0_rp
                     Wm = Wm + vm / 3.0_rp
                     Wp = Wp + ap / 3.0_rp
                     WbF = WbF + vbF / 3.0_rp
                     WbD = WbD + vbD / 3.0_rp

                     ! RK4 stage 4.
                     CALL gyrocenter_drifts(xi + vx * dt, yi + vy * dt, zi + vz * dt, &
                                 vpi + ap * dt, mui + vm * dt, &
                                 q1, q2, q3, t(k) + dt, &
                                 vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx, &
                                 check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph)

                     Wx = Wx + vx / 6.0_rp
                     Wy = Wy + vy / 6.0_rp
                     Wz = Wz + vz / 6.0_rp
                     Wm = Wm + vm / 6.0_rp
                     Wp = Wp + ap / 6.0_rp
                     WbF = WbF + vbF / 6.0_rp
                     WbD = WbD + vbD / 6.0_rp

                     ! Commit back.
                     X(i)      = xi      + dt * Wx
                     Y(i)      = yi      + dt * Wy
                     Z(i)      = zi      + dt * Wz
                     Vp(i)     = vpi     + dt * Wp
                     mu(i)     = mui     + dt * Wm
                     betaF(i)  = betaFi  + dt * WbF
                     betaD(i)  = betaDi  + dt * WbD

                     Ham(i)  = Hi
                     FM(i)   = FMi
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
                        betaFi  = betaF(i)
                        betaDi  = betaD(i)

                        sm1 = rng_uniform()
                        sm2 = rng_uniform()

                        CALL gyrocenter_drifts(xi, yi, zi, vpi, mui, q1, q2, q3, t(k), &
                                    vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx, &
                                    check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph)

                        CALL collision_kicks(sm1, sm2, dt_half, xi, yi, zi, vpi, mui, &
                                             B, vx, vy, vz, vm, ap)

                        X(i)      = X(i)      + dt_half * vx
                        Y(i)      = Y(i)      + dt_half * vy
                        Z(i)      = Z(i)      + dt_half * vz
                        Vp(i)     = Vp(i)     + dt_half * ap
                        mu(i)     = mu(i)     + dt_half * vm
                        betaF(i)  = betaF(i)  + dt_half * vbF
                        betaD(i)  = betaD(i)  + dt_half * vbD

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
               Lagr_corr(k) = quan(51)

               Xtraj(k, :)  = X(1:ntraj)
               Ytraj(k, :)  = Y(1:ntraj)
               Ztraj(k, :)  = Z(1:ntraj)
               Vptraj(k, :) = Vp(1:ntraj)
               mutraj(k, :) = mu(1:ntraj)
               betaFtraj(k, :) = betaF(1:ntraj)
               betaDtraj(k, :) = betaD(1:ntraj)
               FMtraj(k, :) = FM(1:ntraj)
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
            CALL write_trajectories(Xtraj, Ytraj, Ztraj, Vptraj, mutraj, betaFtraj, betaDtraj, FMtraj, Htraj, &
                                    q1traj, q2traj, q3traj, Pctraj, ck1traj, ck2traj, ck3traj)

            CALL write_transport(Obs_GF, Obs_FF, Obs_FD, Obs_DF, Obs_DD)
            CALL write_final_state(X, Y, Z, Vp, mu, betaF, betaD, FM, Ham, q1al, q2al)
            CALL write_lagrangian(Lagr_corr)

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

      DEALLOCATE( X, Y, Z, Vp, mu, Ham, F0, G0, FM, betaF, betaD, Obs_GF, Obs_FF, Obs_FD, Obs_DF, Obs_DD, &
                  Xtraj, Ytraj, Ztraj, Vptraj, mutraj, betaFtraj, betaDtraj, FMtraj, Htraj, &
                  q1traj, q2traj, q3traj, Vcorff, VcorTT, VcorTN, VcorNT, &
                  t, q1al, q2al, q3al, Pcc, ck1, ck2, ck3, Pctraj, &
                  ck1traj, ck2traj, ck3traj, Lagr_ref, Lagr_corr )

      IF ((magnetic_model == 1) .OR. (magnetic_model == 2)) THEN
         DEALLOCATE(Efit_data)
      END IF

   END DO

   CLOSE(333)

END PROGRAM T3ST
