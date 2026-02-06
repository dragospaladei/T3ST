!========================================================================================================
! FAST2 / T3ST :: Turbulent Transport in Tokamaks from Statistical Trajectories
!
! Purpose  :: test-particle (DNS) simulations of ion trajectories in turbulent tokamak-like environments
!          :: Tokamak-like MHD equilibrium (realistic or modeled)
!          :: Turbulence via ensemble of random fields with given statistical properties (ITG / TEM)
!          :: Linear analytical dispersion relation; no magnetic perturbations
!          :: For each random field, compute a trajectory from given initial conditions
!          :: Equations of motion consistent with GK gyro-centers
!          :: Projected onto cylindrical coordinates (X,Y,Z)
!          :: Turbulence in field-aligned (q1,q2,q3) ~ (psi, zeta-nu, nu)
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
   USE constants
   USE random_numbers
   USE omp_lib
   USE drift_kernels         !, only: Drift2
   USE testing_T3ST_mod
   USE rng_mod
   USE, INTRINSIC :: IEEE_ARITHMETIC
   IMPLICIT NONE

!========================================================================================================
! VARIABLE DECLARATION
!========================================================================================================

!----------------- folders & sims -----------------------------------------------------------------------
   CHARACTER(LEN=10000)                      :: CWD          ! current working directory (for exporting data)
   INTEGER                                   :: lenCWD       ! length of CWD
   CHARACTER(LEN=:), ALLOCATABLE             :: address      ! base export directory
   CHARACTER(LEN=35)                         :: runs, nsim   ! run index -> string, and total sims -> string
   INTEGER                                   :: run          ! simulation index
   INTEGER                                   :: simi, simf   ! first/last simulation indices
   CHARACTER(LEN=:), ALLOCATABLE             :: folder       ! per-run output folder

!----------------- CPU time ------------------------------------------------------------------------------
   INTEGER                                   :: count0, count1, count2, count3
   INTEGER                                   :: count_rate, count_max
   CHARACTER(LEN=18)                         :: date_time

!----------------- random seeds -------------------------------------------------------------------------
   INTEGER, ALLOCATABLE, DIMENSION(:)        :: seed
   INTEGER                                   :: ne, clock, ja

!----------------- trajectory variables -----------------------------------------------------------------
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: X, Y, Z, Vp, mu, Ham, Pb, Pcc, ck1, ck2, ck3
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: q1al, q2al, q3al

!   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: Xtot, Ytot, Ztot, Vptot, mutot, Htot
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: q1tot, q2tot, q3tot

   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: Xtr, Ytr, Ztr, Vptr, mutr, Htr, Pctr
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: q1tr, q2tr, q3tr, ck1tr, ck2tr, ck3tr

   REAL(KIND=dp)                             :: normB, Vstar1, Vstar2, Vstar3
   REAL(KIND=dp), DIMENSION(3, 3)            :: gees

!----------------- transport quantities -----------------------------------------------------------------
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: vit, dif
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: Vcorff, VcorTT, VcorTN, VcorNT

!----------------- numerical / auxiliaries --------------------------------------------------------------
   INTEGER                                   :: i, j, k, m, re, lene, k2, i8
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: t
   INTEGER                                   :: filled, ias, eror_flag
   REAL                                      :: progress = 25.0
   INTEGER, PARAMETER                        :: barWidth = 25
   CHARACTER(LEN=80)                         :: bar
   LOGICAL                                   :: has_nan, ok

   REAL(KIND=dp)                             :: dt1, dt2, dt, dt_half
   REAL(KIND=dp)                             :: quan11, quan12, quan13, quan14, quanVtxCorr
   REAL(KIND=dp)                             :: quan21, quan22, quan23, quan24
   REAL(KIND=dp), ALLOCATABLE                :: Lagr_ref(:)     ! Vtx at k=1, per particle
   REAL(KIND=dp), ALLOCATABLE                :: Lagr_corr(:)    ! correlation vs time, size Nt+1

   REAL(dp), POINTER, CONTIGUOUS             :: Qx(:), Qy(:), Qz(:), Qw(:), Qph(:), QL(:), Q0wrp(:)

!========================================================================================================
! ADDRESS DECLARATION FOR EXPORTING DATA
!========================================================================================================
   CALL sim_values

   CALL GetCWD(CWD)
   address = TRIM(ADJUSTL(CWD))
   lenCWD  = LEN(address)

   ! Must have 'source' and 'data' folders at same level.
   address = address(1:lenCWD - 7)//'/data/Sim_'//noofsimstr//'/'  

   OPEN (333, FILE=address//'/logs.dat')     ! logs of simulations

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
      CALL SYSTEM_CLOCK(count0, count_rate, count_max)
      CALL parameters(run)

!========================================================================================================
! FOLDER DECLARATION
!========================================================================================================
      WRITE (runs, *) INT(run)

      IF (run < 10) THEN
         folder = address//'/Run_0000'//TRIM(ADJUSTL(runs))//'/'
      ELSE IF (run < 100) THEN
         folder = address//'/Run_000'//TRIM(ADJUSTL(runs))//'/'
      ELSE IF (run < 1000) THEN
         folder = address//'/Run_00'//TRIM(ADJUSTL(runs))//'/'
      ELSE IF (run < 10000) THEN
         folder = address//'/Run_0'//TRIM(ADJUSTL(runs))//'/'
      ELSE
         folder = address//'/Run_'//TRIM(ADJUSTL(runs))//'/'
      END IF

      ! Create data files: "...\folder\file_XX.*"
      OPEN (11, FILE=folder//'file_01.dat')   ! parameters export

      DO i = 2, 30
         WRITE (runs, *) i
         IF (i < 10) THEN
            OPEN (UNIT=i+10, FILE=folder//'file_0'//TRIM(ADJUSTL(runs))//'.bin', &
                  FORM='unformatted', ACCESS='stream', STATUS='replace')
         ELSE
            OPEN (UNIT=i+10, FILE=folder//'file_'//TRIM(ADJUSTL(runs))//'.bin', &
                  FORM='unformatted', ACCESS='stream', STATUS='replace')
         END IF
      END DO


      ! Console print: simulation number and folder
      WRITE (runs, *) run
      WRITE (nsim, *) (simf - simi + 1)
      lene = LEN(folder)
      CALL printing(runs, nsim, date_time, folder, lene, simf, simi, 333)

!========================================================================================================
! VARIABLE ALLOCATION
!========================================================================================================
      ! Allocate (and deallocate) variables for each run (depending on Np/Nt).
      ALLOCATE (X(Np), Xtr(Nt + 1, ntraj), vit(4, Nt + 1), dif(4, Nt + 1), &
                Vcorff(Nt + 1, Nt + 1), t(Nt + 1))
      ALLOCATE (Y, Z, Vp, mu, Ham, Pb, q1al, q2al, q3al, Pcc, ck1,ck2,ck3,MOLD=X)
      ALLOCATE (Ytr, Ztr, Vptr, mutr, Htr, q1tr, q2tr, q3tr, Pctr, ck1tr, ck2tr, ck3tr, MOLD=Xtr)
      ALLOCATE (VcorTT, VcorTN, VcorNT, MOLD=Vcorff)
      ALLOCATE( Lagr_ref(Np), Lagr_corr(Nt+1) )
!      ALLOCATE (randm1(Nt+1, Np), randm2(Nt+1, Np))
!      ALLOCATE ( Xtot(Np, 2))
!      ALLOCATE (Ytot, Ztot, Vptot, mutot, Htot, q1tot, q2tot, q3tot, MOLD=Xtot)

!========================================================================================================
! EXPORT SOME PARAMETERS (names & values) OF THIS RUN
!========================================================================================================
      WRITE(11,*) [ CHARACTER(LEN=15) :: "t0", "tc", "tt", "tmax", "Np", "Nc", "Nloop", "Nreal", "Nt", "ntraj", &
                   "Nci", "Nce", "Ti", "Te", "B0", "vth", "rhoi", "R0", "Ln", "Li", "Le", "Zeff", "Aeff", "ndens", &
                   "a0", "Omgt0", "Omgtprim", "wi", "C1", "C2", "C3", "safe1", "safe2", "safe3", &
                   "amp", "elong", "device", "shot", "shotslice", "NgridR", "NgridZ", "Phi", "turbprof", "Ai", &
                   "Ae", "lambdax", "lambday", "lambdaz", "lbalonz", "tauc", "k0i", "k0e", "gamma_ZF", "gamma_E","x_corr","t_corr", "X0", "Y0", "Z0", "r00", "q00", &
                   "psi0", "Tw", "Ew", "pitch", "Aw", "Zw", "taucc", "magnetic_model", "USE_larmor", "position_type", &
                   "pitch_type", "energy_type", "USE_coll", "USE_turb", "USE_magnturb", "USE_freq", "USE_polar", &
                   "USE_PC", "USE_real", "USE_corr", "USE_balloon", "USE_tilt", "USE_testing"]

      WRITE(11,*) [ t0, tc, tt, tmax, REAL(Np,dp), REAL(Nc,dp), REAL(Nloop,dp), REAL(Nreal,dp), REAL(Nt,dp), &
                   REAL(ntraj,dp), REAL(Nci,dp), REAL(Nce,dp), Ti, Te, B0, vth, rhoi, R0, Ln, Li, Le, Zeff, Aeff, ndens, &
                   a0, Omgt0, Omgtprim, wi, C1, C2, C3, safe1, safe2, safe3, amp, elong, REAL(device,dp), &
                   REAL(shot,dp), REAL(shotslice,dp), REAL(NgridR,dp), REAL(NgridZ,dp), Phi, turbprof, Ai, Ae, &
                   lambdax, lambday, lambdaz, lbalonz, tauc, k0i, k0e, gamma_ZF, gamma_E, REAL(x_corr,dp), REAL(t_corr,dp), X0, Y0, Z0, r00, q00, psi0, Tw, Ew, pitch, Aw, Zw, &
                   taucc, REAL(magnetic_model, dp), REAL(USE_larmor, dp), REAL(position_type, dp), REAL(pitch_type, dp), &
                   REAL(energy_type, dp), REAL(USE_coll, dp), REAL(USE_turb, dp), REAL(USE_magnturb, dp), REAL(USE_freq, dp), &
                   REAL(USE_polar, dp), REAL(USE_PC, dp), REAL(USE_real, dp), REAL(USE_corr, dp), REAL(USE_balloon,dp), REAL(USE_tilt,dp), REAL(USE_testing,dp)]

!========================================================================================================
! TIME GRID
!========================================================================================================
      dt = (tmax - t0) / REAL(Nt)
      dt_half = dt/2.0_dp
      t  = [(i8*dt + t0, i8=0, Nt)]
      Lagr_ref  = 0.0_dp
      Lagr_corr = 0.0_dp


!========================================================================================================
! REALIZATIONS LOOP
!========================================================================================================
      DO re = 1, Nreal
         DO k2 = 1, Nloop

            CALL dispersion(normB, Vstar1, Vstar2, Vstar3, gees)

            IF (USE_real == OFF) THEN
               CALL wavenum(normB, Vstar1, Vstar2, Vstar3, gees)
               CALL initial_conditions(X, Y, Z, mu, Vp, pb)
               CALL Larmor(USE_larmor, X, Y, Z, mu)
            ELSEIF (USE_real == ON) THEN
               CALL wavenum_single(normB, Vstar1, Vstar2, Vstar3, gees)
               CALL initial_conditions(X, Y, Z, mu, Vp, pb)
               CALL Larmor_single(USE_larmor, X, Y, Z, mu)
            END IF
            
    	    IF(USE_testing == ON) then
	       CALL testing_T3ST(X,Y,Z,Vp,mu,pb,Vstar1,Vstar2,Vstar3, eror_flag)
	       STOP
	    END IF

            WRITE (23) X
            WRITE (24) Y

            !============================================================================================
            ! TIME PROPAGATION
            !============================================================================================

            !$omp parallel default(shared) private(k,i)
            DO k = 1, Nt + 1

               !$omp single
               quan11 = 0.0_dp; quan12 = 0.0_dp; quan13 = 0.0_dp; quan14 = 0.0_dp
               quan21 = 0.0_dp; quan22 = 0.0_dp; quan23 = 0.0_dp; quan24 = 0.0_dp; quanVtxCorr = 0.0_dp
               !$omp end single

               !$omp do schedule(static) private(Qx,Qy,Qz,Qw,Qph,QL,Q0wrp) &
               !$omp   reduction(+: quan11,quan12,quan13,quan14, quan21,quan22,quan23,quan24,quanVtxCorr)
               DO i = 1, Np
                  BLOCK
                     REAL(dp) :: xi, yi, zi, mui, vpi, pbi, mask, maskR, maskZ
                     REAL(dp) :: vx, vy, vz, vm, ap, Wx, Wy, Wz, Wm, Wp
                     REAL(dp) :: q1, q2, q3, Vrt, Vrr, temp
                     REAL(dp) :: Hi, B, Vtx, Vty, sm1, sm2, Pc, check_1, check_2, check_3!, dt_half

                     ! Select per-particle or single-spectrum arrays
                     IF (USE_real == OFF) THEN
                        Qx  => kx(:, i);   Qy  => ky(:, i);   Qz  => kz(:, i)
                        Qw  => w(:, i);    Qph => ph(:, i);   QL  => L(:, i);    Q0wrp => q00wrap(:,i)
		     ELSEIF (USE_real == ON) THEN
                        Qx  => kxs;        Qy  => kys;        Qz  => kzs
                        Qw  => ws;         Qph => phs;        QL  => Ls(:, i);   Q0wrp => q00wraps
                     END IF


                     ! Load to scalars (helps vectorizers)
                     xi  = X(i);  yi  = Y(i);  zi  = Z(i)
                     mui = mu(i); vpi = Vp(i); pbi = pb(i)

		     ! Half_collisional EulerMaruyama step
                     ! Note: we do here a half EM step in order to keep the Stratanovich interpretation in play
		     IF ((USE_coll == ON).and.(t(k).gt.tc)) then
                        sm1 = rng_uniform()
                        sm2 = rng_uniform()
                        CALL Drift2(dt, xi, yi, zi, vpi, mui, q1, q2, q3, t(k), &
                                    vx, vy, vz, ap, vm, Hi, Pc, B, Vtx, Vty, check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph, QL, Q0wrp)
                        CALL collisions_1_MP(sm1, sm2, dt_half, xi, yi, zi, vpi, mui, B, vx, vy, vz, vm, ap)

                        xi  = xi  + dt_half*vx
                        yi  = yi  + dt_half*vy
                        zi  = zi  + dt_half*vz
                        vpi = vpi + dt_half*ap
                        mui = mui + dt_half*vm
                     END IF

                     ! RK4 stage 1
                     CALL Drift2(dt, xi, yi, zi, vpi, mui, q1, q2, q3, t(k), &
                                 vx, vy, vz, ap, vm, Hi, Pc, B, Vtx, Vty, check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph, QL, Q0wrp)
                     ! store reference Vtx at k=1 (per particle)
                     IF (k == 1) Lagr_ref(i) = Vtx

                     Wx = vx/6.0_dp; Wy = vy/6.0_dp; Wz = vz/6.0_dp
                     Wm = vm/6.0_dp; Wp = ap/6.0_dp

                     ! Mask for particles outside the domain
                     mask = 1.0_dp
                     IF (ANY(magnetic_model == [1, 2])) THEN
                        maskR = (xi - minR) * (maxR - xi)
                        maskZ = (yi - minZ) * (maxZ - yi)
                        mask  = MERGE(1.0_dp, 0.0_dp, (maskR >= 0.0) .AND. (maskZ >= 0.0))
                     END IF

                     ! Observables (ignore NaNs)
                     IF (.NOT. ISNAN(q1)) THEN
                        quan11 = quan11 + pbi*mask*q1           ! < x >   (x coordinate, q1)
                        quan12 = quan12 + pbi*mask*vpi          ! < vp >  (parallel velocity)
                        quan13 = quan13 + pbi*mask*Hi           ! < E >   (energy)
                        quan14 = quan14 + pbi*mask*Vtx          ! < Vtx > (turbulent x drift)

                        quan21 = quan21 + pbi*mask*q1**2        ! < x^2 >   (x coordinate, q1)
                        quan22 = quan22 + pbi*mask*vpi**2       ! < vp^2 >  (parallel velocity)
                        quan23 = quan23 + pbi*mask*Hi**2        ! < E^2 >   (energy)
                        quan24 = quan24 + pbi*mask*Vtx**2       ! < Vtx^2 > (turbulent x drift)
	                quanVtxCorr = quanVtxCorr + pbi*mask * (Lagr_ref(i) * Vtx)    ! Lagrangian correlation <v_{q1}(0)v_{q1}(t)>
                     END IF
                     ! RK4 stage 2
                     CALL Drift2(dt, xi+vx*dt/2.0_dp, yi+vy*dt/2.0_dp, zi+vz*dt/2.0_dp, &
                                 vpi+ap*dt/2.0_dp, mui+vm*dt/2.0_dp, q1, q2, q3, t(k)+dt/2.0_dp, &
                                 vx, vy, vz, ap, vm, Hi, Pc, B, Vtx, Vty, check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph, QL, Q0wrp)
                     Wx = Wx + vx/3.0_dp; Wy = Wy + vy/3.0_dp; Wz = Wz + vz/3.0_dp
                     Wm = Wm + vm/3.0_dp; Wp = Wp + ap/3.0_dp

                             ! RK4 stage 3
                     CALL Drift2(dt, xi+vx*dt/2.0_dp, yi+vy*dt/2.0_dp, zi+vz*dt/2.0_dp, &
                                 vpi+ap*dt/2.0_dp, mui+vm*dt/2.0_dp, q1, q2, q3, t(k)+dt/2.0_dp, &
                                 vx, vy, vz, ap, vm, Hi, Pc, B, Vtx, Vty, check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph, QL, Q0wrp)
                     Wx = Wx + vx/3.0_dp; Wy = Wy + vy/3.0_dp; Wz = Wz + vz/3.0_dp
                     Wm = Wm + vm/3.0_dp; Wp = Wp + ap/3.0_dp

                     ! RK4 stage 4
                     CALL Drift2(dt, xi+vx*dt, yi+vy*dt, zi+vz*dt, vpi+ap*dt, mui+vm*dt, &
                                 q1, q2, q3, t(k)+dt, vx, vy, vz, ap, vm, Hi, Pc, B, Vtx, Vty, check_1, check_2, check_3, &
                                 Qx, Qy, Qz, Qw, Qph, QL, Q0wrp)
                     Wx = Wx + vx/6.0_dp; Wy = Wy + vy/6.0_dp; Wz = Wz + vz/6.0_dp
                     Wm = Wm + vm/6.0_dp; Wp = Wp + ap/6.0_dp

                     ! Commit back
                     X(i)   = xi  + dt*Wx
                     Y(i)   = yi  + dt*Wy
                     Z(i)   = zi  + dt*Wz
                     Vp(i)  = vpi + dt*Wp
                     mu(i)  = mui + dt*Wm
                     Ham(i) = Hi
                     Pcc(i) = Pc
                     q1al(i) = q1
                     q2al(i) = q2
                     q3al(i) = q3
                     ck1(i) = check_1
                     ck2(i) = check_2
                     ck3(i) = check_3
                     
                     ! Half_collisional EulerMaruyama step
                     ! Note: we do here the other half EM step in order to keep the Stratanovich interpretation in play;
                     !       thus, with the next loop step it combines into a single step operation
		     IF ((USE_coll == ON).and.(t(k).gt.tc)) then
                     xi  = X(i);  yi  = Y(i);  zi  = Z(i)
                     mui = mu(i); vpi = Vp(i); pbi = pb(i)
                        sm1 = rng_uniform()
                        sm2 = rng_uniform()
                        CALL Drift2(dt, xi, yi, zi, vpi, mui, q1, q2, q3, t(k), &
                                    vx, vy, vz, ap, vm, Hi,Pc, B, Vtx, Vty, check_1, check_2, check_3, Qx, Qy, Qz, Qw, Qph, QL, Q0wrp)
                        CALL collisions_1_MP(sm1, sm2, dt_half, xi, yi, zi, vpi, mui, B, vx, vy, vz, vm, ap)

                        X(i)  = X(i)  + dt_half*vx
                        Y(i)  = Y(i)  + dt_half*vy
                        Z(i)  = Z(i)  + dt_half*vz
                        vp(i) = Vp(i) + dt_half*ap
                        mu(i) = mu(i) + dt_half*vm
                     END IF

                  END BLOCK
               END DO
               !$omp end do

               !$omp single
               vit(1, k) = quan11/Np;  vit(2, k) = quan12/Np;  vit(3, k) = quan13/Np;  vit(4, k) = quan14/Np
               dif(1, k) = quan21/Np;  dif(2, k) = quan22/Np;  dif(3, k) = quan23/Np;  dif(4, k) = quan24/Np; Lagr_corr(k) = quanVtxCorr / Np

               Xtr(k, :)  = X(1:ntraj)
               Ytr(k, :)  = Y(1:ntraj)
               Ztr(k, :)  = Z(1:ntraj)
               Vptr(k, :) = Vp(1:ntraj)
               mutr(k, :) = mu(1:ntraj)
               Htr(k, :)  = Ham(1:ntraj)
               Pctr(k, :) = Pcc(1:ntraj)

	       ck1tr(k, :) = ck1(1:ntraj)
               ck2tr(k, :) = ck2(1:ntraj)
               ck3tr(k, :) = ck3(1:ntraj)

               q1tr(k, :) = q1al(1:ntraj)
               q2tr(k, :) = q2al(1:ntraj)
               q3tr(k, :) = q3al(1:ntraj)
               !$omp end single

               !$OMP BARRIER
            END DO
            !$omp end parallel

            !============================================================================================
            ! EXPORT DATA
            !============================================================================================
            WRITE (12) Xtr
            WRITE (13) Ytr
            WRITE (14) Ztr
            WRITE (15) Vptr
            WRITE (16) mutr
            WRITE (17) Htr
            WRITE (18) q1tr
            WRITE (19) q2tr
            WRITE (20) q3tr

            WRITE (21) vit
            WRITE (22) dif
            WRITE (25) X
            WRITE (26) Y
            WRITE (27) Ham
            WRITE (28) q1al
            WRITE (29) q2al

            WRITE (30) Vcorff
            WRITE (31) VcorTT
            WRITE (32) VcorTN
            WRITE (33) VcorNT

            WRITE (35) Pctr
            WRITE (36) ck1tr
            WRITE (37) ck2tr
            WRITE (38) ck3tr

            WRITE (39) Lagr_corr

           IF (ANY(magnetic_model == [1, 2])) THEN
               WRITE (34) REAL(COUNT( ((X - minR)*(maxR - X) >= 0.0_dp) .AND. &
                                     ((Y - minZ)*(maxZ - Y) >= 0.0_dp) .AND. &
                                     .NOT. IEEE_IS_NAN(X) .AND. .NOT. IEEE_IS_NAN(Y) ), dp)
            ELSEIF (ANY(magnetic_model == [3, 4])) THEN
               WRITE (34) REAL(COUNT( ((X - (1.0-a0))*((1.0+a0) - X) >= 0.0_dp) .AND. &
                                     ((Y - (-a0))*((a0) - Y) >= 0.0_dp) .AND. &
                                     .NOT. IEEE_IS_NAN(X) .AND. .NOT. IEEE_IS_NAN(Y) ), dp)
            END IF

            ! Progress bar (per super-ensemble loop)
            progress = REAL(k2*re) / REAL(Nloop*Nreal)
            filled   = INT(progress * barWidth)

            bar = REPEAT(' ', LEN(bar))
            DO ias = 1, barWidth
               IF (ias <= filled) THEN
                  bar(ias:ias) = '#'
               ELSE
                  bar(ias:ias) = '-'
               END IF
            END DO

            WRITE(*,'(1X,I6,"/",I6," realizations & ",I6,"/",I6," loops [",A,"] ",F6.1,"%",A)', ADVANCE='yes') &
                 re, Nreal, k2, Nloop, bar(1:barWidth), progress*100.0, ACHAR(13)
            CALL flush(6)

         END DO

         PRINT *, ''   ! newline after completion of loops
      END DO  ! Nreal

!========================================================================================================
! CLOSE FILES / TIMING / DEALLOCATIONS
!========================================================================================================
      DO i = 1, 30
         CLOSE (i + 10)
      END DO

      CALL SYSTEM_CLOCK(count3, count_rate, count_max)
      WRITE (*, *) ' '
      WRITE (*, *) 'Total computation time  = ', (1.0*(count3 - count0)) / (1.0*count_rate)
      WRITE (*, *) '------------------------------------------------------------------------------'
      WRITE (333, *) 'Total computation time  = ', (1.0*(count3 - count0)) / (1.0*count_rate)
      WRITE (333, *) '------------------------------------------------------------------------------'

      DEALLOCATE( X, Y, Z, Vp, mu, Ham, Pb, vit, dif, Xtr, Ytr, Ztr, Vptr, mutr, Htr, q1tr, q2tr, q3tr, &
                  Vcorff, VcorTT, VcorTN, VcorNT, t, q1al, q2al, q3al,  Pcc, ck1,ck2,ck3, Pctr, ck1tr, ck2tr, ck3tr)
!      DEALLOCATE(Xtot, Ytot, Ztot, Vptot, q1tot, q2tot, q3tot, mutot, Htot,randm1, randm2)
  
      IF ((magnetic_model == 1) .OR. (magnetic_model == 2)) THEN
         DEALLOCATE (Efit_data)
      END IF

   END DO  ! run loop

   CLOSE (333)

END PROGRAM T3ST

