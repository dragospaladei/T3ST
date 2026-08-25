!  ========================================================================================================================================================
!      ! File       :: Particle_initial_conditions.f90
!      ! Purpose    :: Particle energy, pitch, position, weight, and initial-field setup
!
!      ! Subroutine :: initialize_particles
!      ! Purpose    :: fixes the initial GY coordinates of particles. Consequently their energetic distribution
!  ========================================================================================================================================================
!      ! Record of revisions:
!      !    Date       |   Programmer           |   Description of change
!      !  =======      |   =============        |   =========================
!      !  --/--/2021   |   D. I. Palade         |   Initial version
!      !  --/--/2023   |   L. M. Pomarjanschi   |   Corections and development
!      !  01/01/2024   |   D. I. Palade         |   Remake
!  ========================================================================================================================================================

SUBROUTINE initialize_particles(X, Y, Z, mu, Vp, Einit, F0, G0, FM, betaF, betaD, alpha_1, alpha_2, alpha_3)
   USE constants
   USE random_numbers
   IMPLICIT NONE

   ! I/O variables
   REAL(KIND=wp), DIMENSION(Np), INTENT(OUT) :: X, Y, Z, mu, Vp, Einit, F0, G0, FM, betaF, betaD
   REAL(KIND=wp), DIMENSION(Np), INTENT(OUT) :: alpha_1, alpha_2, alpha_3

   ! Local variables
   REAL(KIND=wp), DIMENSION(Np) :: En, PA, aux, Jac
   REAL(KIND=wp), DIMENSION(Np) :: phi0, phi1, phi2, B, Bx, By, Bz, q1, q2, q3
   REAL(KIND=wp), DIMENSION(Np) :: rhot, delta_q1, temp, dens
   REAL(KIND=wp)                :: m1, m2, Ts_eff

   ! ------------------------------------------------------------------------------------------------------------------------
   ! Pitch angle PA = (-pi, pi) or fixed
   ! ------------------------------------------------------------------------------------------------------------------------

   IF (pitch_type == 1) THEN
      PA = pitch
   ELSEIF (pitch_type == 2) THEN
      CALL random_number(PA)
      PA = (2.0_wp*PA - 1.0_wp)
!      m1 = sum(PA)/Np
!      m2 = sum(PA**2)/Np       ! this was an old attempt to create a distribution with the correct first 2 moments; i found out that it induces correlation and leads to |PA|>1
!      PA = (PA - m1)*sqrt(1.0/3.0/(m2 - m1**2))
   ELSE
      WRITE (*, *) 'Error [initialize_particles]: unknown type for pitch angle distribution specification'
   END IF

   ! ------------------------------------------------------------------------------------------------------------------------
   ! Energy En distribution :: note, these are kinetic energies only.
   ! ------------------------------------------------------------------------------------------------------------------------
   Ts_eff = MAX(Ts, 1.0e-18_wp)   ! avoid pure zero's in energy without mutating the run parameter

   IF (energy_type == 1) THEN
      En = Es
   ELSEIF (energy_type == 2) THEN
      CALL PDF_Boltz(Np, Ts_eff, En)   ! the distribution of real kinetic energies mv^2/2
   ELSE
      WRITE (*, *) 'Error [initialize_particles]: unknown type for energy distribution specification'
   END IF

   Vp = SQRT(2.0_wp*En/As)*PA                      !SQRT(2*(En-Phi*Zw*phi0)/As)*COS(PA)

   ! ------------------------------------------------------------------------------------------------------------------------
   ! Initial positions
   ! ------------------------------------------------------------------------------------------------------------------------
   CALL initialize_particle_positions(X, Y, Z, Vp)

   ! ------------------------------------------------------------------------------------------------------------------------
   ! |B| and phi (used for mu and vp)
   ! ------------------------------------------------------------------------------------------------------------------------

   CALL magnetic_norm(X, Y, Z, B)                  ! the norm of the magnetic field evaluated at initial positions
   ! The former Onlycoord call populated q1/q2/q3 solely for potent_3D, which is
   ! currently disabled. q1 is computed from rhot below.
   ! CALL potent_3D(q1, q2, q3, t0, phi1)          ! turbulent field at initial positions

   mu = En/B*(1.0_wp - PA**2)                      !(En-Phi*Zw*phi0)/normB*SIN(PA)**2.
   ! CALL phi_zero(X,Y,Z,t0,phi2)                  ! the macroscopic electric potential field values evaluated at initial positions

   phi2 = Omgt0**2*Te/(Ti + Te)*(X - 1.0_wp)*(3.0_wp + X)/2.0_wp

   ! there is a question weather it is ok to evaluate the potential without gyro-average...
   ! ------------------------------------------------------------------------------------------------------------------------
   ! Note:: the energies are shifted with the potential value directly in "mu, vp" for simplicity

!  IF (USE_turb == ON) THEN
!     pb = E**(-Zs*(0.0*Phi*phi1 + phi2)/Ts- 0.0*0.5*(Zs*Phi/Ts)**2)
!  ELSE
!     pb = E**(-Zs*phi2/Ts)
!  END IF

!  pb = pb*X

!  pb = pb*Np/sum(pb)

!  Weight = 0.0_wp

   CALL equilibrium_rhot_all(X, Y, Z, rhot)
   q1 = C1*rhot
   delta_q1 = q1 - C1*rhot0

   temp = Ts_eff*exp(delta_q1*a0/C1*Lts)           !Ts*(1.0_wp + delta_q1*a0/C1*Lts)
   dens = 1.0_wp*exp(delta_q1*a0/C1*Lns)           !1.0*(1.0_wp + delta_q1*a0/C1*Lns)

! here we define 3 "weights":
! 1) G0 the marker distribution in the z phase space implemented :: it is a pure Maxwellian in energy (no jacobian resulting from there); the rhot*a0*X = r*R comes from the jacobian of the transformation between uniform (r,theta,varphi) coordinates and real space (x,y,z) - cartesian
! 2) F0 the initial distribution function (in practice its a maxwellian with exponential profiles)
! 3) beta the exponential weight associated with our modified deltaF scheme

   Jac = (rhot*a0*X)
   FM  = dens*exp(-En/temp)/temp**1.5_wp
   F0  = dens*exp(-En/temp)/temp**1.5_wp
   G0  = 1.0_wp*exp(-En/Ts_eff)/Ts_eff**(1.5_wp)/Jac
!  FM  = FM/(sum(FM*Jac)/Np)                       ! normalization of JF_M
!  F0  = F0/(sum(F0*Jac)/Np)                       ! normalization of JF_0
   G0  = G0/(sum(G0*Jac)/Np)                       ! normalization of JG_0

   G0 = G0*sum(F0/G0)/Np

   betaF = log(F0/FM)
   betaD = log(F0/FM)
   Einit = En
   alpha_1 = 0.0_wp
   alpha_2 = 0.0_wp
   alpha_3 = 0.0_wp

!  write(*,*) sum(G0*Jac)/Np,sum(F0*Jac)/Np,sum(FM*Jac)/Np,sum(F0/G0)/Np
!  pause
!  stop

END SUBROUTINE initialize_particles


!      ! Subroutine :: initialize_particle_positions
!      ! Purpose    :: fixes the initial GY coordinates of particles (only spatial)
!  ========================================================================================================================================================
!      ! Record of revisions:
!      !    Date       |   Programmer           |   Description of change
!      !  =======      |   =============        |   =========================
!      !  --/--/2021   |   D. I. Palade         |   Initial version
!      !  --/--/2023   |   L. M. Pomarjanschi   |   Corections and development
!      !  01/01/2024   |   D. I. Palade         |   Remake
!  ========================================================================================================================================================

SUBROUTINE initialize_particle_positions(X, Y, Z, Vp)
   USE constants
   IMPLICIT NONE

   ! I/O variables
   REAL(KIND=wp), DIMENSION(Np), INTENT(OUT) :: X, Y, Z
   REAL(KIND=wp), DIMENSION(Np)              :: Vp
   REAL(KIND=wp), DIMENSION(Np)              :: aux, aux2, aux1, aux3
   REAL(KIND=wp)                             :: rmin, rmax

   rmin = max(0.0_wp, 0.5_wp - 0.5_wp*annulus_width)
   rmax = min(1.0_wp, 0.5_wp + 0.5_wp*annulus_width)

   IF (position_type .eq. 1) THEN
      X = X0
      Y = Y0
      Z = Z0
      !         X = 1.
      !         Y = r00
   ELSEIF (position_type .eq. 2) THEN
      CALL random_number(aux)
      aux = pi*(2.0_wp*aux - 1.0_wp)

      IF (magnetic_model .eq. 1) THEN
         CALL psisurf_efit3(X, Y, Vp)
         Z = aux
      ELSEIF (magnetic_model .eq. 2) THEN
         CALL psisurf_solov2(X, Y)
         Z = aux
      ELSEIF (magnetic_model .eq. 3) THEN
         CALL psisurf_circ(X, Y, Vp)
         Z = aux
      END IF

   ELSEIF (position_type .eq. 3) THEN
      CALL random_number(aux)
      CALL random_number(aux2)
      CALL random_number(aux3)
      aux = pi*(2.0_wp*aux - 1.0_wp)
      aux3 = pi*(2.0_wp*aux3 - 1.0_wp)
      rmin = rmin*a0
      rmax = rmax*a0
!     aux2 = sqrt(aux2*(rmax**2 - rmin**2) + rmin**2)
      aux2 = rmin + aux2*(rmax - rmin)
      Y = aux2*sin(aux)
      X = 1.0_wp + aux2*cos(aux)
      Z = aux3
   ELSE
      WRITE (*, *) 'Error [initialize_particles]: unknown type for GC position distribution specification'
   END IF

END SUBROUTINE initialize_particle_positions


! a new version (2026) made because routine restructuring
SUBROUTINE potent_3D(q_1, q_2, q_3, time, phi0)
   USE constants
   USE random_numbers
   IMPLICIT NONE

   ! I/O variables
   REAL(KIND=wp), DIMENSION(Np), INTENT(IN)  :: q_1, q_2, q_3
   REAL(KIND=wp), INTENT(IN)                 :: time
   REAL(KIND=wp), DIMENSION(Np), INTENT(OUT) :: phi0

   REAL(KIND=wp), DIMENSION(Np) :: gpar, gprim, env

   ! Local variables
   INTEGER :: i

   ! Defaults (kept as in your original snippet)
   gpar  = exp((cos(q_3/C3) - 1.0_wp)/lbalonz**2)*balloon   ! this must be periodic
   gprim = -sin(q_3/C3)/lbalonz**2/C3*balloon               ! note that this is g'/g
   env   = noballoon + gpar                                 ! ballooning envelope multiplier

   IF (USE_real .eq. OFF) THEN
      !$OMP PARALLEL DO &
      !$OMP SHARED (q_1,q_2,q_3,time) &
      !$OMP PRIVATE (i)
      DO i = 1, Np
         phi0(i) = sum(L(:, i)*sin(kx(:, i)*q_1(i) + ky(:, i)*q_2(i) + kz(:, i)*q_3(i) - w(:, i)*time + ph(:, i))) ! sin(k_i x+ph_i) terms for all realizations (i) and all wavenumbers; summing over all waves in a realization
      END DO
      !$OMP END parallel DO
   ELSEIF (USE_real .eq. ON) THEN
      !$OMP PARALLEL DO &
      !$OMP SHARED (q_1,q_2,q_3,time) &
      !$OMP PRIVATE (i)
      DO i = 1, Np
         phi0(i) = sum(Ls(:, i)*sin(kxs(:)*q_1(i) + kys(:)*q_2(i) + kzs(:)*q_3(i) - ws(:)*time + phs(:)))
      END DO
      !$OMP END parallel DO
   END IF

   phi0 = phi0/sqrt(Nc/2.0)                            ! normalization

END SUBROUTINE potent_3D
