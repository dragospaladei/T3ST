!========================================================================================================
! Subroutine :: Circular
! Purpose    :: Circular equilibrium in which B = B_t + B_theta
!               B_t     = B0*R0/R
!               B_theta = B0*R0/R * r/q/sqrt(R0^2-r^2)
!
! Record of revisions:
!    Date       | Programmer    | Description
!  -----------  | ------------ | -----------------------------------------
!  01/01/2024   | D. I. Palade  | Remake
!  01/01/2025   | D. I. Palade  | Final check
!========================================================================================================

SUBROUTINE Circular(X, Y, R)
   USE constants
   IMPLICIT NONE

   !---------------------------------------------------------------------------------
   ! I/O
   !---------------------------------------------------------------------------------
   REAL(KIND=dp), DIMENSION(Np), INTENT(IN)          :: X, Y
   REAL(KIND=dp), DIMENSION(Nqua, Np), INTENT(OUT)   :: R

   !---------------------------------------------------------------------------------
   ! Locals
   !---------------------------------------------------------------------------------
   REAL(KIND=dp), DIMENSION(Np) :: rr, theta, psiprim, psiprim2
   REAL(KIND=dp), DIMENSION(Np) :: rhot, rhotr, rhotz, rhotrr, rhotrz, rhotzz
   REAL(KIND=dp), DIMENSION(Np) :: psi, psir, psiz, psirr, psizz, psirz
   REAL(KIND=dp), DIMENSION(Np) :: chi, chir, chiz, qpsi, qprim, Fpsi, Fprim
   INTEGER :: i

   !---------------------------------------------------------------------------------
   ! Geometry / coordinates
   !---------------------------------------------------------------------------------
   rr    = sqrt((X - 1.0_dp)**2 + Y**2) + eps_rr
   theta = atan2(Y, X - 1.0_dp)

   chi  = 2.0_dp*atan( sqrt((1.0_dp - rr)/(1.0_dp + rr)) * tan(theta/2.0_dp) )
   chir =  sin(theta) * (rr**2 - X) / X / rr / sqrt(1.0_dp - rr**2)
   chiz = -(rr - cos(theta)) / rr / sqrt(1.0_dp - rr**2)

   rhot   = rr/a0
   rhotr  = (X - 1.0_dp)/rr/a0
   rhotz  = Y/rr/a0
   rhotrr = (1.0_dp/a0) * Y**2/rr**3
   rhotrz = (1.0_dp/a0) * Y*(1.0_dp - X)/rr**3
   rhotzz = (1.0_dp/a0) * (1.0_dp - X)**2/rr**3

   !---------------------------------------------------------------------------------
   ! q(r) and psi(r) + derivatives
   !---------------------------------------------------------------------------------
   qpsi  = safe1 + safe2*rhot + safe3*rhot**2
   qprim = safe2 + 2.0_dp*safe3*rhot

   psiprim  = a0*rr/qpsi/sqrt(1.0_dp - rr**2)
   psiprim2 = -a0**2 * (rr/qpsi/sqrt(1.0_dp - rr**2)) * (qprim/a0/qpsi - 1.0_dp/(1.0_dp - rr**2)/rr)

   psir  = psiprim*rhotr
   psiz  = psiprim*rhotz
   psirr = psiprim*rhotrr + psiprim2*rhotr*rhotr
   psirz = psiprim*rhotrz + psiprim2*rhotr*rhotz
   psizz = psiprim*rhotzz + psiprim2*rhotz*rhotz

   Fpsi  = 1.0_dp
   Fprim = 0.0_dp

   ! Closed-form psi(r)
   psi = -(a0**2/safe3) / sqrt(1.0_dp + a0**2*safe1/safe3) * &
         ( atanh( sqrt(1.0_dp - rr**2) / sqrt(1.0_dp + a0**2*safe1/safe3) ) - &
           atanh( 1.0_dp / sqrt(1.0_dp + a0**2*safe1/safe3) ) )

   !---------------------------------------------------------------------------------
   ! Pack results
   !---------------------------------------------------------------------------------
   R( 1, :) = psi
   R( 2, :) = psir
   R( 3, :) = psiz
   R( 4, :) = psirr
   R( 5, :) = psirz
   R( 6, :) = psizz
   R( 7, :) = Fpsi
   R( 8, :) = Fprim
   R( 9, :) = qpsi
   R(10, :) = qprim
   R(11, :) = rhot
   R(12, :) = rhotr
   R(13, :) = rhotz
   R(14, :) = chi
   R(15, :) = chir
   R(16, :) = chiz

END SUBROUTINE Circular


!========================================================================================================
! psisurf_circ2
!========================================================================================================
SUBROUTINE psisurf_circ_old(X, Y, Vp)
   USE constants
   IMPLICIT NONE

   !---------------------------------------------------------------------------------
   ! I/O
   !---------------------------------------------------------------------------------
   REAL(KIND=dp), DIMENSION(Np), INTENT(OUT) :: X, Y
   REAL(KIND=dp), DIMENSION(Np), INTENT(IN)  :: Vp

   !---------------------------------------------------------------------------------
   ! Locals
   !---------------------------------------------------------------------------------
   REAL(KIND=dp), DIMENSION(Np)     :: aux
   REAL,          DIMENSION(2000)   :: aux1, aux2, aux3, psi, F, qpsi, psiprim, psir, psiz, normB
   INTEGER                          :: k, ioc

   CALL random_number(aux)
   CALL random_number(aux3)

   aux3 = a0*aux3
   aux  = pi*(2.0_dp*aux - 1.0_dp)     ! theta in [-pi, pi]

   F      = 1.0
   qpsi   = safe1 + safe2*aux3*(1.0/a0) + safe3*aux3**2*(1.0/a0)**2
   psiprim = aux3/qpsi/sqrt(1.0 - aux3**2)

   DO k = 1, Np
      aux2 = aux3*sin(aux(k))
      aux1 = aux3*cos(aux(k))

      psir  = psiprim*aux1/aux3
      psiz  = psiprim*aux2/aux3
      normB = sqrt(F**2 + psir**2 + psiz**2) / (1.0 + aux2)

      psi = -(a0**2/safe3) / sqrt(1.0 + a0**2*safe1/safe3) * &
            ( atanh( sqrt(1.0 - aux3**2) / sqrt(1.0 + a0**2*safe1/safe3) ) - &
              atanh( 1.0 / sqrt(1.0 + a0**2*safe1/safe3) ) )

      ioc = minloc(abs(psi/psi0 - real(USE_PC,dp)*rhoi/R0*Aw/Zw*F/normB*Vp(k)/psi0 - 1.0), 1)

      X(k) = 1.0 + aux1(ioc)
      Y(k) = aux2(ioc)

      IF (sqrt((X(k) - 1.0)**2 + Y(k)**2) > a0) THEN
         X(k) = X0
         Y(k) = Y0
      END IF
   END DO

END SUBROUTINE psisurf_circ_old


!========================================================================================================
! psisurf_circ
!========================================================================================================
SUBROUTINE psisurf_circ(X, Y, Vp)
   USE constants
   IMPLICIT NONE

   !---------------------------------------------------------------------------------
   ! I/O
   !---------------------------------------------------------------------------------
   REAL(dp), INTENT(OUT) :: X(Np), Y(Np)
   REAL(dp), INTENT(IN)  :: Vp(Np)

   ! Tunable radial grid size (kept 2000 to match original)
   INTEGER, PARAMETER :: Nr = 2000

   ! Cached, r-only tables (computed once and reused)
   REAL(dp), SAVE :: r_g(Nr), inv_base_norm(Nr), psi_over_psi0(Nr)
!   LOGICAL, SAVE  :: inited = .false.

   ! Per-call locals
   REAL(dp) :: th(Np)
   REAL(dp) :: s, c
   REAL(dp) :: alpha
   REAL(dp) :: val, best
   INTEGER  :: k, j, ioc

   ! Initialize r-grid & tables once
 !  IF (.NOT. inited) THEN
      CALL init_tables()
  !    inited = .true.
  ! END IF
   
   ! Random poloidal angles
   CALL random_number(th)
   th = pi*(2.0_dp*th - 1.0_dp)

   alpha = real(USE_PC,dp)*rhoi/R0*Aw/Zw/psi0   ! F=1 ⇒ K/psi0

   ! Main loop over particles
   DO k = 1, Np
      s = sin(th(k))
      c = cos(th(k))
      ! Minimize | psi/psi0 - alpha*Vp*(1 + s*r)/sqrt(1+psiprim^2) - 1 |
      best = huge(1.0_dp)
      ioc  = 1

      DO j = 1, Nr
         val = abs(psi_over_psi0(j) - alpha*Vp(k)*(1.0_dp + s*r_g(j))*inv_base_norm(j) - 1.0_dp)
         IF (val < best) THEN
            best = val
            ioc  = j
         END IF
      END DO

      ! Map (r, theta) -> (X, Y)
      X(k) = 1.0_dp + r_g(ioc)*c
      Y(k) = r_g(ioc)*s


      ! Safety clamp (should never trigger since r_g ∈ [0, a0])
      IF ((X(k) - 1.0_dp)*(X(k) - 1.0_dp) + Y(k)*Y(k) > a0*a0) THEN
         X(k) = X0
         Y(k) = Y0
      END IF
   END DO


CONTAINS

   SUBROUTINE init_tables()
      ! Build r-grid and r-only tables:
      !   r_g(j)             in [0, a0]
      !   inv_base_norm(j) = 1/sqrt(1 + psiprim(r)^2)
      !   psi_over_psi0(j) = psi(r)/psi0
      INTEGER  :: j
      REAL(dp) :: r, rh, inv_a0, qpsi, psip, one_m_r2
      REAL(dp) :: inv_sqrt_beta, cA, cB, beta

      inv_a0 = 1.0_dp/a0

      ! Use midpoints to avoid exactly r=0 and r=a0
      DO j = 1, Nr
         r      = (REAL(j, dp) - 0.5_dp) * a0 / REAL(Nr, dp)
         r_g(j) = r
      END DO

      beta          = 1.0_dp + (a0*a0)*safe1/safe3
      inv_sqrt_beta = 1.0_dp / sqrt(beta)

      cA = -(a0*a0)/safe3 * inv_sqrt_beta
      cB = atanh(inv_sqrt_beta)

      DO j = 1, Nr
         r  = r_g(j)
         rh = r*inv_a0

         ! Horner: safe1 + safe2*(r/a0) + safe3*(r/a0)^2
         qpsi = (safe3*rh + safe2)*rh + safe1

         one_m_r2 = max(1.0_dp - r*r, 1.0e-15_dp)
         psip     = r/(qpsi*sqrt(one_m_r2))                    ! psi'(r)

         inv_base_norm(j) = 1.0_dp / sqrt(1.0_dp + psip*psip)
         ! psi(r)/psi0
         psi_over_psi0(j) = cA*(atanh(sqrt(one_m_r2)*inv_sqrt_beta) - cB) / psi0
         
         


      END DO
   END SUBROUTINE init_tables

END SUBROUTINE psisurf_circ

