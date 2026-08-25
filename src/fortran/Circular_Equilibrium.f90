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
   REAL(KIND=wp), DIMENSION(Np), INTENT(IN)          :: X, Y
   REAL(KIND=wp), DIMENSION(Nqua, Np), INTENT(OUT)   :: R

   !---------------------------------------------------------------------------------
   ! Locals
   !---------------------------------------------------------------------------------
   REAL(KIND=wp), DIMENSION(Np) :: rr, theta, psiprim, psiprim2
   REAL(KIND=wp), DIMENSION(Np) :: rhot, rhotr, rhotz, rhotrr, rhotrz, rhotzz
   REAL(KIND=wp), DIMENSION(Np) :: psi, psir, psiz, psirr, psizz, psirz
   REAL(KIND=wp), DIMENSION(Np) :: chi, chir, chiz, qpsi, qprim, Fpsi, Fprim

   !---------------------------------------------------------------------------------
   ! Geometry / coordinates
   !---------------------------------------------------------------------------------
   rr    = sqrt((X - 1.0_wp)**2 + Y**2) + eps_rr
   theta = atan2(Y, X - 1.0_wp)

   chi  = 2.0_wp*atan( sqrt((1.0_wp - rr)/(1.0_wp + rr)) * tan(theta/2.0_wp) )
   chir =  sin(theta) * (rr**2 - X) / X / rr / sqrt(1.0_wp - rr**2)
   chiz = -(rr - cos(theta)) / rr / sqrt(1.0_wp - rr**2)

   rhot   = rr/a0
   rhotr  = (X - 1.0_wp)/rr/a0
   rhotz  = Y/rr/a0
   rhotrr = (1.0_wp/a0) * Y**2/rr**3
   rhotrz = (1.0_wp/a0) * Y*(1.0_wp - X)/rr**3
   rhotzz = (1.0_wp/a0) * (1.0_wp - X)**2/rr**3

   !---------------------------------------------------------------------------------
   ! q(r) and psi(r) + derivatives
   !---------------------------------------------------------------------------------
   qpsi  = s1 + s2*rhot + s3*rhot**2
   qprim = s2 + 2.0_wp*s3*rhot

   psiprim  = a0*rr/qpsi/sqrt(1.0_wp - rr**2)
   psiprim2 = -a0**2 * (rr/qpsi/sqrt(1.0_wp - rr**2)) * (qprim/a0/qpsi - 1.0_wp/(1.0_wp - rr**2)/rr)

   psir  = psiprim*rhotr
   psiz  = psiprim*rhotz
   psirr = psiprim*rhotrr + psiprim2*rhotr*rhotr
   psirz = psiprim*rhotrz + psiprim2*rhotr*rhotz
   psizz = psiprim*rhotzz + psiprim2*rhotz*rhotz

   Fpsi  = 1.0_wp
   Fprim = 0.0_wp

   ! Closed-form psi(r)
   psi = a0**2/s2*(rr/a0 - s1/s2*log(1.0_wp + s2/s1*rr/a0))

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
   REAL(KIND=wp), DIMENSION(Np), INTENT(OUT) :: X, Y
   REAL(KIND=wp), DIMENSION(Np), INTENT(IN)  :: Vp

   !---------------------------------------------------------------------------------
   ! Locals
   !---------------------------------------------------------------------------------
   REAL(KIND=wp), DIMENSION(Np)     :: aux
   REAL,          DIMENSION(2000)   :: aux1, aux2, aux3, psi, F, qpsi, psiprim, psir, psiz, normB
   INTEGER                          :: k, ioc

   CALL random_number(aux)
   CALL random_number(aux3)

   aux3 = a0*aux3
   aux  = pi*(2.0_wp*aux - 1.0_wp)     ! theta in [-pi, pi]

   F      = 1.0
   qpsi   = s1 + s2*aux3*(1.0/a0) + s3*aux3**2*(1.0/a0)**2
   psiprim = aux3/qpsi/sqrt(1.0 - aux3**2)

   DO k = 1, Np
      aux2 = aux3*sin(aux(k))
      aux1 = aux3*cos(aux(k))

      psir  = psiprim*aux1/aux3
      psiz  = psiprim*aux2/aux3
      normB = sqrt(F**2 + psir**2 + psiz**2) / (1.0 + aux2)

      psi = a0**2/s2*(aux3/a0 - s1/s2*log(1.0_wp + s2/s1*aux3/a0))

      ioc = minloc(abs(psi/psi0 - real(USE_PC,wp)*rhoi/R0*As/Zs*F/normB*Vp(k)/psi0 - 1.0), 1)

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
   REAL(wp), INTENT(OUT) :: X(Np), Y(Np)
   REAL(wp), INTENT(IN)  :: Vp(Np)

   ! Tunable radial grid size (kept 2000 to match original)
   INTEGER, PARAMETER :: Nr = 2000

   ! Cached, r-only tables (computed once and reused)
   REAL(wp), SAVE :: r_g(Nr), inv_base_norm(Nr), psi_over_psi0(Nr)
!   LOGICAL, SAVE  :: inited = .false.

   ! Per-call locals
   REAL(wp) :: th(Np)
   REAL(wp) :: s, c
   REAL(wp) :: alpha
   REAL(wp) :: val, best
   INTEGER  :: k, j, ioc

   ! Initialize r-grid & tables once
 !  IF (.NOT. inited) THEN
      CALL init_tables()
  !    inited = .true.
  ! END IF
   
   ! Random poloidal angles
   CALL random_number(th)
   th = pi*(2.0_wp*th - 1.0_wp)

   alpha = real(USE_PC,wp)*rhoi/R0*As/Zs/psi0   ! F=1 ⇒ K/psi0

   ! Main loop over particles
   DO k = 1, Np
      s = sin(th(k))
      c = cos(th(k))
      ! Minimize | psi/psi0 - alpha*Vp*(1 + s*r)/sqrt(1+psiprim^2) - 1 |
      best = huge(1.0_wp)
      ioc  = 1

      DO j = 1, Nr
         val = abs(psi_over_psi0(j) - alpha*Vp(k)*(1.0_wp + s*r_g(j))*inv_base_norm(j) - 1.0_wp)
         IF (val < best) THEN
            best = val
            ioc  = j
         END IF
      END DO

      ! Map (r, theta) -> (X, Y)
      X(k) = 1.0_wp + r_g(ioc)*c
      Y(k) = r_g(ioc)*s


      ! Safety clamp (should never trigger since r_g ∈ [0, a0])
      IF ((X(k) - 1.0_wp)*(X(k) - 1.0_wp) + Y(k)*Y(k) > a0*a0) THEN
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
      REAL(wp) :: r, rh, inv_a0, qpsi, psip, one_m_r2

      inv_a0 = 1.0_wp/a0

      ! Use mirpoints to avoid exactly r=0 and r=a0
      DO j = 1, Nr
         r      = (REAL(j, wp) - 0.5_wp) * a0 / REAL(Nr, wp)
         r_g(j) = r
      END DO

      DO j = 1, Nr
         r  = r_g(j)
         rh = r*inv_a0

         ! Horner: safe1 + safe2*(r/a0) + safe3*(r/a0)^2
         qpsi = (s3*rh + s2)*rh + s1

         one_m_r2 = max(1.0_wp - r*r, 1.0e-15_wp)
         psip     = r/(qpsi*sqrt(one_m_r2))                    ! psi'(r)

         inv_base_norm(j) = 1.0_wp / sqrt(1.0_wp + psip*psip)
         ! psi(r)/psi0
         psi_over_psi0(j) = &
            (a0**2/s2*(r/a0 - s1/s2*log(1.0_wp + s2/s1*r/a0))) / psi0
         
         


      END DO
   END SUBROUTINE init_tables

END SUBROUTINE psisurf_circ
