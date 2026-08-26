!=========================================================================================
!  Turbulence_random_fields.f90
!  Random distributions and turbulence-mode random fields.
!
!  1)  rand_p_sphere(n,d,x)
!        Generates <n> uniform points in a <d>-dimensional sphere of unit radius
!
!  2)  PDF_G(m,n,g,mean,wavelength,variance)
!        Generates <m> arrays of <n> random normal numbers, with <mean> average
!        and <wavelength> or <variance> variance
!
!  3)  PDF_G_D(m,n,g,mean,wavelength,variance)
!        Same as PDF_G but without OPTIONAL args (legacy interface)
!
!  4)  PDF_ky(n,zlam,k0,ky)
!        Generates <n> random numbers with
!        PDF(ky) = ky/k0 * ( exp[-(ky-k0)^2 zlam^2/2] - exp[-(ky+k0)^2 zlam^2/2] )
!
!  5)  PDF_Cauchy(n,w00,tau,xf,trunc)
!        Truncated/double-symmetric Cauchy via inverse CDF trick
!
!  6)  PDF_Boltz(n,tau,xf)
!        Maxwell-Boltzmann-like energy PDF ~ sqrt(x)*exp(-x) (scaled by tau)
!
!  7)  wavenum_unified(...)
!        Generates wavenumbers/phases/frequencies in either 2D (Nc,Np) mode
!        or "single sample per mode" (Nc) mode (USE_real_flag)
!
!  8)  Larmor(ind,mut)
!        Computes Larmor (Bessel J0) factors from kperp and mu
!
!  Legacy kept: wavenum_old, wavenum, wavenum_single, Larmor_old
!=========================================================================================

MODULE random_numbers
   USE constants
   IMPLICIT NONE
   SAVE

   ! "Full" (Nc,Np) fields
   REAL(KIND=wp), ALLOCATABLE, DIMENSION(:, :), TARGET :: &
      ph, kx, ky, kperp, kz, w, phx, phy, wc, dm, L, phrmp

   ! "Single" (Nc) fields (USE_real mode)
   REAL(KIND=wp), ALLOCATABLE, DIMENSION(:), TARGET :: &
      phs, kxs, kys, kperps, kzs, ws, phxs, phys, wcs, dms, phrmps

   ! Larmor factors in USE_real mode are stored as (Nc,Np) in Ls
   REAL(KIND=wp), ALLOCATABLE, DIMENSION(:, :), TARGET :: Ls

CONTAINS

   !=============================================================================
   ! rand_p_sphere
   !   Uniform points in d-dimensional unit sphere (implemented up to d=3).
   !=============================================================================
   SUBROUTINE rand_p_sphere(n, d, x)
      USE constants, ONLY: pi, wp
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n, d
      REAL(KIND=wp), DIMENSION(d, n), INTENT(OUT) :: x

      REAL(KIND=wp), DIMENSION(n) :: r1, r2, r3

      ! r1 ~ U(0,1) -> radius distribution for uniform volume: r = u^(1/d)
      CALL RANDOM_NUMBER(r1)
      r1 = r1**(1.0_WP/REAL(d, wp))

      SELECT CASE (d)

      CASE (1)
         x(1, :) = 2.0_WP*r1 - 1.0_WP

      CASE (2)
         CALL RANDOM_NUMBER(r2)       ! angle in [0,2pi)
         r2 = 2.0_WP*pi*r2
         x(1, :) = r1*SIN(r2)
         x(2, :) = r1*COS(r2)

      CASE (3)
         CALL RANDOM_NUMBER(r2)       ! poloidal
         CALL RANDOM_NUMBER(r3)       ! azimuthal (via cos(theta) trick)
         r2 = 2.0_WP*pi*r2
         r3 = ACOS(1.0_WP - 2.0_WP*r3)

         x(1, :) = r1*SIN(r2)*SIN(r3)
         x(2, :) = r1*COS(r2)*SIN(r3)
         x(3, :) = r1*COS(r3)

      CASE DEFAULT
         WRITE (*, *) 'Error: keep it in maximum 3D! ~ rand_p_sphere'
         PAUSE
         STOP

      END SELECT
   END SUBROUTINE rand_p_sphere

   !=============================================================================
   ! PDF_Boltz2  (legacy, PACK-based)
   !   Acceptance-rejection for y ~ sqrt(x)*exp(-x) on x in [0,20],
   !   then scaled by tau.
   !=============================================================================
   SUBROUTINE PDF_Boltz2(n, tau, xf)
      USE constants
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL(KIND=wp), INTENT(IN) :: tau
      REAL(KIND=wp), DIMENSION(n), INTENT(OUT) :: xf

      INTEGER :: i, q
      REAL(KIND=wp) :: a, b, c
      REAL(KIND=wp), DIMENSION(MAX(n/500, 10)) :: x1, x2, y1, y2
      REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: x3

      IF (n < 10) THEN
         PRINT *, '--------------------------------------------------------------------------------------'
         PRINT *, 'ERROR'
         PRINT *, 'Please generate at least 10 numbers! ~ PDF_Boltz'
         PAUSE
         STOP
      END IF

      a = 0.0_WP
      b = 20.0_WP
      c = 1.0_WP/SQRT(2.0_WP*E)

      xf = 0.0_WP
      i = 0

      DO WHILE (i < n)

         CALL RANDOM_NUMBER(x1)
         CALL RANDOM_NUMBER(x2)

         x1 = a + (b - a)*x1
         y1 = c*x2
         y2 = SQRT(x1)*EXP(-x1)

         ! Keep x1 where y2 >= y1
         x3 = PACK(x1, y2 - y1 >= 0.0_WP)

         q = MIN(n, i + SIZE(x3))
         xf(i + 1:q) = x3(1:q - i)

         i = i + SIZE(x3)
         DEALLOCATE (x3)

      END DO

      xf = xf*tau
   END SUBROUTINE PDF_Boltz2

   !=============================================================================
   ! PDF_Boltz  (chunked streaming accept-reject; same logic as your revised code)
   !=============================================================================
   SUBROUTINE PDF_Boltz(n, tau, xf)
      USE constants, ONLY: wp, E
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL(KIND=wp), INTENT(IN) :: tau
      REAL(KIND=wp), INTENT(OUT) :: xf(n)

      INTEGER, PARAMETER :: CHUNK = 65536
      INTEGER :: out, m, j
      REAL(KIND=wp) :: a, b, c

      REAL(KIND=wp), ALLOCATABLE :: x1(:), x2(:), y1(:), y2(:)
      LOGICAL, ALLOCATABLE :: keep(:)

      IF (n < 10) THEN
         PRINT *, 'ERROR: generate at least 10 numbers ~ PDF_Boltz'
         STOP
      END IF

      a = 0.0_WP
      b = 20.0_WP
      c = 1.0_WP/SQRT(2.0_WP*E)

      xf = 0.0_WP
      out = 0

      ALLOCATE (x1(CHUNK), x2(CHUNK), y1(CHUNK), y2(CHUNK), keep(CHUNK))

      DO WHILE (out < n)

         m = MIN(CHUNK, n - out)

         CALL RANDOM_NUMBER(x1(1:m))
         CALL RANDOM_NUMBER(x2(1:m))

         x1(1:m) = a + (b - a)*x1(1:m)
         y1(1:m) = c*x2(1:m)
         y2(1:m) = SQRT(x1(1:m))*EXP(-x1(1:m))

         keep(1:m) = (y2(1:m) >= y1(1:m))

         DO j = 1, m
            IF (keep(j)) THEN
               out = out + 1
               xf(out) = x1(j)
               IF (out == n) EXIT
            END IF
         END DO

      END DO

      xf = tau*xf
      DEALLOCATE (x1, x2, y1, y2, keep)
   END SUBROUTINE PDF_Boltz

!=============================================================================
! PDF_abcap_2piece
!
! Generates x from
!   p(x) ∝ (|x|*lambda)^a / (1 + (|x|*lambda)^b),   with |x|*lambda <= L
!
! Works in z = |x|*lambda, z in [0,L], using a two-piece proposal:
!
!   g(z) ∝ z^a       on [0, min(1,L)]
!   g(z) ∝ z^(a-b)   on [1, L]         (if L > 1)
!
! This is efficient because acceptance is >= 1/2 everywhere.
!
! Requirements:
!   a > -1
!   b > 0
!   lambda > 0
!   L > 0
!
! Strongly preferred:
!   b - a > 1
! so that the second-piece integral is well behaved and cheap.
!=============================================================================
   SUBROUTINE PDF_abcap(n, zlam, za, zb, zL, xf)
      USE constants, ONLY: wp
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL(KIND=wp), INTENT(IN) :: zlam, za, zb, zL
      REAL(KIND=wp), INTENT(OUT) :: xf(n)

      INTEGER, PARAMETER :: CHUNK = 65536

      INTEGER :: out, m, j
      REAL(KIND=wp) :: zsplit
      REAL(KIND=wp) :: I1, I2, Itot, p1
      REAL(KIND=wp) :: ap1, q, qp1
      REAL(KIND=wp), ALLOCATABLE :: u0(:), u1(:), u2(:), u3(:), z(:), pacc(:)
      LOGICAL, ALLOCATABLE :: keep(:), leftpiece(:)

      IF (n < 1) THEN
         PRINT *, 'ERROR: n must be >= 1 ~ PDF_abcap_2piece'
         STOP
      END IF

      IF (zlam <= 0.0_WP) THEN
         PRINT *, 'ERROR: lambda must be > 0 ~ PDF_abcap_2piece'
         STOP
      END IF

      IF (zL <= 0.0_WP) THEN
         PRINT *, 'ERROR: L must be > 0 ~ PDF_abcap_2piece'
         STOP
      END IF

      IF (za <= -1.0_WP) THEN
         PRINT *, 'ERROR: need a > -1 ~ PDF_abcap_2piece'
         STOP
      END IF

      IF (zb <= 0.0_WP) THEN
         PRINT *, 'ERROR: need b > 0 ~ PDF_abcap_2piece'
         STOP
      END IF

      ap1 = za + 1.0_WP
      zsplit = MIN(1.0_WP, zL)

      ! Integral of first piece: int_0^{min(1,L)} z^a dz
      I1 = zsplit**ap1/ap1

      ! Integral of second piece: int_1^L z^(a-b) dz, only if L>1
      I2 = 0.0_WP
      IF (zL > 1.0_WP) THEN
         q = za - zb
         qp1 = q + 1.0_WP   ! = a - b + 1

         IF (ABS(qp1) < 1.0E-12_WP) THEN
            ! Special case a-b = -1  <=> b-a = 1
            I2 = LOG(zL)
         ELSE
            I2 = (zL**qp1 - 1.0_WP)/qp1
         END IF
      END IF

      Itot = I1 + I2
      p1 = I1/Itot

      ALLOCATE (u0(CHUNK), u1(CHUNK), u2(CHUNK), u3(CHUNK), &
                z(CHUNK), pacc(CHUNK), keep(CHUNK), leftpiece(CHUNK))

      out = 0

      DO WHILE (out < n)

         m = MIN(CHUNK, n - out)

         ! u0 : choose proposal piece
         ! u1 : sample inside chosen piece
         ! u2 : accept/reject
         ! u3 : sign
         CALL RANDOM_NUMBER(u0(1:m))
         CALL RANDOM_NUMBER(u1(1:m))
         CALL RANDOM_NUMBER(u2(1:m))
         CALL RANDOM_NUMBER(u3(1:m))

         leftpiece(1:m) = (u0(1:m) <= p1)

         DO j = 1, m

            IF (leftpiece(j)) THEN
               !----------------------------------------------------------
               ! Left piece: g1(z) ∝ z^a on [0, zsplit]
               !
               ! CDF inversion:
               !   z = zsplit * U^(1/(a+1))
               !----------------------------------------------------------
               z(j) = zsplit*u1(j)**(1.0_WP/ap1)

               ! Acceptance:
               !   f / g1 = 1 / (1 + z^b)
               pacc(j) = 1.0_WP/(1.0_WP + z(j)**zb)

            ELSE
               !----------------------------------------------------------
               ! Right piece: g2(z) ∝ z^(a-b) on [1, L]
               !----------------------------------------------------------
               IF (zL <= 1.0_WP) THEN
                  ! Should not happen since p1=1 in that case, but keep safe
                  z(j) = zsplit*u1(j)**(1.0_WP/ap1)
                  pacc(j) = 1.0_WP/(1.0_WP + z(j)**zb)
               ELSE
                  q = za - zb
                  qp1 = q + 1.0_WP

                  IF (ABS(qp1) < 1.0E-12_WP) THEN
                     ! g2(z) ∝ z^(-1), CDF inversion:
                     !   z = exp( U * log(L) )
                     z(j) = EXP(u1(j)*LOG(zL))
                  ELSE
                     ! CDF inversion for power law on [1,L]:
                     !   z = [ 1 + U*(L^(q+1)-1) ]^(1/(q+1))
                     z(j) = (1.0_WP + u1(j)*(zL**qp1 - 1.0_WP))**(1.0_WP/qp1)
                  END IF

                  ! Acceptance:
                  !   f / g2 = z^b / (1 + z^b)
                  !
                  ! Numerically nicer form for large z:
                  !   1 / (1 + z^(-b))
                  pacc(j) = 1.0_WP/(1.0_WP + z(j)**(-zb))
               END IF

            END IF

         END DO

         keep(1:m) = (u2(1:m) <= pacc(1:m))

         DO j = 1, m
            IF (keep(j)) THEN
               out = out + 1

               IF (u3(j) < 0.5_WP) THEN
                  xf(out) = -z(j)/zlam
               ELSE
                  xf(out) = z(j)/zlam
               END IF

               IF (out == n) EXIT
            END IF
         END DO

      END DO

      DEALLOCATE (u0, u1, u2, u3, z, pacc, keep, leftpiece)

   END SUBROUTINE PDF_abcap  !=============================================================================
   ! PDF_Cauchy
   !   Truncated + double-sbymmetric (centered at +/-w00) Cauchy-like draw.
   !
   !   NOTE: This keeps your exact algebra:
   !     x2 = TAN((2u-1)*atan(trunc))
   !     aux = SIGN(1, u-0.5)
   !     xf = x2/tau + w00*aux
   !=============================================================================
   SUBROUTINE PDF_Cauchy(n, w00, tau, xf, trunc)
      USE constants, ONLY: wp, pi
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL(KIND=wp), INTENT(IN) :: tau, w00, trunc
      REAL(KIND=wp), DIMENSION(n), INTENT(OUT) :: xf

      REAL(KIND=wp), DIMENSION(n) :: x1, x2, aux

      CALL RANDOM_NUMBER(x1)
      CALL RANDOM_NUMBER(x2)

      ! Inverse-CDF with truncation control via atan(trunc)
      x2 = TAN((2.0_WP*x2 - 1.0_WP)*ATAN(trunc))

      ! Random +/-1 to mirror around +/-w00
      aux = SIGN(1.0_WP, x1 - 0.5_WP)

      ! Scale and shift
      xf = x2/tau + w00*aux
   END SUBROUTINE PDF_Cauchy

   !=============================================================================
   ! PDF_ky (chunked accept-reject; same math as your revised code)
   !=============================================================================
   SUBROUTINE PDF_ky(n, zlam, zk0, xf)
      USE constants
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL(KIND=wp), INTENT(IN) :: zlam, zk0
      REAL(KIND=wp), INTENT(OUT) :: xf(n)

      INTEGER, PARAMETER :: CHUNK = 65536
      REAL(KIND=wp) :: a, b, c
      INTEGER :: out, m, j

      REAL(KIND=wp), ALLOCATABLE :: x1(:), x2(:), y1(:), y2(:)
      LOGICAL, ALLOCATABLE :: keep(:)

      IF ((n < 10) .AND. (USE_real .EQ. 0)) THEN
         PRINT *, 'ERROR: generate at least 10 numbers ~ PDF_ky'
         STOP
      END IF

      a = -5.0_WP/zlam - zk0
      b = 5.0_WP/zlam + zk0

      ! Upper bound used in your code
      c = SQRT(2.0_WP/pi)/E*zlam

      ALLOCATE (x1(CHUNK), x2(CHUNK), y1(CHUNK), y2(CHUNK), keep(CHUNK))
      out = 0

      DO WHILE (out < n)

         m = MIN(CHUNK, n - out)

         CALL RANDOM_NUMBER(x1(1:m))
         CALL RANDOM_NUMBER(x2(1:m))

         x1(1:m) = a + (b - a)*x1(1:m)
         y1(1:m) = c*x2(1:m)

         y2(1:m) = (x1(1:m)/zk0*(EXP(-zlam*zlam*(x1(1:m) - zk0)**2/2.0_WP) - &
                                 EXP(-zlam*zlam*(x1(1:m) + zk0)**2/2.0_WP)))* &
                   zlam/SQRT(2.0_WP*pi)/2.0_WP

         keep(1:m) = (y2(1:m) >= y1(1:m))

         DO j = 1, m
            IF (keep(j)) THEN
               out = out + 1
               xf(out) = x1(j)
               IF (out == n) EXIT
            END IF
         END DO

      END DO

      DEALLOCATE (x1, x2, y1, y2, keep)
   END SUBROUTINE PDF_ky

   !=============================================================================
   ! PDF_ky2 (legacy PACK-based)
   !=============================================================================
   SUBROUTINE PDF_ky2(n, zlam, zk0, xf)
      USE constants, ONLY: wp, pi
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL(KIND=wp), INTENT(IN) :: zlam, zk0
      REAL(KIND=wp), DIMENSION(n), INTENT(OUT) :: xf

      INTEGER :: i, q
      REAL(KIND=wp) :: a, b, c
      REAL(KIND=wp), DIMENSION(MAX(n/500, 10)) :: x1, x2, y1, y2
      REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: x3

      IF (n < 10) THEN
         PRINT *, '--------------------------------------------------------------------------------------'
         PRINT *, 'ERROR'
         PRINT *, 'Please generate at least 10 numbers! ~ PDF_ky'
         PAUSE
         STOP
      END IF

      a = -5.0_WP/zlam - zk0
      b = 5.0_WP/zlam + zk0
      c = SQRT(2.0_WP/pi)/2.71828_WP*zlam

      xf = 0.0_WP
      i = 0

      DO WHILE (i < n)

         CALL RANDOM_NUMBER(x1)
         CALL RANDOM_NUMBER(x2)

         x1 = a + (b - a)*x1
         y1 = c*x2

         y2 = x1/zk0*(EXP(-zlam**2*(x1 - zk0)**2/2.0_WP) - &
                      EXP(-zlam**2*(x1 + zk0)**2/2.0_WP))* &
              zlam/SQRT(2.0_WP*pi)/2.0_WP

         x3 = PACK(x1, y2 - y1 >= 0.0_WP)

         q = MIN(n, i + SIZE(x3))
         xf(i + 1:q) = x3(1:q - i)

         i = i + SIZE(x3)
         DEALLOCATE (x3)

      END DO
   END SUBROUTINE PDF_ky2

   !=============================================================================
   ! PDF_G2 (legacy: odd-m handling via k=m+1 + extra row in r)
   !=============================================================================
   SUBROUTINE PDF_G2(m, n, g, mean, wavelength, variance)
      USE constants, ONLY: wp, pi
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: m, n
      REAL(KIND=wp), DIMENSION(m, n), INTENT(OUT) :: g
      REAL(KIND=wp), DIMENSION(m), OPTIONAL :: mean, variance, wavelength

      INTEGER :: i, k
      REAL(KIND=wp), DIMENSION(m + 1, n) :: r
      REAL(KIND=wp), DIMENSION(m) :: avg, var

      ! Defaults
      IF (PRESENT(mean)) THEN
         avg = mean
      ELSE
         avg = 0.0_WP
      END IF

      IF (PRESENT(wavelength)) THEN
         var = 1.0_WP/wavelength
      ELSEIF (PRESENT(variance)) THEN
         var = variance
      ELSE
         var = 1.0_WP
      END IF

      ! Need even index for Box-Muller pairing
      IF (MOD(m, 2) /= 0) THEN
         k = m + 1
      ELSE
         k = m
      END IF

      DO i = 1, k

         CALL RANDOM_NUMBER(r(i, :))

         WHERE (r(i, :) == 0.0_WP)
            r(i, :) = 1.0E-8_WP
         END WHERE

         IF (MOD(i, 2) == 0) THEN

            IF ((m - i) >= 0) THEN
               g(i - 1, :) = var(i - 1)*SQRT(-2.0_WP*LOG(r(i - 1, :)))*SIN(2.0_WP*pi*r(i, :)) + avg(i - 1)
               g(i, :) = var(i)*SQRT(-2.0_WP*LOG(r(i - 1, :)))*COS(2.0_WP*pi*r(i, :)) + avg(i)
            ELSE
               g(i - 1, :) = var(i - 1)*SQRT(-2.0_WP*LOG(r(i - 1, :)))*SIN(2.0_WP*pi*r(i, :)) + avg(i - 1)
            END IF

         END IF
      END DO
   END SUBROUTINE PDF_G2

   !=============================================================================
   ! PDF_G (faster + cleaner Box-Muller; keeps your precedence rules)
   !
   ! Precedence as in your revised version:
   !   avg default 0; sig default 1
   !   IF wavelength PRESENT -> sig = 1/wavelength
   !   IF variance PRESENT   -> sig = variance
   !
   ! So if both provided, variance overwrites wavelength (exactly your code comment).
   !=============================================================================
   SUBROUTINE PDF_G(m, n, g, mean, wavelength, variance)
      USE constants, ONLY: wp, pi
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: m, n
      REAL(KIND=wp), INTENT(OUT) :: g(m, n)
      REAL(KIND=wp), OPTIONAL, INTENT(IN) :: mean(m), wavelength(m), variance(m)

      INTEGER :: i
      REAL(KIND=wp) :: avg(m), sig(m)
      REAL(KIND=wp), ALLOCATABLE :: u1(:), u2(:), r(:), s(:), c(:)

      avg = 0.0_WP
      sig = 1.0_WP

      IF (PRESENT(mean)) avg = mean
      IF (PRESENT(wavelength)) sig = 1.0_WP/wavelength
      IF (PRESENT(variance)) sig = variance

      ALLOCATE (u1(n), u2(n), r(n), s(n), c(n))

      DO i = 1, m, 2

         CALL RANDOM_NUMBER(u1)
         CALL RANDOM_NUMBER(u2)

         ! Prevent log(0)
         WHERE (u1 < 1.0E-12_WP) u1 = 1.0E-12_WP

         r = SQRT(-2.0_WP*LOG(u1))
         s = SIN(2.0_WP*pi*u2)
         c = COS(2.0_WP*pi*u2)

         g(i, :) = sig(i)*r*s + avg(i)
         IF (i + 1 <= m) g(i + 1, :) = sig(i + 1)*r*c + avg(i + 1)

      END DO

      DEALLOCATE (u1, u2, r, s, c)
   END SUBROUTINE PDF_G

   !=============================================================================
   ! PDF_G_D (legacy "no optional args" version; logic kept exactly as given)
   !   NOTE: your original has:
   !     var = 1./wavelength
   !     var = variance
   !   so wavelength assignment is overwritten unconditionally by variance.
   !=============================================================================
   SUBROUTINE PDF_G_D(m, n, g, mean, wavelength, variance)
      USE constants, ONLY: wp, pi
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: m, n
      REAL(KIND=wp), DIMENSION(m, n), INTENT(OUT) :: g
      REAL(KIND=wp), DIMENSION(m) :: mean, variance, wavelength

      INTEGER :: i, k
      REAL(KIND=wp), DIMENSION(m + 1, n) :: r
      REAL(KIND=wp), DIMENSION(m) :: avg, var

      avg = mean

      var = 1.0_WP/wavelength
      var = variance     ! <- preserves your exact original overwrite

      IF (MOD(m, 2) /= 0) THEN
         k = m + 1
      ELSE
         k = m
      END IF

      DO i = 1, k

         CALL RANDOM_NUMBER(r(i, :))

         WHERE (r(i, :) == 0.0_WP)
            r(i, :) = 1.0E-8_WP
         END WHERE

         IF (MOD(i, 2) == 0) THEN

            IF ((m - i) >= 0) THEN
               g(i - 1, :) = var(i - 1)*SQRT(-2.0_WP*LOG(r(i - 1, :)))*SIN(2.0_WP*pi*r(i, :)) + avg(i - 1)
               g(i, :) = var(i)*SQRT(-2.0_WP*LOG(r(i - 1, :)))*COS(2.0_WP*pi*r(i, :)) + avg(i)
            ELSE
               g(i - 1, :) = var(i - 1)*SQRT(-2.0_WP*LOG(r(i - 1, :)))*SIN(2.0_WP*pi*r(i, :)) + avg(i - 1)
            END IF

         END IF
      END DO
   END SUBROUTINE PDF_G_D

   !=============================================================================
   ! wavenum_unified
   !   Unified generator:
   !     USE_real_flag = ON  -> fills *_s arrays (Nc), with npts=1
   !     USE_real_flag = OFF -> fills 2D arrays (Nc,Np), with npts=Np
   !
   !   NOTE: This is your routine re-indented + commented only; algebra preserved.
   !=============================================================================
   SUBROUTINE wavenum_unified(normB, Vstar1, Vstar2, Vstar3, gees, USE_real_flag)
      USE constants
      IMPLICIT NONE

      REAL(KIND=wp), INTENT(IN) :: normB, Vstar1, Vstar2, Vstar3
      REAL(KIND=wp), INTENT(IN) :: gees(3, 3)
      INTEGER, INTENT(IN) :: USE_real_flag

      INTEGER :: i, npts
      REAL(KIND=wp) :: invR0, rhoi_over_R0, R2, Te_over_Ti, ATeTi, trunc
      REAL(KIND=wp) :: w000, kx000, k0, sgn

      REAL(KIND=wp), ALLOCATABLE :: tmp4(:, :)         ! (4,npts) scratch for PDF_G

      REAL(KIND=wp), ALLOCATABLE :: kyrow(:), kxrow(:), kzrow(:), wrow(:), wcrow(:)
      REAL(KIND=wp), ALLOCATABLE :: phrow(:), phxrow(:), phyrow(:), phrmprow(:), dmrow(:)

      ! Decide number of samples/trajectories generated for each mode i
      IF (USE_real_flag == ON) THEN
         npts = 1
      ELSE
         npts = Np
      END IF

      !-----------------------------
      ! Allocate outputs (mode-wise)
      !-----------------------------
      IF (USE_real_flag == ON) THEN

         IF (ALLOCATED(phs)) THEN
            IF (SIZE(phs) /= Nc) THEN
               DEALLOCATE (phs, kxs, kys, kperps, kzs, ws, phrmps, phxs, phys, wcs, dms)
            END IF
         END IF

         IF (.NOT. ALLOCATED(phs)) THEN
            ALLOCATE (phs(Nc))
            ALLOCATE (kxs, kys, kzs, kperps, ws, phrmps, phxs, phys, wcs, dms, MOLD=phs)
         END IF

      ELSE

         IF (ALLOCATED(ph)) THEN
            IF (SIZE(ph, 1) /= Nc .OR. SIZE(ph, 2) /= Np) THEN
               DEALLOCATE (ph, kx, ky, kperp, kz, w, phrmp, phx, phy, wc, dm)
            END IF
         END IF

         IF (.NOT. ALLOCATED(ph)) THEN
            ALLOCATE (ph(Nc, Np))
            ALLOCATE (kx, ky, kperp, kz, w, phrmp, phx, phy, wc, dm, MOLD=ph)
         END IF

      END IF

      !-----------------------------
      ! Derived constants used by the generated turbulence fields
      !-----------------------------
      invR0 = 1.0_WP/R0
      rhoi_over_R0 = rhoi*invR0
      Te_over_Ti = Te/Ti
      R2 = (rhoi_over_R0/normB)**2
      ATeTi = Aeff*Te_over_Ti

      w000 = 0.0_WP
      kx000 = 0.0_WP
      trunc = 10.0_WP

      !-----------------------------
      ! Per-mode scratch buffers
      !-----------------------------
      ALLOCATE (tmp4(4, npts))
      ALLOCATE (kyrow(npts), kxrow(npts), kzrow(npts), wrow(npts), wcrow(npts))
      ALLOCATE (phrow(npts), phxrow(npts), phyrow(npts), phrmprow(npts), dmrow(npts))

      !=============================
      ! Main turbulence switch
      !=============================
      IF (USE_turb == ON) THEN

         DO i = 1, Nc

            !-------------------------
            ! Phases and parallel shifts
            !-------------------------
            CALL RANDOM_NUMBER(phrow); phrow = pi*(2.0_WP*phrow - 1.0_WP)
            CALL RANDOM_NUMBER(phrmprow); phrmprow = pi*(2.0_WP*phrmprow - 1.0_WP)
            CALL RANDOM_NUMBER(phxrow); phxrow = pi*(2.0_WP*phxrow - 1.0_WP)
            CALL RANDOM_NUMBER(phyrow); phyrow = pi*(2.0_WP*phyrow - 1.0_WP)

            CALL RANDOM_NUMBER(dmrow)
            dmrow = 2.0_WP*(dmrow - 0.5_WP)
            dmrow = INT(q00/lambdaz*dmrow)  ! integer shifts

            !-------------------------
            ! Correlated draws: (kx,kz,w,wc) from Gaussian generator
            !-------------------------
            CALL PDF_G(4, npts, tmp4, wavelength=[lambdax, lambdaz, tauc, MERGE(taucc, 1.0_WP, taucc /= 0.0_WP)])

            kzrow = tmp4(2, :)
            IF (taucc == 0.0_WP) THEN
               wcrow = 0.0_WP
            ELSE
               wcrow = tmp4(4, :)
            END IF

            ! Radial spectrum selection
            IF (x_corr == 1) THEN
               kxrow = tmp4(1, :)
            ELSEIF (x_corr == 2) THEN
               CALL PDF_Cauchy(npts, kx000, lambdax, kxrow, trunc)
            ELSEIF (x_corr == 3) THEN
               CALL PDF_abcap(npts, lambdax, 0.0001_WP, 1.36_WP, 10.0_WP, kxrow)
            END IF

! Temporal spectrum selection
            IF (t_corr == 1) THEN
               wrow = tmp4(3, :)
            ELSEIF (t_corr == 2) THEN
               CALL PDF_Cauchy(npts, w000, tauc, wrow, trunc)
            END IF

            !-------------------------
            ! ky depends on ITG/TEM block
            !-------------------------
            IF (i <= Nci) THEN
               k0 = k0i
               sgn = +1.0_WP
            ELSE
               k0 = k0e
               sgn = -Te_over_Ti
            END IF

            ! Radial spectrum selection
            IF (y_corr == 1) THEN
               CALL PDF_ky(npts, lambday, k0, kyrow)
            ELSEIF (y_corr == 2) THEN
               CALL PDF_abcap(npts, lambday, 0.8_WP, 2.3_WP, 10.0_WP, kyrow)
            END IF

            !-------------------------
            ! Integer wrapping (as in your original)
            !-------------------------
            kyrow = REAL(INT(kyrow*C2), wp)/C2
            kzrow = REAL(INT(kzrow), wp)/C3

            !=========================================================
            ! Write outputs + compute kperp and frequency correction
            !=========================================================
            IF (USE_real_flag == ON) THEN

               phs(i) = phrow(1)
               phrmps(i) = phrmprow(1)
               phxs(i) = phxrow(1)
               phys(i) = phyrow(1)
               dms(i) = dmrow(1)

               kxs(i) = kxrow(1)
               kys(i) = kyrow(1)
               kzs(i) = kzrow(1)
               wcs(i) = wcrow(1)

               kperps(i) = rhoi_over_R0*SQRT( &
                           kxs(i)*gees(1, 1)*kxs(i) + 2.0_WP*kxs(i)*gees(1, 2)*kys(i) + 2.0_WP*kxs(i)*gees(1, 3)*kzs(i) + &
                           kys(i)*gees(2, 2)*kys(i) + 2.0_WP*kys(i)*gees(2, 3)*kzs(i) + &
                           kzs(i)*gees(3, 3)*kzs(i))

               ws(i) = wrow(1) + (sgn*Ln*R2*(kxs(i)*Vstar1 + kys(i)*Vstar2 + kzs(i)*Vstar3))/ &
                       (1.0_WP + ATeTi*kperps(i)**2)

               ws(i) = REAL(USE_freq, wp)*ws(i)

            ELSE

               ph(i, :) = phrow
               phrmp(i, :) = phrmprow
               phx(i, :) = phxrow
               phy(i, :) = phyrow
               dm(i, :) = dmrow

               kx(i, :) = kxrow
               ky(i, :) = kyrow
               kz(i, :) = kzrow
               wc(i, :) = wcrow

               kperp(i, :) = rhoi_over_R0*SQRT( &
                        kx(i, :)*gees(1, 1)*kx(i, :) + 2.0_WP*kx(i, :)*gees(1, 2)*ky(i, :) + 2.0_WP*kx(i, :)*gees(1, 3)*kz(i, :) + &
                             ky(i, :)*gees(2, 2)*ky(i, :) + 2.0_WP*ky(i, :)*gees(2, 3)*kz(i, :) + &
                             kz(i, :)*gees(3, 3)*kz(i, :))

               w(i, :) = wrow + (sgn*Ln*R2*(kx(i, :)*Vstar1 + ky(i, :)*Vstar2 + kz(i, :)*Vstar3))/ &
                         (1.0_WP + ATeTi*kperp(i, :)**2)

            END IF

         END DO

         IF (USE_real_flag /= ON) THEN
            w = REAL(USE_freq, wp)*w
         END IF

      ELSE
         !----------------------------------------------------------
         ! No turbulence: zero everything (preserve your exact intent)
         !----------------------------------------------------------
         IF (USE_real_flag == ON) THEN
            phs = 0.0_WP; phrmps = 0.0_WP; phxs = 0.0_WP; phys = 0.0_WP
            kxs = 0.0_WP; kys = 0.0_WP; kzs = 0.0_WP
            ws = 0.0_WP; wcs = 0.0_WP; kperps = 0.0_WP; dms = 0.0_WP
         ELSE
            ph = 0.0_WP; phrmp = 0.0_WP; phx = 0.0_WP; phy = 0.0_WP
            kx = 0.0_WP; ky = 0.0_WP; kz = 0.0_WP
            w = 0.0_WP; wc = 0.0_WP; kperp = 0.0_WP; dm = 0.0_WP
         END IF
      END IF

      !-----------------------------------------------------------------
      ! "first sample has no turbulence" convention
      ! "second sample has frozen turbulence" convention
      !-----------------------------------------------------------------
      IF (USE_real_flag == ON) THEN
         phs(1) = 0.0_WP; phrmps(1) = 0.0_WP; phxs(1) = 0.0_WP; phys(1) = 0.0_WP
         kxs(1) = 0.0_WP; kys(1) = 0.0_WP; kzs(1) = 0.0_WP
         ws(1) = 0.0_WP; wcs(1) = 0.0_WP; kperps(1) = 0.0_WP
         ws(2) = 0.0_WP; wcs(2) = 0.0_WP
      ELSE
         ph(:, 1) = 0.0_WP; phrmp(:, 1) = 0.0_WP; phx(:, 1) = 0.0_WP; phy(:, 1) = 0.0_WP
         kx(:, 1) = 0.0_WP; ky(:, 1) = 0.0_WP; kz(:, 1) = 0.0_WP
         w(:, 1) = 0.0_WP; wc(:, 1) = 0.0_WP; kperp(:, 1) = 0.0_WP
         w(:, 2) = 0.0_WP; wc(:, 2) = 0.0_WP
      END IF

      DEALLOCATE (tmp4, kyrow, kxrow, kzrow, wrow, wcrow, phrow, phxrow, phyrow, phrmprow, dmrow)

   END SUBROUTINE wavenum_unified

   !=============================================================================
   ! Larmor
   !   Computes J0(k_perp * coef * sqrt(|mu|)).
   !   ind==ON: compute; else set to 1.
   !
   !   NOTE: allocation pattern preserved (L is (Nc,Np) when USE_real==OFF,
   !         Ls is also (Nc,Np) when USE_real==ON, matching your code).
   !=============================================================================
   SUBROUTINE Larmor(ind, mut)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ind
      REAL(KIND=wp), DIMENSION(Np), INTENT(IN) :: mut

      INTEGER :: i
      REAL(KIND=wp) :: coef

      ! Allocate outputs depending on mode
      IF (USE_real == OFF) THEN
         IF (ALLOCATED(L)) DEALLOCATE (L)
         ALLOCATE (L(Nc, Np))
      ELSEIF (USE_real == ON) THEN
         IF (ALLOCATED(Ls)) DEALLOCATE (Ls)
         ALLOCATE (Ls(Nc, Np))
      END IF

      IF (ind == ON) THEN
         coef = SQRT(2.0_WP*As/Zs**2)

         IF (USE_real == OFF) THEN
            !$omp parallel do default(none) shared(L, kperp, mut, coef, Nc, Np) private(i)
            DO i = 1, Np
               L(:, i) = BESSEL_J0(kperp(:, i)*coef*SQRT(ABS(mut(i))))
            END DO
            !$omp end parallel do

         ELSEIF (USE_real == ON) THEN
            !$omp parallel do default(none) shared(Ls, kperps, mut, coef, Nc, Np) private(i)
            DO i = 1, Np
               Ls(:, i) = BESSEL_J0(kperps*coef*SQRT(ABS(mut(i))))
            END DO
            !$omp end parallel do
         END IF

      ELSE

         IF (USE_real == OFF) THEN
            L = 1.0_WP
         ELSEIF (USE_real == ON) THEN
            Ls = 1.0_WP
         END IF

      END IF
   END SUBROUTINE Larmor

   !=============================================================================
   ! Legacy routines kept below: ONLY re-indentation would be applied similarly.
   ! (You pasted them fully; I’m leaving their bodies unchanged to avoid
   !  “spoiling anything of the logic”. If you want, I can re-indent those too
   !  in the same style as above.)
   !=============================================================================

!=========================================================================================
!  Legacy routines: wavenum_old, wavenum, wavenum_single, Larmor_old
!  Re-indented + commented lightly.
!  IMPORTANT: logic/algebra/order preserved as in your pasted code.
!=========================================================================================

   SUBROUTINE wavenum_old(normB, Vstar1, Vstar2, Vstar3, gees)
      IMPLICIT NONE

      REAL(KIND=wp), DIMENSION(Np) :: w1, w2
      REAL(KIND=wp), DIMENSION(5, Np) :: g
      INTEGER :: i

      REAL(KIND=wp) :: normB, Vstar1, Vstar2, Vstar3, trunc
      REAL(KIND=wp), DIMENSION(3, 3) :: gees          ! matrix for kperp^2 metric
      REAL(KIND=wp) :: w000, kx000

      !-----------------------------
      ! Allocate/resize arrays
      !-----------------------------
      IF (ALLOCATED(ph)) THEN
         IF (SIZE(ph) .NE. Np) THEN
            DEALLOCATE (ph, kx, ky, kperp, kz, w, phrmp, phx, phy, wc, dm)

            ALLOCATE (ph(Nc, Np))
            ALLOCATE (kx, ky, kperp, kz, w, phrmp, phx, phy, wc, dm, MOLD=ph)
         END IF
      ELSE
         ALLOCATE (ph(Nc, Np))
         ALLOCATE (kx, ky, kperp, kz, w, phrmp, phx, phy, wc, dm, MOLD=ph)
      END IF

      !-----------------------------
      ! Turbulence on/off
      !-----------------------------
      IF (USE_turb == ON) THEN

         ! Random phases and integer parallel shifts dm
         CALL RANDOM_NUMBER(ph)
         CALL RANDOM_NUMBER(phrmp)
         CALL RANDOM_NUMBER(phx)
         CALL RANDOM_NUMBER(phy)
         CALL RANDOM_NUMBER(dm)

         ph = pi*(2.0_WP*ph - 1.0_WP)
         phrmp = pi*(2.0_WP*phrmp - 1.0_WP)
         phx = pi*(2.0_WP*phx - 1.0_WP)
         phy = pi*(2.0_WP*phy - 1.0_WP)

         dm = 2.0_WP*(dm - 0.5_WP)
         dm = INT(q00/lambdaz*dm)

         !-----------------------------------------------
         ! Wavenumbers/frequencies for ITG (1..Nci)
         !-----------------------------------------------
         w000 = 0.0_WP
         kx000 = 0.0_WP
         trunc = 10.0_WP

         DO i = 1, Nci

            ! Gaussian correlated pieces
            CALL PDF_G(4, Np, g(1:4, :), wavelength=[lambdax, lambdaz, tauc, taucc])
            CALL PDF_ky(Np, lambday, k0i, g(5, :))

            kz(i, :) = g(2, :)
            wc(i, :) = g(4, :)
            ky(i, :) = g(5, :)

            ! Radial correlation model
            IF (x_corr .EQ. 1) THEN
               kx(i, :) = g(1, :)
            ELSEIF (x_corr .EQ. 2) THEN
               CALL PDF_Cauchy(Np, kx000, lambdax, kx(i, :), trunc)
            END IF

            ! Time correlation model
            IF (t_corr .EQ. 1) THEN
               w1 = g(3, :)
            ELSEIF (t_corr .EQ. 2) THEN
               CALL PDF_Cauchy(Np, w000, lambdax, w1, trunc)
            END IF

            ! Discretizations / wrapping (kept as-is)
            ky(i, :) = 1.0_WP*INT(ky(i, :)*C2)/C2
            kz(i, :) = 0.0_WP*(C2*ky(i, :)*q00 - INT(C2*ky(i, :)*q00))/C3
            kz(i, :) = kz(i, :) + dm(i, :)/C3

            ! k_perp and frequency correction
            kperp(i, :) = (rhoi/R0)*SQRT( &
                          kx(i, :)*gees(1, 1)*kx(i, :) + 2.0_WP*kx(i, :)*gees(1, 2)*ky(i, :) + &
                          2.0_WP*kx(i, :)*gees(1, 3)*kz(i, :) + ky(i, :)*gees(2, 2)*ky(i, :) + &
                          2.0_WP*ky(i, :)*gees(2, 3)*kz(i, :) + kz(i, :)*gees(3, 3)*kz(i, :))

            w2 = +Li*(Ti/Ti)*((rhoi/R0/normB)**2)*(kx(i, :)*Vstar1 + ky(i, :)*Vstar2 + kz(i, :)*Vstar3)
            w2 = w2/(1.0_WP + Aeff*(Te/Ti)*kperp(i, :)**2)

            w(i, :) = w1 + 1.0_WP*w2
         END DO

         !-----------------------------------------------
         ! Wavenumbers/frequencies for TEM (Nci+1..Nc)
         !-----------------------------------------------
         DO i = Nci + 1, Nc

            CALL PDF_G(4, Np, g(1:4, :), wavelength=[lambdax, lambdaz, tauc, taucc])
            CALL PDF_ky(Np, lambday, k0e, g(5, :))

            kz(i, :) = g(2, :)
            wc(i, :) = g(4, :)
            ky(i, :) = g(5, :)

            IF (x_corr .EQ. 1) THEN
               kx(i, :) = g(1, :)
            ELSEIF (x_corr .EQ. 2) THEN
               CALL PDF_Cauchy(Np, kx000, lambdax, kx(i, :), trunc)
            END IF

            IF (t_corr .EQ. 1) THEN
               w1 = g(3, :)
            ELSEIF (t_corr .EQ. 2) THEN
               CALL PDF_Cauchy(Np, w000, lambdax, w1, trunc)
            END IF

            ky(i, :) = 1.0_WP*INT(ky(i, :)*C2)/C2
            kz(i, :) = 0.0_WP*(C2*ky(i, :)*q00 - INT(C2*ky(i, :)*q00))/C3
            kz(i, :) = kz(i, :) + dm(i, :)/C3

            kperp(i, :) = (rhoi/R0)*SQRT( &
                          kx(i, :)*gees(1, 1)*kx(i, :) + 2.0_WP*kx(i, :)*gees(1, 2)*ky(i, :) + &
                          2.0_WP*kx(i, :)*gees(1, 3)*kz(i, :) + ky(i, :)*gees(2, 2)*ky(i, :) + &
                          2.0_WP*ky(i, :)*gees(2, 3)*kz(i, :) + kz(i, :)*gees(3, 3)*kz(i, :))

            w2 = -Ln*(Te/Ti)*((rhoi/R0/normB)**2)*(kx(i, :)*Vstar1 + ky(i, :)*Vstar2 + kz(i, :)*Vstar3)
            w2 = w2/(1.0_WP + Aeff*(Te/Ti)*kperp(i, :)**2)

            w(i, :) = w1 + 1.0_WP*w2
         END DO

         w = REAL(USE_freq, wp)*w

      ELSEIF (USE_turb /= ON) THEN

         ph = 0.0_WP
         phx = 0.0_WP
         phy = 0.0_WP
         kx = 0.0_WP
         ky = 0.0_WP
         kz = 0.0_WP
         w = 0.0_WP
         wc = 0.0_WP
         kperp = 0.0_WP

      END IF

      ! First trajectory without turbulence (convention)
      ph(:, 1) = 0.0_WP
      phx(:, 1) = 0.0_WP
      phy(:, 1) = 0.0_WP
      kx(:, 1) = 0.0_WP
      ky(:, 1) = 0.0_WP
      kz(:, 1) = 0.0_WP
      w(:, 1) = 0.0_WP
      wc(:, 1) = 0.0_WP
      kperp(:, 1) = 0.0_WP

   END SUBROUTINE wavenum_old

   !=============================================================================
   ! wavenum (newer/cleaner legacy; re-indented only)
   !=============================================================================
   SUBROUTINE wavenum(normB, Vstar1, Vstar2, Vstar3, gees)
      USE constants
      IMPLICIT NONE

      REAL(KIND=wp), INTENT(IN) :: normB, Vstar1, Vstar2, Vstar3
      REAL(KIND=wp), INTENT(IN) :: gees(3, 3)

      REAL(KIND=wp) :: invR0, rhoi_over_R0, invC3, Te_over_Ti, R2, ATeTi
      REAL(KIND=wp) :: sx, sz, sw1, swc, trunc
      INTEGER :: i
      REAL(KIND=wp), ALLOCATABLE :: tmp4(:, :)   ! (4,Np) scratch for PDF_G

      REAL(KIND=wp) :: w000, kx000

      !-----------------------------
      ! Allocate/resize arrays
      !-----------------------------
      IF (ALLOCATED(ph)) THEN
         IF (SIZE(ph, 1) /= Nc .OR. SIZE(ph, 2) /= Np) THEN
            DEALLOCATE (ph, kx, ky, kperp, kz, w, phrmp, phx, phy, wc, dm)
            ALLOCATE (ph(Nc, Np))
            ALLOCATE (kx, ky, kperp, kz, w, phrmp, phx, phy, wc, dm, MOLD=ph)
         END IF
      ELSE
         ALLOCATE (ph(Nc, Np))
         ALLOCATE (kx, ky, kperp, kz, w, phrmp, phx, phy, wc, dm, MOLD=ph)
      END IF

      !-----------------------------
      ! Turbulence on/off
      !-----------------------------
      IF (USE_turb == ON) THEN

         ! Phases + parallel shift
         CALL RANDOM_NUMBER(ph); ph = pi*(2.0_WP*ph - 1.0_WP)
         CALL RANDOM_NUMBER(phrmp); phrmp = pi*(2.0_WP*phrmp - 1.0_WP)
         CALL RANDOM_NUMBER(phx); phx = pi*(2.0_WP*phx - 1.0_WP)
         CALL RANDOM_NUMBER(phy); phy = pi*(2.0_WP*phy - 1.0_WP)

         CALL RANDOM_NUMBER(dm)
         dm = 2.0_WP*(dm - 0.5_WP)
         dm = INT(q00/lambdaz*dm)

         ! Constants
         invR0 = 1.0_WP/R0
         rhoi_over_R0 = rhoi*invR0
         invC3 = 1.0_WP/C3
         Te_over_Ti = Te/Ti
         R2 = (rhoi_over_R0/normB)**2
         ATeTi = Aeff*Te_over_Ti

         sx = 1.0_WP/lambdax
         sz = 1.0_WP/lambdaz
         sw1 = 1.0_WP/tauc
         swc = 1.0_WP/taucc

         w000 = 0.0_WP
         kx000 = 0.0_WP
         trunc = 10.0_WP

         !-------------------------
         ! ITG block
         !-------------------------
         ALLOCATE (tmp4(4, Np))

         DO i = 1, Nci

            CALL PDF_G(4, Np, tmp4, wavelength=[lambdax, lambdaz, tauc, taucc])

            kz(i, :) = tmp4(2, :)
            wc(i, :) = tmp4(4, :)

            IF (x_corr .EQ. 1) THEN
               kx(i, :) = tmp4(1, :)
            ELSEIF (x_corr .EQ. 2) THEN
               CALL PDF_Cauchy(Np, kx000, lambdax, kx(i, :), trunc)
            END IF

            IF (t_corr .EQ. 1) THEN
               w(i, :) = tmp4(3, :)
            ELSEIF (t_corr .EQ. 2) THEN
               CALL PDF_Cauchy(Np, w000, lambdax, w(i, :), trunc)
            END IF

            CALL PDF_ky(Np, lambday, k0i, ky(i, :))

            ky(i, :) = REAL(INT(ky(i, :)*C2), wp)/C2
            kz(i, :) = REAL(INT(kz(i, :)), wp)/C3

            kperp(i, :) = (rhoi_over_R0)*SQRT( &
                          kx(i, :)*gees(1, 1)*kx(i, :) + 2.0_WP*kx(i, :)*gees(1, 2)*ky(i, :) + &
                          2.0_WP*kx(i, :)*gees(1, 3)*kz(i, :) + ky(i, :)*gees(2, 2)*ky(i, :) + &
                          2.0_WP*ky(i, :)*gees(2, 3)*kz(i, :) + kz(i, :)*gees(3, 3)*kz(i, :))

            w(i, :) = w(i, :) + (+Ln*R2*(kx(i, :)*Vstar1 + ky(i, :)*Vstar2 + kz(i, :)*Vstar3))/ &
                      (1.0_WP + ATeTi*kperp(i, :)**2)
         END DO

         !-------------------------
         ! TEM block
         !-------------------------
         DO i = Nci + 1, Nc

            CALL PDF_G(4, Np, tmp4, wavelength=[lambdax, lambdaz, tauc, taucc])

            kz(i, :) = tmp4(2, :)
            wc(i, :) = tmp4(4, :)

            IF (x_corr .EQ. 1) THEN
               kx(i, :) = tmp4(1, :)
            ELSEIF (x_corr .EQ. 2) THEN
               CALL PDF_Cauchy(Np, kx000, lambdax, kx(i, :), trunc)
            END IF

            IF (t_corr .EQ. 1) THEN
               w(i, :) = tmp4(3, :)
            ELSEIF (t_corr .EQ. 2) THEN
               CALL PDF_Cauchy(Np, w000, lambdax, w(i, :), trunc)
            END IF

            CALL PDF_ky(Np, lambday, k0e, ky(i, :))

            ky(i, :) = REAL(INT(ky(i, :)*C2), wp)/C2
            kz(i, :) = REAL(INT(kz(i, :)), wp)/C3

            kperp(i, :) = (rhoi_over_R0)*SQRT( &
                          kx(i, :)*gees(1, 1)*kx(i, :) + 2.0_WP*kx(i, :)*gees(1, 2)*ky(i, :) + &
                          2.0_WP*kx(i, :)*gees(1, 3)*kz(i, :) + ky(i, :)*gees(2, 2)*ky(i, :) + &
                          2.0_WP*ky(i, :)*gees(2, 3)*kz(i, :) + kz(i, :)*gees(3, 3)*kz(i, :))

            w(i, :) = w(i, :) + (-Ln*Te_over_Ti*R2*(kx(i, :)*Vstar1 + ky(i, :)*Vstar2 + kz(i, :)*Vstar3))/ &
                      (1.0_WP + ATeTi*kperp(i, :)**2)
         END DO

         DEALLOCATE (tmp4)

         w = REAL(USE_freq, wp)*w

      ELSE

         ph = 0.0_WP
         phx = 0.0_WP
         phy = 0.0_WP
         kx = 0.0_WP
         ky = 0.0_WP
         kz = 0.0_WP
         w = 0.0_WP
         wc = 0.0_WP
         kperp = 0.0_WP

      END IF

      ! First trajectory without turbulence
      ph(:, 1) = 0.0_WP; phx(:, 1) = 0.0_WP; phy(:, 1) = 0.0_WP
      kx(:, 1) = 0.0_WP; ky(:, 1) = 0.0_WP; kz(:, 1) = 0.0_WP
      w(:, 1) = 0.0_WP; wc(:, 1) = 0.0_WP; kperp(:, 1) = 0.0_WP

   END SUBROUTINE wavenum

   !=============================================================================
   ! wavenum_single (Nc-only mode; re-indented only)
   !=============================================================================
   SUBROUTINE wavenum_single(normB, Vstar1, Vstar2, Vstar3, gees)
      IMPLICIT NONE

      REAL(KIND=wp), DIMENSION(Nc) :: w1s, w2s
      REAL(KIND=wp), DIMENSION(5, Nc) :: g
      INTEGER :: i

      REAL(KIND=wp) :: normB, Vstar1, Vstar2, Vstar3
      REAL(KIND=wp), DIMENSION(3, 3) :: gees

      REAL(KIND=wp) :: w000, kx000, trunc

      !-----------------------------
      ! Allocate/resize arrays
      !-----------------------------
      IF (ALLOCATED(phs)) THEN
         IF (SIZE(phs) .NE. Nc) THEN
            DEALLOCATE (phs, kxs, kys, kperps, kzs, ws, phrmps, phxs, phys, wcs, dms)

            ALLOCATE (phs(Nc))
            ALLOCATE (kxs, kys, kperps, kzs, ws, phrmps, phxs, phys, wcs, dms, MOLD=phs)
         END IF
      ELSE
         ALLOCATE (phs(Nc))
         ALLOCATE (kxs, kys, kperps, kzs, ws, phrmps, phxs, phys, wcs, dms, MOLD=phs)
      END IF

      w000 = 0.0_WP
      kx000 = 0.0_WP
      trunc = 10.0_WP

      !-----------------------------
      ! Turbulence on/off
      !-----------------------------
      IF (USE_turb == ON) THEN

         CALL RANDOM_NUMBER(phs)
         CALL RANDOM_NUMBER(phrmps)
         CALL RANDOM_NUMBER(phxs)
         CALL RANDOM_NUMBER(phys)
         CALL RANDOM_NUMBER(dms)

         dms = 2.0_WP*(dms - 0.5_WP)
         dms = INT(q00/lambdaz*dms)

         phs = pi*(2.0_WP*phs - 1.0_WP)
         phrmps = pi*(2.0_WP*phrmps - 1.0_WP)
         phxs = pi*(2.0_WP*phxs - 1.0_WP)
         phys = pi*(2.0_WP*phys - 1.0_WP)

         ! Correlated (kx,kz,w,wc) as Nc samples
         CALL PDF_G(4, Nc, g(1:4, :), wavelength=[lambdax, lambdaz, tauc, taucc])

         kzs = g(2, :)
         wcs = g(4, :)

         IF (x_corr .EQ. 1) THEN
            kxs = g(1, :)
         ELSEIF (x_corr .EQ. 2) THEN
            CALL PDF_Cauchy(Np, kx000, lambdax, kxs, trunc)
         END IF

         IF (t_corr .EQ. 1) THEN
            w1s = g(3, :)
         ELSEIF (t_corr .EQ. 2) THEN
            CALL PDF_Cauchy(Np, w000, lambdax, w1s, trunc)
         END IF

         !-------------------------
         ! ITG: ky drawn with k0i
         !-------------------------
         CALL PDF_ky(Nc, lambday, k0i, g(5, :))

         DO i = 1, Nci
            kys(i) = g(5, i)
            kys(i) = 1.0_WP*INT(kys(i)*C2)/C2

            kz(i, :) = REAL(INT(kz(i, :), wp))/C3   ! kept as in original paste

            kperps(i) = (rhoi/R0)*SQRT( &
                        kxs(i)*gees(1, 1)*kxs(i) + 2.0_WP*kxs(i)*gees(1, 2)*kys(i) + &
                        2.0_WP*kxs(i)*gees(1, 3)*kzs(i) + kys(i)*gees(2, 2)*kys(i) + &
                        2.0_WP*kys(i)*gees(2, 3)*kzs(i) + kzs(i)*gees(3, 3)*kzs(i))

            w2s(i) = +Ln*(Ti/Ti)*((rhoi/R0/normB)**2)*(kxs(i)*Vstar1 + kys(i)*Vstar2 + kzs(i)*Vstar3)
            w2s(i) = w2s(i)/(1.0_WP + Aeff*(Te/Ti)*kperps(i)**2)

            ws(i) = w1s(i) + 1.0_WP*w2s(i)
         END DO

         !-------------------------
         ! TEM: ky drawn with k0e
         !-------------------------
         CALL PDF_ky(Nc, lambday, k0e, g(5, :))

         DO i = Nci + 1, Nc
            kys(i) = g(5, i)
            kys(i) = 1.0_WP*INT(kys(i)*C2)/C2

            kz(i, :) = REAL(INT(kz(i, :), wp))/C3   ! kept as in original paste

            kperps(i) = (rhoi/R0)*SQRT( &
                        kxs(i)*gees(1, 1)*kxs(i) + 2.0_WP*kxs(i)*gees(1, 2)*kys(i) + &
                        2.0_WP*kxs(i)*gees(1, 3)*kzs(i) + kys(i)*gees(2, 2)*kys(i) + &
                        2.0_WP*kys(i)*gees(2, 3)*kzs(i) + kzs(i)*gees(3, 3)*kzs(i))

            w2s(i) = -Ln*(Te/Ti)*((rhoi/R0/normB)**2)*(kxs(i)*Vstar1 + kys(i)*Vstar2 + kzs(i)*Vstar3)
            w2s(i) = w2s(i)/(1.0_WP + Aeff*(Te/Ti)*kperps(i)**2)

            ws(i) = w1s(i) + 1.0_WP*w2s(i)
         END DO

         ws = REAL(USE_freq, wp)*ws

      ELSEIF (USE_turb /= ON) THEN

         phs = 0.0_WP
         phxs = 0.0_WP
         phys = 0.0_WP
         kxs = 0.0_WP
         kys = 0.0_WP
         kzs = 0.0_WP
         ws = 0.0_WP
         wcs = 0.0_WP
         kperps = 0.0_WP

      END IF

      ! First sample without turbulence
      phs(1) = 0.0_WP; phxs(1) = 0.0_WP; phys(1) = 0.0_WP
      kxs(1) = 0.0_WP; kys(1) = 0.0_WP; kzs(1) = 0.0_WP
      ws(1) = 0.0_WP; wcs(1) = 0.0_WP
      kperps(1) = 0.0_WP; 
   END SUBROUTINE wavenum_single

   !=============================================================================
   ! Larmor_old (legacy variant)
   !=============================================================================
   SUBROUTINE Larmor_old(ind, X, Y, Z, mut)
      IMPLICIT NONE

      REAL(KIND=wp), DIMENSION(Np) :: X, Y, Z, mut, B
      INTEGER :: i, ind

      IF (ALLOCATED(L)) THEN
         DEALLOCATE (L)
      END IF
      ALLOCATE (L(Nc, Np))

      IF (ind == ON) THEN
         ! Initial |B| and rhot are evaluated by equilibrium_initialization_fields.

         !$OMP PARALLEL DO &
         !$OMP SHARED (mut,B) &
         !$OMP PRIVATE (i)
         DO i = 1, Np
            ! L(:,i) = BESSEL_J0(kperp(:,i)*SQRT(2.0*Aw*ABS(mut(i))/Zw**2/B(i)))
            L(:, i) = BESSEL_J0(kperp(:, i)*SQRT(2.0_WP*As*ABS(mut(i))/Zs**2))
         END DO
         !$OMP END PARALLEL DO

      ELSE
         L = 1.0_WP
      END IF

   END SUBROUTINE Larmor_old

END MODULE random_numbers
