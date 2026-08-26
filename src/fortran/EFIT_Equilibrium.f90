!========================================================================================================
! Subroutine :: interpol (EFIT2 / psisurf_efit2 / psisurf_efit3)
! Purpose    :: interpolate quantities known on a 2D grid (EFIT-like data)
!
! Record of revisions:
!    Date       | Programmer          | Description
!  -----------  | ------------------- | -----------------------------------------
!  --/--/2021   | D. I. Palade        | Initial version
!  --/--/2023   | L. M. Pomarjanschi  | Corrections and development
!  01/01/2024   | D. I. Palade        | Remake
!========================================================================================================

!========================================================================================================
! EFIT2: interpolate all Nqua quantities using bilinear interpolation on the EFIT grid
!========================================================================================================
SUBROUTINE EFIT2(npoints, X, Y, R)
   USE constants
   IMPLICIT NONE

   !---------------------------------------------------------------------------------
   ! I/O
   !---------------------------------------------------------------------------------
   INTEGER, INTENT(IN) :: npoints
   REAL(KIND=wp), DIMENSION(npoints), INTENT(IN) :: X, Y
   REAL(KIND=wp), DIMENSION(Nqua, npoints), INTENT(OUT) :: R

   !---------------------------------------------------------------------------------
   ! Locals
   !---------------------------------------------------------------------------------
   INTEGER, DIMENSION(npoints) :: poz1, poz2
   INTEGER :: i, j
   REAL(KIND=wp), DIMENSION(npoints) :: F1, F2, F3, F4
   REAL(KIND=wp), DIMENSION(npoints) :: X1, Y1, Xef, Yef

   !---------------------------------------------------------------------------------
   ! Cell indices and local coordinates
   !---------------------------------------------------------------------------------
   poz1 = modulo(int((X - minR)/stepR), NgridR)
   poz2 = modulo(int((Y - minZ)/stepZ), NgridZ)

   X1 = minR + poz1*stepR
   Y1 = minZ + poz2*stepZ

   Xef = (X - X1)/stepR
   Yef = (Y - Y1)/stepZ

   !---------------------------------------------------------------------------------
   ! Bilinear interpolation for each quantity i
   !---------------------------------------------------------------------------------
   DO i = 1, Nqua
      DO j = 1, npoints
         F1(j) = Efit_data(NgridR*NgridZ*(i - 1) + poz1(j)*NgridR + poz2(j))
         F2(j) = Efit_data(NgridR*NgridZ*(i - 1) + (poz1(j) + 1)*NgridR + poz2(j))
         F3(j) = Efit_data(NgridR*NgridZ*(i - 1) + poz1(j)*NgridR + (poz2(j) + 1))
         F4(j) = Efit_data(NgridR*NgridZ*(i - 1) + (poz1(j) + 1)*NgridR + (poz2(j) + 1))
      END DO

      R(i, :) = F1 + (F2 - F1)*Xef + (F3 - F1)*Yef + Xef*Yef*(F1 + F4 - F2 - F3)
   END DO

END SUBROUTINE EFIT2

!========================================================================================================
! psisurf_efit3: contour-based sampling of phi(R,Z)=0 where
!                phi = psi - psi0 - K*(F/normB)*Vp(k)
!========================================================================================================
SUBROUTINE psisurf_efit3(X, Y, Vp)
   USE constants
   IMPLICIT NONE

   !---------------------------------------------------------------------------------
   ! Arguments
   !---------------------------------------------------------------------------------
   REAL(KIND=wp), DIMENSION(Np), INTENT(OUT) :: X, Y
   REAL(KIND=wp), DIMENSION(Np), INTENT(IN) :: Vp

   !---------------------------------------------------------------------------------
   ! Locals
   !---------------------------------------------------------------------------------
   INTEGER :: i, j, k, npts, idx
   REAL(KIND=wp), DIMENSION(NgridR, NgridZ) :: psi, F, psir, psiz, normB, phi3
   REAL(KIND=wp), DIMENSION(NgridR) :: Rvals
   REAL(KIND=wp), DIMENSION(NgridZ) :: Zvals

   REAL(KIND=wp) :: dR_loc, dZ_loc, Kconst, urand
   REAL(KIND=wp), ALLOCATABLE :: Rc(:), Zc(:)

   !---------------------------------------------------------------------------------
   ! Grid vectors
   !---------------------------------------------------------------------------------
   dR_loc = stepR
   dZ_loc = stepZ

   Rvals = [(minR + (i - 1)*dR_loc, i=1, NgridR)]
   Zvals = [(minZ + (j - 1)*dZ_loc, j=1, NgridZ)]

   !---------------------------------------------------------------------------------
   ! Rebuild psi and related quantities from Efit_data
   !---------------------------------------------------------------------------------
   DO j = 1, NgridZ
      DO i = 1, NgridR
         idx = (i - 1)*NgridR + j
         psi(i, j) = Efit_data(idx)
         F(i, j) = Efit_data(idx + 6*NgridR*NgridZ)
         psir(i, j) = Efit_data(idx + 1*NgridR*NgridZ)
         psiz(i, j) = Efit_data(idx + 2*NgridR*NgridZ)
         normB(i, j) = sqrt((F(i, j)**2 + psir(i, j)**2 + psiz(i, j)**2)/(Rvals(i)**2))
      END DO
   END DO

   ! Constant prefactor K = (USE_PC*rhoi/R0)*(As/Zs)
   Kconst = REAL(USE_PC, wp)*rhoi/R0*As/Zs

   !---------------------------------------------------------------------------------
   ! Parallel loop over particles: build phi3 and sample from phi3=0 contour intersections
   !---------------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k, phi3, Rc, Zc, npts, urand, idx) SCHEDULE(DYNAMIC)
   DO k = 1, Np

      phi3 = psi - psi0 - Kconst*(F/normB)*Vp(k)

      CALL extract_contour_zero(phi3, Rvals, Zvals, Rc, Zc, npts)

      IF (npts < 1) THEN
         X(k) = Rvals(MAX(1, NgridR/2))
         Y(k) = Zvals(MAX(1, NgridZ/2))
         CYCLE
      END IF

      CALL random_number(urand)
      idx = INT(1 + urand*REAL(npts - 1, wp))   ! random intersection

      X(k) = Rc(idx)
      Y(k) = Zc(idx)

      DEALLOCATE (Rc, Zc)
   END DO
!$OMP END PARALLEL DO

CONTAINS

   !---------------------------------------------------------------------------------
   ! Extract intersections of phi=0 with grid edges; returns points (unordered).
   !---------------------------------------------------------------------------------
   SUBROUTINE extract_contour_zero(phi3, Rvals, Zvals, Rc, Zc, npts)
      USE constants
      IMPLICIT NONE

      REAL(KIND=wp), INTENT(IN) :: phi3(NgridR, NgridZ), Rvals(NgridR), Zvals(NgridZ)
      REAL(KIND=wp), ALLOCATABLE, INTENT(OUT) :: Rc(:), Zc(:)
      INTEGER, INTENT(OUT) :: npts

      INTEGER :: i, j, icount
      REAL(wp) :: f00, f10, f01, f11, t, rint, zint
      REAL(wp), ALLOCATABLE :: rtemp(:), ztemp(:)

      ALLOCATE (rtemp(NgridR*NgridZ*4))
      ALLOCATE (ztemp(NgridR*NgridZ*4))
      icount = 0

      DO j = 1, NgridZ - 1
         DO i = 1, NgridR - 1
            f00 = phi3(i, j)
            f10 = phi3(i + 1, j)
            f01 = phi3(i, j + 1)
            f11 = phi3(i + 1, j + 1)

            ! Bottom edge: (i,j) -> (i+1,j)
            IF (f00*f10 < 0.0_WP) THEN
               t = f00/(f00 - f10)
               rint = Rvals(i) + t*(Rvals(i + 1) - Rvals(i))
               zint = Zvals(j)
               icount = icount + 1
               rtemp(icount) = rint
               ztemp(icount) = zint
            END IF

            ! Top edge: (i,j+1) -> (i+1,j+1)
            IF (f01*f11 < 0.0_WP) THEN
               t = f01/(f01 - f11)
               rint = Rvals(i) + t*(Rvals(i + 1) - Rvals(i))
               zint = Zvals(j + 1)
               icount = icount + 1
               rtemp(icount) = rint
               ztemp(icount) = zint
            END IF

            ! Left edge: (i,j) -> (i,j+1)
            IF (f00*f01 < 0.0_WP) THEN
               t = f00/(f00 - f01)
               rint = Rvals(i)
               zint = Zvals(j) + t*(Zvals(j + 1) - Zvals(j))
               icount = icount + 1
               rtemp(icount) = rint
               ztemp(icount) = zint
            END IF

            ! Right edge: (i+1,j) -> (i+1,j+1)
            IF (f10*f11 < 0.0_WP) THEN
               t = f10/(f10 - f11)
               rint = Rvals(i + 1)
               zint = Zvals(j) + t*(Zvals(j + 1) - Zvals(j))
               icount = icount + 1
               rtemp(icount) = rint
               ztemp(icount) = zint
            END IF

         END DO
      END DO

      IF (icount > 1) THEN
         npts = icount
         ALLOCATE (Rc(npts), Zc(npts))
         Rc = rtemp(1:npts)
         Zc = ztemp(1:npts)
      ELSE
         npts = 0
         ALLOCATE (Rc(0), Zc(0))
      END IF

      DEALLOCATE (rtemp, ztemp)
   END SUBROUTINE extract_contour_zero

   !---------------------------------------------------------------------------------
   ! The routines below are kept (reformatted) for completeness. They are not called
   ! from the current psisurf_efit3 implementation.
   !---------------------------------------------------------------------------------

   SUBROUTINE order_points_by_angle(Rc, Zc)
      USE constants
      IMPLICIT NONE
      REAL(KIND=wp), INTENT(INOUT) :: Rc(:), Zc(:)

      INTEGER :: n, i, j, itmp
      REAL(wp) :: cx, cz, tmpa, tmpr, tmpz
      REAL(wp), ALLOCATABLE :: ang(:)
      INTEGER, ALLOCATABLE :: idx(:)

      n = SIZE(Rc)
      IF (n <= 2) RETURN

      cx = SUM(Rc)/n
      cz = SUM(Zc)/n

      ALLOCATE (ang(n), idx(n))
      DO i = 1, n
         ang(i) = atan2(Zc(i) - cz, Rc(i) - cx)
         idx(i) = i
      END DO

      ! Insertion sort by angle (in-place)
      DO i = 2, n
         tmpa = ang(i)
         itmp = idx(i)
         tmpr = Rc(i)
         tmpz = Zc(i)

         j = i - 1
         DO WHILE (j >= 1 .AND. ang(j) > tmpa)
            ang(j + 1) = ang(j)
            idx(j + 1) = idx(j)
            Rc(j + 1) = Rc(j)
            Zc(j + 1) = Zc(j)
            j = j - 1
         END DO

         ang(j + 1) = tmpa
         idx(j + 1) = itmp
         Rc(j + 1) = tmpr
         Zc(j + 1) = tmpz
      END DO

      DEALLOCATE (ang, idx)
   END SUBROUTINE order_points_by_angle

   SUBROUTINE interp_point_along_contour(Rc, Zc, s, starget, Rout, Zout)
      USE constants
      IMPLICIT NONE

      REAL(KIND=wp), INTENT(IN) :: Rc(:), Zc(:), s(:), starget
      REAL(KIND=wp), INTENT(OUT) :: Rout, Zout

      INTEGER :: i

      DO i = 2, SIZE(s)
         IF (starget <= s(i)) THEN
            Rout = Rc(i - 1) + (Rc(i) - Rc(i - 1))*(starget - s(i - 1))/(s(i) - s(i - 1))
            Zout = Zc(i - 1) + (Zc(i) - Zc(i - 1))*(starget - s(i - 1))/(s(i) - s(i - 1))
            RETURN
         END IF
      END DO

      Rout = Rc(SIZE(Rc))
      Zout = Zc(SIZE(Zc))
   END SUBROUTINE interp_point_along_contour

END SUBROUTINE psisurf_efit3
