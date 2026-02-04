!========================================================================================================
! Subroutine :: interpol (EFIT2 / EFITsa / psisurf_efit2 / psisurf_efit3)
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
SUBROUTINE EFIT2(X, Y, R)
   USE constants
   IMPLICIT NONE

   !---------------------------------------------------------------------------------
   ! I/O
   !---------------------------------------------------------------------------------
   REAL(KIND=dp), DIMENSION(Np),      INTENT(IN)  :: X, Y
   REAL(KIND=dp), DIMENSION(Nqua, Np), INTENT(OUT) :: R

   !---------------------------------------------------------------------------------
   ! Locals
   !---------------------------------------------------------------------------------
   INTEGER,       DIMENSION(Np) :: poz1, poz2
   INTEGER                    :: i, j
   REAL(KIND=dp), DIMENSION(Np) :: F1, F2, F3, F4
   REAL(KIND=dp), DIMENSION(Np) :: X1, Y1, Xef, Yef

   !---------------------------------------------------------------------------------
   ! Cell indices and local coordinates
   !---------------------------------------------------------------------------------
   poz1 = modulo(int((X - minR)/stepR), NgridR)
   poz2 = modulo(int((Y - minZ)/stepZ), NgridZ)

   X1  = minR + poz1*stepR
   Y1  = minZ + poz2*stepZ

   Xef = (X - X1)/stepR
   Yef = (Y - Y1)/stepZ

   !---------------------------------------------------------------------------------
   ! Bilinear interpolation for each quantity i
   !---------------------------------------------------------------------------------
   DO i = 1, Nqua
      DO j = 1, Np
         F1(j) = Efit_data(NgridR*NgridZ*(i - 1) +  poz1(j)   *NgridR +  poz2(j))
         F2(j) = Efit_data(NgridR*NgridZ*(i - 1) + (poz1(j)+1)*NgridR +  poz2(j))
         F3(j) = Efit_data(NgridR*NgridZ*(i - 1) +  poz1(j)   *NgridR + (poz2(j)+1))
         F4(j) = Efit_data(NgridR*NgridZ*(i - 1) + (poz1(j)+1)*NgridR + (poz2(j)+1))
      END DO

      R(i, :) = F1 + (F2 - F1)*Xef + (F3 - F1)*Yef + Xef*Yef*(F1 + F4 - F2 - F3)
   END DO

END SUBROUTINE EFIT2


!========================================================================================================
! EFITsa: EFIT interpolation but overwrite (chi, chir, chiz) with circular-model expressions
!
! Note: using this model might lead to problems; observed that the theta-derivative of chi in
!       EFIT data is not always the same sign as in the circular model.
!========================================================================================================
SUBROUTINE EFITsa(X, Y, R)
   USE constants
   IMPLICIT NONE

   !---------------------------------------------------------------------------------
   ! I/O
   !---------------------------------------------------------------------------------
   REAL(KIND=dp), DIMENSION(Np),      INTENT(IN)  :: X, Y
   REAL(KIND=dp), DIMENSION(Nqua, Np), INTENT(OUT) :: R

   !---------------------------------------------------------------------------------
   ! Locals
   !---------------------------------------------------------------------------------
   INTEGER,       DIMENSION(Np) :: poz1, poz2
   INTEGER                    :: i, j
   REAL(KIND=dp), DIMENSION(Np) :: F1, F2, F3, F4
   REAL(KIND=dp), DIMENSION(Np) :: X1, Y1, Xef, Yef
   REAL(KIND=dp), DIMENSION(Np) :: rr, theta, chi, chir, chiz

   !---------------------------------------------------------------------------------
   ! Cell indices and local coordinates
   !---------------------------------------------------------------------------------
   poz1 = modulo(int((X - minR)/stepR), NgridR)
   poz2 = modulo(int((Y - minZ)/stepZ), NgridZ)

   X1  = minR + poz1*stepR
   Y1  = minZ + poz2*stepZ

   Xef = (X - X1)/stepR
   Yef = (Y - Y1)/stepZ

   !---------------------------------------------------------------------------------
   ! Interpolate everything except the last 3 quantities (chi,chir,chiz)
   !---------------------------------------------------------------------------------
   DO i = 1, Nqua - 3
      DO j = 1, Np
         F1(j) = Efit_data(NgridR*NgridZ*(i - 1) +  poz1(j)   *NgridR +  poz2(j))
         F2(j) = Efit_data(NgridR*NgridZ*(i - 1) + (poz1(j)+1)*NgridR +  poz2(j))
         F3(j) = Efit_data(NgridR*NgridZ*(i - 1) +  poz1(j)   *NgridR + (poz2(j)+1))
         F4(j) = Efit_data(NgridR*NgridZ*(i - 1) + (poz1(j)+1)*NgridR + (poz2(j)+1))
      END DO

      R(i, :) = F1 + (F2 - F1)*Xef + (F3 - F1)*Yef + Xef*Yef*(F1 + F4 - F2 - F3)
   END DO

   !---------------------------------------------------------------------------------
   ! Overwrite chi, chir, chiz with circular expressions
   !---------------------------------------------------------------------------------
   rr    = sqrt((X - 1.0_dp)**2 + Y**2) + 1.0e-7_dp
   theta = atan2(Y, X - 1.0_dp)

   chi  = 2.0_dp*atan( sqrt((1.0_dp - rr)/(1.0_dp + rr)) * tan(theta/2.0_dp) )
   chir =  sin(theta) * (rr**2 - X) / X / rr / sqrt(1.0_dp - rr**2)
   chiz = -(rr - cos(theta)) / rr / sqrt(1.0_dp - rr**2)

   R(Nqua - 2, :) = chi
   R(Nqua - 1, :) = chir
   R(Nqua    , :) = chiz

END SUBROUTINE EFITsa


!========================================================================================================
! psisurf_efit2: choose (X,Y) on a psi-surface satisfying a condition involving Vp
!========================================================================================================
SUBROUTINE psisurf_efit_old(X, Y, Vp)
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
   REAL(KIND=dp), DIMENSION(Np)      :: aux, keep
   INTEGER,       DIMENSION(NgridR)  :: aux1

   LOGICAL                          :: mask(NgridR)
   REAL(KIND=dp), DIMENSION(NgridR) :: Rvals, Zvals, diff, targeta
   INTEGER                          :: i, j, k, ioc, joc
   REAL(KIND=dp), DIMENSION(NgridR) :: F, psi, psir, psiz, normB

   ! Random line angle
   CALL random_number(aux)
   aux = pi*(2.0_dp*aux - 1.0_dp)

   DO k = 1, Np

      ! Build indices along a straight line in (R,Z) via tan(theta)
      DO i = 1, NgridR
         j = int(-minZ/stepZ + tan(aux(k))*(minR - 1.0_dp + i*stepR)/stepZ)
         aux1(i) = i*NgridR + j
      END DO

      ! Sample fields along that line
      DO i = 1, NgridR
         IF ((aux1(i) > (NgridR*NgridZ)) .OR. (aux1(i) <= 0)) THEN
            aux1(i) = int(Ngrid/2)
         END IF

         psi(i)  = Efit_data(aux1(i))
         F(i)    = Efit_data(aux1(i) + 6*NgridR*NgridZ)
         psir(i) = Efit_data(aux1(i) + 1*NgridR*NgridZ)
         psiz(i) = Efit_data(aux1(i) + 2*NgridR*NgridZ)

         normB(i) = sqrt((F(i)**2 + psir(i)**2 + psiz(i)**2) / (minR + i*stepR)**2)
      END DO

      ! R and Z samples along that line
      Rvals = [(minR + (i-1)*stepR, i=1, NgridR)]
      Zvals = tan(aux(k))*(Rvals - 1.0_dp)

      ! Target psi along the line
      targeta = real(USE_PC,dp)*rhoi/R0*Aw/Zw*F/normB*Vp(k) + psi0

      ! Mask inside plasma region
      mask = ((Rvals - 1.0_dp)**2 + Zvals**2 <= (a0**2))

      ! Difference array and minimum restricted to mask
      diff = abs(psi - targeta)
      ioc  = minloc(merge(diff, huge(100.0_dp), mask), 1)

      joc = int(-minZ/stepZ + tan(aux(k))*(minR - 1.0_dp + ioc*stepR)/stepZ)

      keep(k) = ioc
      X(k)    = minR + ioc*stepR
      Y(k)    = minZ + joc*stepZ
   END DO

END SUBROUTINE psisurf_efit_old


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
   REAL(KIND=dp), DIMENSION(Np), INTENT(OUT) :: X, Y
   REAL(KIND=dp), DIMENSION(Np), INTENT(IN)  :: Vp

   !---------------------------------------------------------------------------------
   ! Locals
   !---------------------------------------------------------------------------------
   INTEGER :: i, j, k, npts, idx
   REAL(KIND=dp), DIMENSION(NgridR, NgridZ) :: psi, F, psir, psiz, normB, phi3
   REAL(KIND=dp), DIMENSION(NgridR)         :: Rvals
   REAL(KIND=dp), DIMENSION(NgridZ)         :: Zvals

   REAL(KIND=dp) :: dR_loc, dZ_loc, Kconst, urand
   REAL(KIND=dp), ALLOCATABLE :: Rc(:), Zc(:)

   !---------------------------------------------------------------------------------
   ! Grid vectors
   !---------------------------------------------------------------------------------
   dR_loc = stepR
   dZ_loc = stepZ

   Rvals = [(minR + (i-1)*dR_loc, i=1, NgridR)]
   Zvals = [(minZ + (j-1)*dZ_loc, j=1, NgridZ)]

   !---------------------------------------------------------------------------------
   ! Rebuild psi and related quantities from Efit_data
   !---------------------------------------------------------------------------------
   DO j = 1, NgridZ
      DO i = 1, NgridR
         idx = (i-1)*NgridR + j
         psi(i, j)   = Efit_data(idx)
         F(i, j)     = Efit_data(idx + 6*NgridR*NgridZ)
         psir(i, j)  = Efit_data(idx + 1*NgridR*NgridZ)
         psiz(i, j)  = Efit_data(idx + 2*NgridR*NgridZ)
         normB(i, j) = sqrt((F(i, j)**2 + psir(i, j)**2 + psiz(i, j)**2) / (Rvals(i)**2))
      END DO
   END DO

   ! Constant prefactor K = (USE_PC*rhoi/R0)*(Aw/Zw)
   Kconst = real(USE_PC,dp)*rhoi/R0*Aw/Zw

   !---------------------------------------------------------------------------------
   ! Parallel loop over particles: build phi3 and sample from phi3=0 contour intersections
   !---------------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k, phi3, Rc, Zc, npts, urand, idx) SCHEDULE(DYNAMIC)
   DO k = 1, Np

      phi3 = psi - psi0 - Kconst * (F/normB) * Vp(k)

      CALL extract_contour_zero(phi3, Rvals, Zvals, Rc, Zc, npts)

      IF (npts < 1) THEN
         X(k) = Rvals(MAX(1, NgridR/2))
         Y(k) = Zvals(MAX(1, NgridZ/2))
         CYCLE
      END IF

      CALL random_number(urand)
      idx = INT(1 + urand * REAL(npts-1, dp))   ! random intersection

      X(k) = Rc(idx)
      Y(k) = Zc(idx)

      DEALLOCATE(Rc, Zc)
   END DO
!$OMP END PARALLEL DO

CONTAINS

   !---------------------------------------------------------------------------------
   ! Extract intersections of phi=0 with grid edges; returns points (unordered).
   !---------------------------------------------------------------------------------
   SUBROUTINE extract_contour_zero(phi3, Rvals, Zvals, Rc, Zc, npts)
      USE constants
      IMPLICIT NONE

      REAL(KIND=dp), INTENT(IN)  :: phi3(NgridR, NgridZ), Rvals(NgridR), Zvals(NgridZ)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: Rc(:), Zc(:)
      INTEGER, INTENT(OUT) :: npts

      INTEGER  :: i, j, icount
      REAL(dp) :: f00, f10, f01, f11, t, rint, zint
      REAL(dp), ALLOCATABLE :: rtemp(:), ztemp(:)

      ALLOCATE(rtemp(NgridR*NgridZ*4))
      ALLOCATE(ztemp(NgridR*NgridZ*4))
      icount = 0

      DO j = 1, NgridZ - 1
         DO i = 1, NgridR - 1
            f00 = phi3(i  , j  )
            f10 = phi3(i+1, j  )
            f01 = phi3(i  , j+1)
            f11 = phi3(i+1, j+1)

            ! Bottom edge: (i,j) -> (i+1,j)
            IF (f00*f10 < 0.0_dp) THEN
               t    = f00 / (f00 - f10)
               rint = Rvals(i) + t*(Rvals(i+1) - Rvals(i))
               zint = Zvals(j)
               icount = icount + 1
               rtemp(icount) = rint
               ztemp(icount) = zint
            END IF

            ! Top edge: (i,j+1) -> (i+1,j+1)
            IF (f01*f11 < 0.0_dp) THEN
               t    = f01 / (f01 - f11)
               rint = Rvals(i) + t*(Rvals(i+1) - Rvals(i))
               zint = Zvals(j+1)
               icount = icount + 1
               rtemp(icount) = rint
               ztemp(icount) = zint
            END IF

            ! Left edge: (i,j) -> (i,j+1)
            IF (f00*f01 < 0.0_dp) THEN
               t    = f00 / (f00 - f01)
               rint = Rvals(i)
               zint = Zvals(j) + t*(Zvals(j+1) - Zvals(j))
               icount = icount + 1
               rtemp(icount) = rint
               ztemp(icount) = zint
            END IF

            ! Right edge: (i+1,j) -> (i+1,j+1)
            IF (f10*f11 < 0.0_dp) THEN
               t    = f10 / (f10 - f11)
               rint = Rvals(i+1)
               zint = Zvals(j) + t*(Zvals(j+1) - Zvals(j))
               icount = icount + 1
               rtemp(icount) = rint
               ztemp(icount) = zint
            END IF

         END DO
      END DO

      IF (icount > 1) THEN
         npts = icount
         ALLOCATE(Rc(npts), Zc(npts))
         Rc = rtemp(1:npts)
         Zc = ztemp(1:npts)
      ELSE
         npts = 0
         ALLOCATE(Rc(0), Zc(0))
      END IF

      DEALLOCATE(rtemp, ztemp)
   END SUBROUTINE extract_contour_zero


   !---------------------------------------------------------------------------------
   ! The routines below are kept (reformatted) for completeness. They are not called
   ! from the current psisurf_efit3 implementation.
   !---------------------------------------------------------------------------------

   SUBROUTINE order_points_by_angle(Rc, Zc)
      USE constants
      IMPLICIT NONE
      REAL(KIND=dp), INTENT(INOUT) :: Rc(:), Zc(:)

      INTEGER :: n, i, j, itmp
      REAL(dp) :: cx, cz, tmpa, tmpr, tmpz
      REAL(dp), ALLOCATABLE :: ang(:)
      INTEGER,  ALLOCATABLE :: idx(:)

      n = SIZE(Rc)
      IF (n <= 2) RETURN

      cx = SUM(Rc)/n
      cz = SUM(Zc)/n

      ALLOCATE(ang(n), idx(n))
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
            ang(j+1) = ang(j)
            idx(j+1) = idx(j)
            Rc(j+1)  = Rc(j)
            Zc(j+1)  = Zc(j)
            j = j - 1
         END DO

         ang(j+1) = tmpa
         idx(j+1) = itmp
         Rc(j+1)  = tmpr
         Zc(j+1)  = tmpz
      END DO

      DEALLOCATE(ang, idx)
   END SUBROUTINE order_points_by_angle


   SUBROUTINE interp_point_along_contour(Rc, Zc, s, starget, Rout, Zout)
      USE constants
      IMPLICIT NONE

      REAL(KIND=dp), INTENT(IN)  :: Rc(:), Zc(:), s(:), starget
      REAL(KIND=dp), INTENT(OUT) :: Rout, Zout

      INTEGER :: i

      DO i = 2, SIZE(s)
         IF (starget <= s(i)) THEN
            Rout = Rc(i-1) + (Rc(i) - Rc(i-1))*(starget - s(i-1))/(s(i) - s(i-1))
            Zout = Zc(i-1) + (Zc(i) - Zc(i-1))*(starget - s(i-1))/(s(i) - s(i-1))
            RETURN
         END IF
      END DO

      Rout = Rc(SIZE(Rc))
      Zout = Zc(SIZE(Zc))
   END SUBROUTINE interp_point_along_contour

END SUBROUTINE psisurf_efit3

