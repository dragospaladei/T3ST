!========================================================================================================
! Subroutine :: dispersion
! Purpose    :: Computes the quantities needed for the evaluation of mode frequencies
!               w(k) = Vstar*k / (1 + ceva*kperp^2)
!
! Record of revisions:
!    Date       | Programmer          | Description
!  -----------  | ------------------- | -----------------------------------------
!  --/--/2021   | D. I. Palade        | Initial version
!  --/--/2023   | L. M. Pomarjanschi  | Corrections and development
!  01/01/2024   | D. I. Palade        | Remake
!========================================================================================================

SUBROUTINE dispersion(normB, Vstar1, Vstar2, Vstar3, gees)
   USE constants
   IMPLICIT NONE

   !---------------------------------------------------------------------------------
   ! I/O
   !---------------------------------------------------------------------------------
   REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT) :: gees
   REAL(KIND=dp), INTENT(OUT)                 :: Vstar1, Vstar2, Vstar3, normB

   !---------------------------------------------------------------------------------
   ! Locals
   !---------------------------------------------------------------------------------
   REAL(KIND=dp), DIMENSION(Np)       :: hx, hy, hz, q1, q2, q3
   REAL(KIND=dp), DIMENSION(Np)       :: X1, Y1, Z1
   REAL(KIND=dp), DIMENSION(Np)       :: Bx, By, Bz, gradBx, gradBy, gradBz
   REAL(KIND=dp), DIMENSION(Np)       :: gradu2x, gradu2y, gradu2z
   REAL(KIND=dp), DIMENSION(Np)       :: rotbx, rotby, rotbz
   REAL(KIND=dp), DIMENSION(Np)       :: rotux, rotuy, rotuz
   REAL(KIND=dp), DIMENSION(Np)       :: phi0x, phi0y, phi0z
   REAL(KIND=dp), DIMENSION(Np)       :: B

   REAL(KIND=dp), DIMENSION(3, 3, Np) :: GG
   REAL(KIND=dp), DIMENSION(3)        :: lame, beta
   REAL(KIND=dp), DIMENSION(3, 3)     :: gege, kron

   INTEGER :: i, j, k, p

   !---------------------------------------------------------------------------------
   ! Reference point for geometry evaluation
   !---------------------------------------------------------------------------------
   X1 = X0
   Y1 = Y0
   Z1 = Z0

   hx = 1.0_dp
   hy = 1.0_dp
   hz = X1

   CALL Equil(X1, Y1, Z1, Bx, By, Bz, gradBx, gradBy, gradBz,                          &
              gradu2x, gradu2y, gradu2z, rotbx, rotby, rotbz,                         &
              rotux, rotuy, rotuz, phi0x, phi0y, phi0z, B, GG, q1, q2, q3)

   !---------------------------------------------------------------------------------
   ! Scalars at the reference point (index 1)
   !---------------------------------------------------------------------------------
   normB = B(1)                        ! |B| from covariant components
   lame  = (/ hx(1), hy(1), hz(1) /)
   beta  = (/ Bx(1), By(1), Bz(1) /)

   gege = GG(:, :, 1)

   kron = RESHAPE((/ 1.0_dp, 0.0_dp, 0.0_dp, &
                    0.0_dp, 1.0_dp, 0.0_dp, &
                    0.0_dp, 0.0_dp, 1.0_dp /), SHAPE(kron))

   !---------------------------------------------------------------------------------
   ! Diamagnetic velocity projections (Vstar)
   !---------------------------------------------------------------------------------
   Vstar1 = 0.0_dp

   Vstar2 =  Bx(1)*hx(1) * (GG(2, 2, 1)*GG(3, 1, 1) - GG(3, 2, 1)*GG(2, 1, 1)) / (hy(1)*hz(1))
   Vstar2 =  Vstar2 + By(1)*hy(1) * (GG(3, 2, 1)*GG(1, 1, 1) - GG(1, 2, 1)*GG(3, 1, 1)) / (hx(1)*hz(1))
   Vstar2 =  Vstar2 + Bz(1)*hz(1) * (GG(1, 2, 1)*GG(2, 1, 1) - GG(2, 2, 1)*GG(1, 1, 1)) / (hx(1)*hy(1))

   Vstar3 =  Bx(1)*hx(1) * (GG(2, 3, 1)*GG(3, 1, 1) - GG(3, 3, 1)*GG(2, 1, 1)) / (hy(1)*hz(1))
   Vstar3 =  Vstar3 + By(1)*hy(1) * (GG(3, 3, 1)*GG(1, 1, 1) - GG(1, 3, 1)*GG(3, 1, 1)) / (hx(1)*hz(1))
   Vstar3 =  Vstar3 + Bz(1)*hz(1) * (GG(1, 3, 1)*GG(2, 1, 1) - GG(2, 3, 1)*GG(1, 1, 1)) / (hx(1)*hy(1))

   !---------------------------------------------------------------------------------
   ! Matrix for kperp^2 evaluation: k_perp^2 = k_i * gees(i,j) * k_j
   !---------------------------------------------------------------------------------
   gees = 0.0_dp

   DO i = 1, 3
      DO j = 1, 3
         DO k = 1, 3
            DO p = 1, 3
               gees(i, j) = gees(i, j) + gege(k, i) * gege(p, j) * &
                            (kron(k, p)/(lame(k)*lame(p)) - beta(k)*beta(p)/normB**2)
            END DO
         END DO
      END DO
   END DO

END SUBROUTINE dispersion

