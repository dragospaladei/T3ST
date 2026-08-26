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
   REAL(KIND=wp), DIMENSION(3, 3), INTENT(OUT) :: gees
   REAL(KIND=wp), INTENT(OUT) :: Vstar1, Vstar2, Vstar3, normB

   !---------------------------------------------------------------------------------
   ! Locals
   !---------------------------------------------------------------------------------
   REAL(KIND=wp), DIMENSION(Np) :: hx, hy, hz
!  REAL(KIND=wp), DIMENSION(Np)       :: q1, q2, q3
   REAL(KIND=wp), DIMENSION(Np) :: X1, Y1, Z1
   REAL(KIND=wp), DIMENSION(Np) :: Bx, By, Bz
!  REAL(KIND=wp), DIMENSION(Np)       :: gradBx, gradBy, gradBz
!  REAL(KIND=wp), DIMENSION(Np)       :: gradu2x, gradu2y, gradu2z
!  REAL(KIND=wp), DIMENSION(Np)       :: rotbx, rotby, rotbz
!  REAL(KIND=wp), DIMENSION(Np)       :: rotux, rotuy, rotuz
!  REAL(KIND=wp), DIMENSION(Np)       :: phi0x, phi0y, phi0z
   REAL(KIND=wp), DIMENSION(Np) :: B

   REAL(KIND=wp), DIMENSION(3, 3, Np) :: GG
   REAL(KIND=wp), DIMENSION(3) :: lame, beta_vec
   REAL(KIND=wp), DIMENSION(3, 3) :: gege, kron

   INTEGER :: i, j, k, p

   !---------------------------------------------------------------------------------
   ! Reference point for geometry evaluation
   !---------------------------------------------------------------------------------
   X1 = X0
   Y1 = Y0
   Z1 = Z0

   hx = 1.0_WP
   hy = 1.0_WP
   hz = X1

   CALL Equil(X1, Y1, Z1, Bx, By, Bz, B, GG)

   !---------------------------------------------------------------------------------
   ! Scalars at the reference point (index 1)
   !---------------------------------------------------------------------------------
   normB = B(1)                        ! |B| from covariant components
   lame = (/hx(1), hy(1), hz(1)/)
   beta_vec = (/Bx(1), By(1), Bz(1)/)

   gege = GG(:, :, 1)

   kron = RESHAPE((/1.0_WP, 0.0_WP, 0.0_WP, &
                    0.0_WP, 1.0_WP, 0.0_WP, &
                    0.0_WP, 0.0_WP, 1.0_WP/), SHAPE(kron))

   !---------------------------------------------------------------------------------
   ! Diamagnetic velocity projections (Vstar)
   !---------------------------------------------------------------------------------
   Vstar1 = 0.0_WP

   Vstar2 = Bx(1)*hx(1)*(GG(2, 2, 1)*GG(3, 1, 1) - GG(3, 2, 1)*GG(2, 1, 1))/(hy(1)*hz(1))
   Vstar2 = Vstar2 + By(1)*hy(1)*(GG(3, 2, 1)*GG(1, 1, 1) - GG(1, 2, 1)*GG(3, 1, 1))/(hx(1)*hz(1))
   Vstar2 = Vstar2 + Bz(1)*hz(1)*(GG(1, 2, 1)*GG(2, 1, 1) - GG(2, 2, 1)*GG(1, 1, 1))/(hx(1)*hy(1))

   Vstar3 = Bx(1)*hx(1)*(GG(2, 3, 1)*GG(3, 1, 1) - GG(3, 3, 1)*GG(2, 1, 1))/(hy(1)*hz(1))
   Vstar3 = Vstar3 + By(1)*hy(1)*(GG(3, 3, 1)*GG(1, 1, 1) - GG(1, 3, 1)*GG(3, 1, 1))/(hx(1)*hz(1))
   Vstar3 = Vstar3 + Bz(1)*hz(1)*(GG(1, 3, 1)*GG(2, 1, 1) - GG(2, 3, 1)*GG(1, 1, 1))/(hx(1)*hy(1))

   !---------------------------------------------------------------------------------
   ! Matrix for kperp^2 evaluation: k_perp^2 = k_i * gees(i,j) * k_j
   !---------------------------------------------------------------------------------
   gees = 0.0_WP

   DO i = 1, 3
      DO j = 1, 3
         DO k = 1, 3
            DO p = 1, 3
               gees(i, j) = gees(i, j) + gege(k, i)*gege(p, j)* &
                            (kron(k, p)/(lame(k)*lame(p)) - beta_vec(k)*beta_vec(p)/normB**2)
            END DO
         END DO
      END DO
   END DO

END SUBROUTINE dispersion
