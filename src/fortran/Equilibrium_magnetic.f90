    !      ! Subroutine :: coordinates(q1,q2,q3,Qr,DQr,u1,u2,u3,hh3,Alfa,Beta)
    !      ! Purpose    :: Computes old (toroidal) coordinates (u1,u2,u3) from new (field-aligned) positions (q1,q2,q3)==(x,y,z)
    !                      Computes co/contravariant basis vectors (Alfa, Beta) of new basis (q1,q2,q3)
    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  01/01/2024   |     D. I. Palade       |   Remake
    !  ========================================================================================================================================================

    SUBROUTINE Equil(X, Y, Z, Bx, By, Bz, B, G)
       ! Former unused outputs:
       ! gradBx, gradBy, gradBz, gradu2x, gradu2y, gradu2z,
       ! rotbx, rotby, rotbz, rotux, rotuy, rotuz,
       ! phi0x, phi0y, phi0z, q1, q2, q3
       USE constants
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=wp), DIMENSION(Np), INTENT(IN) :: X, Y, Z
!      REAL(KIND=wp), DIMENSION(Np), INTENT(OUT)  :: phi0x, phi0y, phi0z
!      REAL(KIND=wp), DIMENSION(Np), INTENT(OUT)  :: gradBx, gradBy, gradBz, gradu2x, gradu2y, gradu2z
!      REAL(KIND=wp), DIMENSION(Np), INTENT(OUT)  :: rotbx, rotby, rotbz, rotux, rotuy, rotuz
       REAL(KIND=wp), DIMENSION(Np), INTENT(OUT) :: Bx, By, Bz, B
       REAL(KIND=wp), DIMENSION(3, 3, Np), INTENT(OUT) :: G
!      REAL(KIND=wp), DIMENSION(Np), INTENT(OUT)  :: q1, q2, q3

       ! Local variables
       REAL(KIND=wp), DIMENSION(Np) :: hx, hy, hz                                          ! lame coefficients
       REAL(KIND=wp), DIMENSION(Np) :: psir, psiz, rhotr, rhotz
       REAL(KIND=wp), DIMENSION(Np) :: chi, chir, chiz, qpsi, qprim, Fpsi
!      REAL(KIND=wp), DIMENSION(Np)                    :: psi, psirr, psizz, psirz, rhot, Fprim
       REAL(KIND=wp), DIMENSION(Nqua, Np) :: R                                                   ! G_ij = grad(qj).dr/dxi;

       hx = 1.0                                                      ! lame coefficient for x = R -> hx = 1
       hy = 1.0                                                      ! lame coefficient for y = Z -> hy = 1
       hz = X                                                        ! lame coefficient for z = zeta -> hz = R = X

       ! ------------------------------------------------------------------------------------------------------------------------
       ! R contains: d_r\psi, d_z\psi, F(psi), rot(b), grad(B)
       ! ------------------------------------------------------------------------------------------------------------------------

       CALL equilibrium_all(Np, X, Y, R)

!      psi = R(1, :)
       psir = R(2, :)
       psiz = R(3, :)
!      psirr = R(4, :)
!      psirz = R(5, :)
!      psizz = R(6, :)
       Fpsi = R(7, :)
!      Fprim = R(8, :)
       qpsi = R(9, :)
       qprim = R(10, :)
!      rhot = R(11, :)
       rhotr = R(12, :)
       rhotz = R(13, :)
       chi = R(14, :)
       chir = R(15, :)
       chiz = R(16, :)

       ! ------------------------------------------------------------------------------------------------------------------------
       ! Scaled field-aligned coordinates, a la GENE
       ! ------------------------------------------------------------------------------------------------------------------------
       !Atentie : pentru a evita jump - ul de la \[Chi] =  2 \[Pi], adica lipsa periodicitatii poloidale,
       !am ales sa definesc y = Cy (\[Xi] -  q0*\[Chi])!!! aceasta este echivalent cu un safety factor \
       !constant; sau magnetic shear null
       !       q2 = C2*(Z - qpsi*chi)
       ! instead of
       !      G(1,2,:) = - C2*(qprim*chi*psir + qpsi*chir)
       !      G(2,2,:) = - C2*(qprim*chi*psiz + qpsi*chiz)
       !   Ba nu, am revenit la qpsi; pentru ca altfel nu se alinia corect driftul neoclasic cu directia y
       !       q1 = C1*(psi - psi0)
!      q1 = C1*rhot
!      q2 = C2*(Z - qpsi*chi)
!      q3 = C3*chi

       ! ------------------------------------------------------------------------------------------------------------------------
       ! The contra-variant elements of grad(qj) :: G_ij = grad(q_j).dr/dxi;
       ! ------------------------------------------------------------------------------------------------------------------------

       !            G(1,1,:) = C1*psir                                            ! d_R\psi
       !     G(2,1,:) = C1*psiz
       G(1, 1, :) = C1*rhotr                                            ! d_R\psi
       G(2, 1, :) = C1*rhotz
       G(3, 1, :) = 0.0                                                ! d_Z\psi

       G(1, 3, :) = C3*chir
       G(2, 3, :) = C3*chiz
       G(3, 3, :) = 0.0

       !    A correction must be done:  since q2 = C2*(Z - q00*chi), then
       G(1, 2, :) = -C2*qpsi*chir - C2*chi*qprim*rhotr
       G(2, 2, :) = -C2*qpsi*chiz - C2*chi*qprim*rhotz
       G(3, 2, :) = C2

       ! ------------------------------------------------------------------------------------------------------------------------
       ! Note that we assume that the toroidal angular frequency is constant ...
       ! also, <R^2> is independent of R,Z: that is strange...
       ! ------------------------------------------------------------------------------------------------------------------------

       ! ------------------------------------------------------------------------------------------------------------------------
       ! Magnetic field components
       ! ------------------------------------------------------------------------------------------------------------------------

       Bx = -psiz/X                                                ! covariant component of B on x  (d_z\psi/R)
       By = psir/X                                                ! covariant component of B on y  (d_rpsi/R)
       Bz = Fpsi/X**2                                             ! covariant component of B on z  (F(psi)/R^2)

       B = sqrt(Bx**2*hx**2 + By**2*hy**2 + Bz**2*hz**2)
       ! ------------------------------------------------------------------------------------------------------------------------
       ! Gradients and curls
       ! ------------------------------------------------------------------------------------------------------------------------

!      gradBx = -(Fpsi**2 + psir**2 + psiz**2 - Fpsi*Fprim*X*psir - X*psir*psirr - psiz*X*psirz)/B/X**3
!      gradBy = (Fpsi*Fprim*psiz + psir*psirz + psiz*psizz)/B/X**2
!      gradBz = 0.0

!      rotbx = (Fprim*psiz*(psir**2+psiz**2) - Fpsi*(psir*psirz + psiz*psizz))/B**3/X**3
!      rotby = -(Fpsi**3 + Fpsi*(psir**2+psiz**2) - X*Fpsi*(psirz*psiz+psir*psirr) &
!              + X*Fprim*psir*(psir**2+psiz**2))/B**3/X**4
       !      rotbz = - (Fpsi*Fprim*(psir**2+psiz**2) - Fpsi**2*(psir**2+psiz**2) - &
       !      &psirr*psiz*psiz - psizz*psir*psir + 2.0*psirz*psir*psiz)/B**3/X**4  ! curl of u proj on z
!      rotbz = -(Fpsi*Fprim*(psir**2 + psiz**2) - Fpsi**2*(psirr + psizz) - &
!               psirr*psiz*psiz - psizz*psir*psir + 2.0*psirz*psir*psiz)/B**3/X**4

!      rotux = Omgt0*(-X*psiz*Omgtprim)
!      rotuy = Omgt0*(2.0 + 2.0*Omgtprim*(psi - psi0) + X*psir*Omgtprim)
!      rotuz = 0.0

!      gradu2z = 1.0 + Omgtprim*(psi - psi0)
!      gradu2x = Omgt0**2*X**2*(psir*Omgtprim*gradu2z + gradu2z**2/X)
!      gradu2y = Omgt0**2*X**2*(psiz*Omgtprim*gradu2z)
!      gradu2z = 0.0

!      phi0z = 1.0 + Omgtprim*(psi - psi0)
!      phi0x = Te/(Ti + Te)*(gradu2x + Omgt0**2*phi0z**2)
!      phi0y = Te/(Ti + Te)*gradu2y*2.0*(X - 1.0)
!      phi0z = 0.0

    END SUBROUTINE Equil

!========================================================================================================
! Model dispatch and reduced magnetic-equilibrium queries
!========================================================================================================
    !  ========================================================================================================================================================
    !      ! File       :: Equilibrium_magnetic.f90
    !      ! Purpose    :: Dispatch magnetic-equilibrium calls to EFIT, Solovev, or circular models
    !      !            :: according to magnetic_model.
    !
    !      ! Subroutine :: equilibrium_all
    !      ! Purpose    :: Generate the poloidal flux function, F(psi), rot(b) and grad(B)
    !      !            :: Different choices are possible, indexed
    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  01/01/2024   |     D. I. Palade       |   Initial version
    !  ========================================================================================================================================================

    SUBROUTINE equilibrium_all(npoints, X, Y, R)
       USE constants
       IMPLICIT NONE

       ! I/O variables
       INTEGER, INTENT(IN) :: npoints
       REAL(KIND=wp), DIMENSION(npoints), INTENT(IN) :: X, Y
       REAL(KIND=wp), DIMENSION(Nqua, npoints), INTENT(OUT) :: R

       IF (magnetic_model .EQ. 1) THEN
          CALL EFIT2(npoints, X, Y, R)                                      ! realistic equilibrium:: interpolation from Efit files
       ELSEIF (magnetic_model .EQ. 2) THEN
          CALL Solovev(npoints, X, Y, R)                                    ! Solovev equilibrium ....
       ELSEIF (magnetic_model .EQ. 3) THEN
          CALL Circular(npoints, X, Y, R)                                   ! Circular equilibria ....
       ELSE
          WRITE (*, *) 'Error: The index for the magnetic model does not correspond to any choice'
       END IF

    END SUBROUTINE equilibrium_all

    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  23/08/2026   |     chatgpt            |   Initial version
    !  ========================================================================================================================================================
    SUBROUTINE equilibrium_q_rhot_psi(R, Z, qval, rhot, psival)
       USE constants
       IMPLICIT NONE

       REAL(KIND=wp), INTENT(IN) :: R, Z
       REAL(KIND=wp), INTENT(OUT) :: qval, rhot, psival

       REAL(KIND=wp) :: position_R(1), position_Z(1)
       REAL(KIND=wp) :: equilibrium(Nqua, 1)

       position_R(1) = R
       position_Z(1) = Z

       CALL equilibrium_all(1, position_R, position_Z, equilibrium)

       psival = equilibrium(1, 1)
       qval = equilibrium(9, 1)
       rhot = equilibrium(11, 1)

    END SUBROUTINE equilibrium_q_rhot_psi

!========================================================================================================
! Evaluate the magnetic-field norm and toroidal-flux radius needed during
! particle initialization. The full model is evaluated once, keeping B and
! rhot mathematically consistent while avoiding two separate equilibrium passes.
!========================================================================================================
    SUBROUTINE equilibrium_initialization_fields(X, Y, B, rhot)
       USE constants
       IMPLICIT NONE

       REAL(KIND=wp), DIMENSION(Np), INTENT(IN) :: X, Y
       REAL(KIND=wp), DIMENSION(Np), INTENT(OUT) :: B, rhot

       REAL(KIND=wp), DIMENSION(Nqua, Np) :: equilibrium

       CALL equilibrium_all(Np, X, Y, equilibrium)

       rhot = equilibrium(11, :)
       B = SQRT(equilibrium(7, :)**2 + equilibrium(2, :)**2 + equilibrium(3, :)**2)/X

    END SUBROUTINE equilibrium_initialization_fields

    !  ========================================================================================================================================================
