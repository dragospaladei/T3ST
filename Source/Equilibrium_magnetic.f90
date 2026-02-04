    !      ! Subroutine :: coordinates(q1,q2,q3,Qr,DQr,u1,u2,u3,hh3,Alfa,Beta)
    !      ! Purpose    :: Computes old (toroidal) coordinates (u1,u2,u3) from new (field-aligned) positions (q1,q2,q3)==(x,y,z)
    !                      Computes co/contravariant basis vectors (Alfa, Beta) of new basis (q1,q2,q3)
    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  01/01/2024   |     D. I. Palade       |   Remake
    !  ========================================================================================================================================================

    SUBROUTINE Equil(X, Y, Z, Bx, By, Bz, gradBx, gradBy, gradBz, gradu2x, gradu2y, gradu2z,&
        &rotbx, rotby, rotbz, rotux, rotuy, rotuz, phi0x, phi0y, phi0z, B, G, q1, q2, q3)
       USE constants
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN)   :: X, Y, Z                                             ! new basis coordinates (input)
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT)  :: phi0x, phi0y, phi0z                                 ! The gradient of the neoclassical field Phi_0 (for quasineutrality in rotating plasmas)
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT)  :: gradBx, gradBy, gradBz, gradu2x, gradu2y, gradu2z   ! Grad B, grad u**2
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT)  :: rotbx, rotby, rotbz, rotux, rotuy, rotuz            ! curl(b), curl(u) components in the x,y,z system
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT)  :: Bx, By, Bz, B                                       ! Bx,  By, Bz, B components in the x,y,z system
       REAL(KIND=dp), DIMENSION(3, 3, Np), INTENT(OUT)  :: G                                                   ! G_ij = grad(qj).dr/dxi;
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT)  :: q1, q2, q3                                            ! Scaled field-aligned coordinates, a la GENE

       ! Local variables
       REAL(KIND=dp), DIMENSION(Np)                    :: hx, hy, hz                                          ! lame coefficients
       REAL(KIND=dp), DIMENSION(Np)                    :: psi, psir, psiz, psirr, psizz, psirz, rhot, rhotr, rhotz                ! miscelaneous
       REAL(KIND=dp), DIMENSION(Np)                    :: chi, chir, chiz, qpsi, qprim, Fpsi, Fprim              ! miscelaneous
       REAL(KIND=dp), DIMENSION(Nqua, Np)             :: R                                                   ! G_ij = grad(qj).dr/dxi;

       hx = 1.0                                                      ! lame coefficient for x = R -> hx = 1
       hy = 1.0                                                      ! lame coefficient for y = Z -> hy = 1
       hz = X                                                        ! lame coefficient for z = zeta -> hz = R = X

       ! ------------------------------------------------------------------------------------------------------------------------
       ! R contains: d_r\psi, d_z\psi, F(psi), rot(b), grad(B)
       ! ------------------------------------------------------------------------------------------------------------------------

       CALL MagnModel(X, Y, R)

       psi = R(1, :)
       psir = R(2, :)
       psiz = R(3, :)
       psirr = R(4, :)
       psirz = R(5, :)
       psizz = R(6, :)
       Fpsi = R(7, :)
       Fprim = R(8, :)
       qpsi = R(9, :)
       qprim = R(10, :)
       rhot = R(11, :)
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
       q1 = C1*rhot
       q2 = C2*(Z - qpsi*chi)
       q3 = C3*chi

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

       B = sqrt(Bx**2*hx**2 + By**2*hy**2 + Bz**2*hz**2)             ! the norm of the magnetic field (with the covariant components)
       ! ------------------------------------------------------------------------------------------------------------------------
       ! Gradients and curls
       ! ------------------------------------------------------------------------------------------------------------------------

       gradBx = -(Fpsi**2 + psir**2 + psiz**2 - Fpsi*Fprim*X*psir - X*psir*psirr - psiz*X*psirz)/B/X**3   ! grad B proj. on x
       gradBy = (Fpsi*Fprim*psiz + psir*psirz + psiz*psizz)/B/X**2   ! grad B proj. on y
       gradBz = 0.0                                                  ! grad B proj. on z

        rotbx = (Fprim*psiz*(psir**2+psiz**2) - Fpsi*(psir*psirz + psiz*psizz))/B**3/X**3                                                     ! curl of u proj on x
        rotby = - (Fpsi**3 + Fpsi*(psir**2+psiz**2) - X*Fpsi*(psirz*psiz+psir*psirr) + X*Fprim*psir*(psir**2+psiz**2))/B**3/X**4              ! curl of u proj on y
       !      rotbz = - (Fpsi*Fprim*(psir**2+psiz**2) - Fpsi**2*(psir**2+psiz**2) - &
       !      &psirr*psiz*psiz - psizz*psir*psir + 2.0*psirz*psir*psiz)/B**3/X**4  ! curl of u proj on z
       rotbz = -(Fpsi*Fprim*(psir**2 + psiz**2) - Fpsi**2*(psirr + psizz) - &
           &psirr*psiz*psiz - psizz*psir*psir + 2.0*psirz*psir*psiz)/B**3/X**4  ! curl of u proj on z

       rotux = Omgt0*(-X*psiz*Omgtprim)                              ! curl of u proj on x
       rotuy = Omgt0*(2.0 + 2.0*Omgtprim*(psi - psi0) + X*psir*Omgtprim)                         ! curl of u proj on y
       rotuz = 0.0                                                   ! curl of u proj on z

       gradu2z = 1.0 + Omgtprim*(psi - psi0)
       gradu2x = Omgt0**2*X**2*(psir*Omgtprim*gradu2z + gradu2z**2/X)! grad u**2 proj. on x
       gradu2y = Omgt0**2*X**2*(psiz*Omgtprim*gradu2z)               ! grad u**2 proj. on y
       gradu2z = 0.0                                                 ! grad u**2 proj. on z

       phi0z = 1.0 + Omgtprim*(psi - psi0)
       phi0x = Te/(Ti + Te)*(gradu2x + Omgt0**2*phi0z**2)              ! grad Phi0 proj. on x
       phi0y = Te/(Ti + Te)*gradu2y*2.0*(X - 1.0)                        ! grad Phi0 proj. on y
       phi0z = 0.0                                                   ! grad Phi0 proj. on z

    END SUBROUTINE Equil

    !  ========================================================================================================================================================
