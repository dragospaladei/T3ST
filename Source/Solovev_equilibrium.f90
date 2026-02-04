    !      ! Subroutine :: Solovev
    !      ! Purpose    :: Solovev equilibrium
    !      ! F(psi) = 1 + ceva*psi
    !      ! psi ~ Z^2*(R^2-gama*R0^2)+alfa^2/4*(R^2-R0^2)^2
    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  01/01/2024   |     D. I. Palade       |   Remake
    !  ========================================================================================================================================================

    SUBROUTINE Solovev(X, Z, R)
       USE constants
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN)   :: X, Z                                 ! new basis coordinates (input)
       REAL(KIND=dp), DIMENSION(Nqua, Np), INTENT(OUT)  :: R                                    ! new basis coordinates (input)

       ! Local variables
       ! Local variables
       REAL(KIND=dp), DIMENSION(Np)                       :: psi, psir, psiz, psirr, psizz, psirz  ! miscelaneous
       REAL(KIND=dp), DIMENSION(Np)                       :: chi, chir, chiz, qpsi, qprim, Fpsi, Fprim! miscelaneous
       REAL(KIND=dp), DIMENSION(Np)                       :: rr, theta, rhot, rhotr, rhotz, psiradius

       ! \psi    :: the poloidal flux
       psi = (amp*((alfa**2*(-1.0 + X**2)**2)/4.+(-gama + X**2)*Z**2))/(2.*(1 + alfa**2))
       ! d_R\psi :: R derivative of the poloidal flux (psi)
       psir = (amp*(alfa**2*X*(-1.0 + X**2) + 2.0*X*Z**2))/(2.*(1.0 + alfa**2))
       ! d_Z\psi :: Z derivative of the poloidal flux (psi)
       psiz = (amp*(-gama + X**2)*Z)/(1.0 + alfa**2)
       ! d_RR\psi :: R,R derivative of the poloidal flux (psi)
       psirr = (amp*(2.0*alfa**2*X**2 + alfa**2*(-1.0 + X**2) + 2.0*Z**2))/(2.*(1.0 + alfa**2))
       ! d_RZ\psi :: R,Z derivative of the poloidal flux (psi)
       psirz = (2.0*amp*X*Z)/(1.0 + alfa**2)
       ! d_ZZ\psi :: Z,Z derivative of the poloidal flux (psi)
       psizz = (amp*(-gama + X**2))/(1.0 + alfa**2)
       ! F(\psi)

       Fpsi = sqrt(1.0 + 2.0*amp*gama*psi/(1.0 + alfa**2))
       Fprim = (amp*gama)/(1.0 + alfa**2)*Fpsi

    !!! taken from circular ...
       ! straight poloidal angle = arctan(theta...)

       psiradius = (amp*((alfa**2*(-1.0 + (1.0 + a0)**2)**2)/4.))/(2.*(1 + alfa**2))
       rr = sqrt(psi/psiradius)
       theta = atan2(Z, X - 1.0)      ! atan(z,rr)
       chi = 2.0*atan(sqrt((1.0 - rr)/(1.0 + rr))*tan(theta/2.0))
       chir = sin(theta)*(rr**2 - X)/X/rr/Sqrt(1.0 - rr**2)
       chiz = -(rr - cos(theta))/rr/Sqrt(1.0 - rr**2)
       rhot = rr
       rhotr = rr/2.0*psir/psi
       rhotz = rr/2.0*psiz/psi

       ! q(r)    :: the safety factor
       qpsi = safe1 + safe2*rr*(1.0/a0) + safe3*rr**2*(1.0/a0)**2
       ! q'(r)   :: radial derivative of the safety factor
       qprim = safe2*(1.0/a0) + 2.0*safe3*rr*(1.0/a0)**2 ! a0 a fost deja scalat

       !       R = transpose(reshape([psi, psir, psiz, psirr, psirz, psizz, Fpsi, Fprim, qpsi, qprim, rhot, rhotr, rhotz, chi, chir, chiz], [Np,Nqua]))

   R(1, :) = psi
   R(2, :) = psir
   R(3, :) = psiz
   R(4, :) = psirr
   R(5, :) = psirz
   R(6, :) = psizz
   R(7, :) = Fpsi
   R(8, :) = Fprim
   R(9, :) = qpsi
   R(10, :) = qprim
   R(11, :) = rhot
   R(12, :) = rhotr
   R(13, :) = rhotz
   R(14, :) = chi
   R(15, :) = chir
   R(16, :) = chiz

    END SUBROUTINE Solovev

    SUBROUTINE psisurf_solov2(X, Y)
       USE constants
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT)  :: X, Y                                    ! new basis coordinates (input)

       ! Local variables
       REAL(KIND=dp), DIMENSION(Np)     :: aux                                    ! new basis coordinates (input)
       REAL(KIND=dp), DIMENSION(NgridR) :: F1                                    ! new basis coordinates (input)
       REAl, DIMENSION(2000) :: aux1, aux2, aux3, psi                                   ! new basis coordinates (input)
       INTEGER :: k, ioc

       CALL random_number(aux)
       CALL random_number(aux3)
       aux = pi*(2.0*aux - 1.0)            ! theta
       aux3 = a0*aux3

       do k = 1, Np
          aux2 = aux3*Sin(aux(k))
          aux1 = aux3*Cos(aux(k))

          psi = (amp*((alfa**2*(-1.0 + (aux1 + 1.0)**2)**2)/4.+(-gama + (aux1 + 1.0)**2)*aux2**2))/(2.*(1 + alfa**2))

          ioc = Minloc(abs(psi - psi0), 1)
          X(k) = 1.0 + aux1(ioc)
          Y(k) = aux2(ioc)

          !        if(sqrt((X(k)-1.0)**2 + Y(k)**2).gt.a0) then
          !            X(k) = X0
          !            Y(k) = Y0
          !        else
          !        endif

       end do
    END SUBROUTINE psisurf_solov2

