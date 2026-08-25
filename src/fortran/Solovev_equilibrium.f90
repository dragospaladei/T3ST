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
       REAL(KIND=wp), DIMENSION(Np), INTENT(IN)   :: X, Z                                 ! new basis coordinates (input)
       REAL(KIND=wp), DIMENSION(Nqua, Np), INTENT(OUT)  :: R                                    ! new basis coordinates (input)

       ! This is the array-oriented counterpart of the Solovev branch in
       ! gyrocenter_drifts. Keep the two implementations mathematically synchronized.
       REAL(KIND=wp), DIMENSION(Np) :: psi, psir, psiz, psirr, psizz, psirz
       REAL(KIND=wp), DIMENSION(Np) :: chi, chir, chiz, qpsi, qprim, Fpsi, Fprim
       REAL(KIND=wp), DIMENSION(Np) :: rhot, rhotr, rhotz

       REAL(KIND=wp) :: xi, zi, xi2, zi2
       REAL(KIND=wp) :: psi_sep, delta_sep, delta2, delta_geom
       REAL(KIND=wp) :: h2, hsqrt, rpsi, rpsi_tab, rpsir, rpsiz
       REAL(KIND=wp) :: xsol, fsol, w0sol, w0sol_r, Ttor, Jedge
       REAL(KIND=wp) :: vL, vR, Cqsol, Frpsi, qrpsi
       REAL(KIND=wp) :: eta, etar, etaz, ceta, seta
       REAL(KIND=wp) :: chi_rpsi, chi_eta
       REAL(KIND=wp) :: bk, bkr, sink, cosk, snext, cnext
       INTEGER :: i, k, isol, idxL, idxR, ibL, ibR

       delta_sep = 1.0_wp - gama
       psi_sep = amp*alfa**2*delta_sep**2/(8.0_wp*(1.0_wp + alfa**2))
       Cqsol = (1.0_wp + alfa**2)/(alfa*amp)
       Jedge = Sol_data(SOL_JEDGE_IDX)

       DO i = 1, Np
          xi = X(i)
          zi = Z(i)
          xi2 = xi*xi
          zi2 = zi*zi

          psi(i) = amp*(0.25_wp*alfa**2*(xi2 - 1.0_wp)**2 + (xi2 - gama)*zi2) &
                   /(2.0_wp*(1.0_wp + alfa**2))
          psir(i) = amp*(alfa**2*xi*(xi2 - 1.0_wp) + 2.0_wp*xi*zi2) &
                    /(2.0_wp*(1.0_wp + alfa**2))
          psiz(i) = amp*(xi2 - gama)*zi/(1.0_wp + alfa**2)
          psirr(i) = amp*(2.0_wp*alfa**2*xi2 + alfa**2*(xi2 - 1.0_wp) + 2.0_wp*zi2) &
                     /(2.0_wp*(1.0_wp + alfa**2))
          psirz(i) = 2.0_wp*amp*xi*zi/(1.0_wp + alfa**2)
          psizz(i) = amp*(xi2 - gama)/(1.0_wp + alfa**2)

          Fpsi(i) = sqrt(1.0_wp + 2.0_wp*amp*gama*psi(i)/(1.0_wp + alfa**2))
          Fprim(i) = (amp*gama)/(1.0_wp + alfa**2)/Fpsi(i)

          h2 = max(xi2 - gama, eps_root)
          hsqrt = sqrt(h2)
          delta2 = (xi2 - 1.0_wp)**2 + 4.0_wp*zi2*h2/(alfa*alfa)
          delta_geom = sqrt(max(delta2, 0.0_wp))
          rpsi = delta_geom/delta_sep
          rpsi_tab = min(max(rpsi, 0.0_wp), SOL_RMAX)

          xsol = rpsi_tab*SOL_INV_DR
          isol = min(int(xsol), Nsol - 2)
          fsol = xsol - real(isol, wp)
          idxL = isol + 1
          idxR = idxL + 1

          vL = Sol_data(SOL_W0_OFF + idxL)
          vR = Sol_data(SOL_W0_OFF + idxR)
          IF (isol == 0) THEN
             w0sol = vL + (vR - vL)*fsol*fsol
             w0sol_r = 2.0_wp*(vR - vL)*fsol*SOL_INV_DR
          ELSE
             w0sol = vL + (vR - vL)*fsol
             w0sol_r = (vR - vL)*SOL_INV_DR
          END IF

          vL = Sol_data(SOL_T_OFF + idxL)
          vR = Sol_data(SOL_T_OFF + idxR)
          IF (isol == 0) THEN
             Ttor = vL + (vR - vL)*fsol*fsol
          ELSE
             Ttor = vL + (vR - vL)*fsol
          END IF

          qpsi(i) = Cqsol*Fpsi(i)*w0sol
          Frpsi = 2.0_wp*psi_sep*rpsi*Fprim(i)
          qrpsi = Cqsol*(Frpsi*w0sol + Fpsi(i)*w0sol_r)
          rhot(i) = sqrt(max(Ttor, 0.0_wp))

          IF ((rhot(i) > eps_rr) .AND. (abs(Jedge) > eps_B)) THEN
             rhotr(i) = qpsi(i)*psir(i)/(2.0_wp*psi_sep*Jedge*rhot(i))
             rhotz(i) = qpsi(i)*psiz(i)/(2.0_wp*psi_sep*Jedge*rhot(i))
          ELSE
             rhotr(i) = 0.0_wp
             rhotz(i) = 0.0_wp
          END IF

          IF ((rpsi > eps_rr) .AND. (rhot(i) > eps_rr) .AND. (abs(qpsi(i)) > eps_B)) THEN
             qprim(i) = Jedge*rhot(i)*qrpsi/(rpsi*qpsi(i))
          ELSE
             qprim(i) = 0.0_wp
          END IF

          IF (delta_geom > eps_rr) THEN
             ceta = (xi2 - 1.0_wp)/delta_geom
             seta = 2.0_wp*zi*hsqrt/(alfa*delta_geom)
             eta = atan2(seta, ceta)
             etar = -2.0_wp*xi*zi*(xi2 + 1.0_wp - 2.0_wp*gama)/(alfa*hsqrt*delta2)
             etaz = 2.0_wp*(xi2 - 1.0_wp)*hsqrt/(alfa*delta2)
          ELSE
             eta = 0.0_wp
             ceta = 1.0_wp
             seta = 0.0_wp
             etar = 0.0_wp
             etaz = 0.0_wp
          END IF

          IF (rpsi > eps_rr) THEN
             rpsir = psir(i)/(2.0_wp*psi_sep*rpsi)
             rpsiz = psiz(i)/(2.0_wp*psi_sep*rpsi)
          ELSE
             rpsir = 0.0_wp
             rpsiz = 0.0_wp
          END IF

          chi(i) = eta
          chi_rpsi = 0.0_wp
          chi_eta = 1.0_wp
          sink = seta
          cosk = ceta

          DO k = 1, Ksol
             ibL = SOL_B_OFF + (k - 1)*Nsol + idxL
             ibR = ibL + 1
             vL = Sol_data(ibL)
             vR = Sol_data(ibR)

             IF (isol == 0) THEN
                bk = vR*fsol**k
                bkr = real(k, wp)*vR*fsol**(k - 1)*SOL_INV_DR
             ELSE
                bk = vL + (vR - vL)*fsol
                bkr = (vR - vL)*SOL_INV_DR
             END IF

             chi(i) = chi(i) + bk*sink/real(k, wp)
             chi_rpsi = chi_rpsi + bkr*sink/real(k, wp)
             chi_eta = chi_eta + bk*cosk

             IF (k < Ksol) THEN
                snext = sink*ceta + cosk*seta
                cnext = cosk*ceta - sink*seta
                sink = snext
                cosk = cnext
             END IF
          END DO

          chir(i) = chi_rpsi*rpsir + chi_eta*etar
          chiz(i) = chi_rpsi*rpsiz + chi_eta*etaz
       END DO

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
       REAL(KIND=wp), DIMENSION(Np), INTENT(OUT)  :: X, Y                                    ! new basis coordinates (input)

       ! Local variables
       REAL(KIND=wp), DIMENSION(Np)     :: aux                                    ! new basis coordinates (input)
       REAL(KIND=wp), DIMENSION(NgridR) :: F1                                    ! new basis coordinates (input)
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
