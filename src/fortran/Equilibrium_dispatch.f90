    !  ========================================================================================================================================================
    !      ! File       :: Equilibrium_dispatch.f90
    !      ! Purpose    :: Dispatch magnetic-equilibrium calls to EFIT, Solovev, or circular models
    !      !            :: according to magnetic_model.
    !
    !      ! Subroutine :: MagnModel
    !      ! Purpose    :: Generate the poloidal flux function, F(psi), rot(b) and grad(B)
    !      !            :: Different choices are possible, indexed
    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  01/01/2024   |     D. I. Palade       |   Initial version
    !  ========================================================================================================================================================

    SUBROUTINE MagnModel(X, Y, R)
       USE constants
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=wp), DIMENSION(Np), INTENT(IN)   :: X, Y                  ! new basis coordinates (input)
       REAL(KIND=wp), DIMENSION(Nqua, Np), INTENT(OUT)  :: R                     ! new basis coordinates (input)

       IF (magnetic_model .eq. 1) THEN
          CALL EFIT2(X, Y, R)                                               ! realistic equilibrium:: interpolation from Efit files
       ELSEIF (magnetic_model .eq. 2) THEN
          CALL Solovev(X, Y, R)                                             ! Solovev equilibrium ....
       ELSEIF (magnetic_model .eq. 3) THEN
          CALL Circular(X, Y, R)                                            ! Circular equilibria ....
       ELSE
          WRITE (*, *) 'Error: The index for the magnetic model does not correspond to any choice'
       END IF

    END SUBROUTINE MagnModel

    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  23/08/2026   |     chatgpt            |   Initial version
    !  ========================================================================================================================================================
SUBROUTINE equilibrium_q_rhot_psi(R, Z, qval, rhot,psival)

   USE constants
   IMPLICIT NONE

   REAL(KIND=wp), INTENT(IN)  :: R, Z
   REAL(KIND=wp), INTENT(OUT) :: qval, rhot, psival

   ! locals
   REAL(KIND=wp) :: rr
   REAL(KIND=wp) :: X1, Y1, Xef, Yef
   REAL(KIND=wp) :: F1, F2, F3, F4
   REAL(KIND=wp) :: efit_q, efit_rhot, efit_psi
   INTEGER :: poz1, poz2

   ! Solovev locals
   REAL(KIND=wp) :: delta_sep, delta2, delta_geom
   REAL(KIND=wp) :: rpsi, rpsi_tab
   REAL(KIND=wp) :: xsol, fsol
   REAL(KIND=wp) :: w0, Ttor, Fpsi
   REAL(KIND=wp) :: vL, vR
   INTEGER :: isol, idxL, idxR

   SELECT CASE (magnetic_model)

   !===============================================================
   ! EFIT
   !===============================================================
   CASE (1)

      poz1 = modulo(int((R - minR)/stepR), NgridR)
      poz2 = modulo(int((Z - minZ)/stepZ), NgridZ)

      X1 = minR + poz1*stepR
      Y1 = minZ + poz2*stepZ

      Xef = (R - X1)/stepR
      Yef = (Z - Y1)/stepZ

      !------------------------------------------------------------
      ! psi = EFIT quantity 1
      !------------------------------------------------------------

      F1 = Efit_data(0*Ngrid + poz1*NgridR + poz2)
      F2 = Efit_data(0*Ngrid + (poz1+1)*NgridR + poz2)
      F3 = Efit_data(0*Ngrid + poz1*NgridR + (poz2+1))
      F4 = Efit_data(0*Ngrid + (poz1+1)*NgridR + (poz2+1))

      psival = F1 + (F2-F1)*Xef + (F3-F1)*Yef &
             + Xef*Yef*(F1 + F4 - F2 - F3)

      !------------------------------------------------------------
      ! q = EFIT quantity 9
      !------------------------------------------------------------

      F1 = Efit_data(8*Ngrid + poz1*NgridR + poz2)
      F2 = Efit_data(8*Ngrid + (poz1+1)*NgridR + poz2)
      F3 = Efit_data(8*Ngrid + poz1*NgridR + (poz2+1))
      F4 = Efit_data(8*Ngrid + (poz1+1)*NgridR + (poz2+1))

      qval = F1 + (F2-F1)*Xef + (F3-F1)*Yef &
           + Xef*Yef*(F1 + F4 - F2 - F3)

      !------------------------------------------------------------
      ! rhot = EFIT quantity 11
      !------------------------------------------------------------

      F1 = Efit_data(10*Ngrid + poz1*NgridR + poz2)
      F2 = Efit_data(10*Ngrid + (poz1+1)*NgridR + poz2)
      F3 = Efit_data(10*Ngrid + poz1*NgridR + (poz2+1))
      F4 = Efit_data(10*Ngrid + (poz1+1)*NgridR + (poz2+1))

      rhot = F1 + (F2-F1)*Xef + (F3-F1)*Yef &
           + Xef*Yef*(F1 + F4 - F2 - F3)


   !===============================================================
   ! Solovev
   !===============================================================
   CASE (2)

      ! exact psi
      psival = amp*( &
           0.25_wp*alfa**2*(R**2 - 1.0_wp)**2 &
         + (R**2 - gama)*Z**2 ) &
         /(2.0_wp*(1.0_wp + alfa**2))

      ! physical sqrt-poloidal-flux label
      delta_sep = 1.0_wp - gama

      delta2 = (R**2 - 1.0_wp)**2 &
             + 4.0_wp*Z**2*(R**2 - gama)/(alfa*alfa)

      delta_geom = sqrt(max(delta2,0.0_wp))

      rpsi = delta_geom/delta_sep

      ! table domain
      rpsi_tab = min(max(rpsi,0.0_wp),SOL_RMAX)

      xsol = rpsi_tab*SOL_INV_DR
      isol = min(int(xsol),Nsol-2)

      fsol = xsol - real(isol,wp)

      idxL = isol + 1
      idxR = idxL + 1

      !------------------------------------------------------------
      ! w0
      !------------------------------------------------------------

      vL = Sol_data(SOL_W0_OFF + idxL)
      vR = Sol_data(SOL_W0_OFF + idxR)

      if (isol == 0) then
         w0 = vL + (vR-vL)*fsol*fsol
      else
         w0 = vL + (vR-vL)*fsol
      end if

      ! q
      Fpsi = sqrt(1.0_wp &
           + 2.0_wp*amp*gama*psival/(1.0_wp + alfa**2))

      qval = (1.0_wp + alfa**2)/(alfa*amp) * Fpsi*w0

      !------------------------------------------------------------
      ! rhot^2 = T
      !------------------------------------------------------------

      vL = Sol_data(SOL_T_OFF + idxL)
      vR = Sol_data(SOL_T_OFF + idxR)

      if (isol == 0) then
         Ttor = vL + (vR-vL)*fsol*fsol
      else
         Ttor = vL + (vR-vL)*fsol
      end if

      rhot = sqrt(max(Ttor,0.0_wp))


   !===============================================================
   ! Circular
   !===============================================================
   CASE (3)

      rr = sqrt((R - 1.0_wp)**2 + Z**2)

      rhot = rr/a0

      qval = (s3*rhot + s2)*rhot + s1
      psival = a0**2/s2*(rr/a0 - s1/s2*log(1.0_wp + s2/s1*rr/a0))

   CASE DEFAULT

      qval   = 0.0_wp
      rhot   = 0.0_wp
      psival = 0.0_wp

   END SELECT

END SUBROUTINE equilibrium_q_rhot_psi

   !====================================================================
   !====================================================================
   !====================================================================
SUBROUTINE equilibrium_rhot_all(X7, Y7, Z7, q17)

   USE constants
   IMPLICIT NONE

   !--------------------------------------------------------------------
   ! I/O
   !--------------------------------------------------------------------
   REAL(KIND=wp), DIMENSION(Np), INTENT(IN)  :: X7, Y7, Z7
   REAL(KIND=wp), DIMENSION(Np), INTENT(OUT) :: q17

   !--------------------------------------------------------------------
   ! General locals
   !--------------------------------------------------------------------
   REAL(KIND=wp), DIMENSION(Np) :: rr

   REAL(KIND=wp) :: xi, yi
   REAL(KIND=wp) :: X1, Y1, Xef, Yef
   REAL(KIND=wp) :: F1, F2, F3, F4, efit_vals

   INTEGER :: i, poz1, poz2, jax

   !--------------------------------------------------------------------
   ! Solovev locals
   !--------------------------------------------------------------------
   REAL(KIND=wp) :: delta_sep
   REAL(KIND=wp) :: delta2, delta_geom
   REAL(KIND=wp) :: rpsi, rpsi_tab
   REAL(KIND=wp) :: xsol, fsol
   REAL(KIND=wp) :: TL, TR, Ttor

   INTEGER :: isol, idxL, idxR


   !====================================================================
   ! Circular equilibrium
   !====================================================================

   IF (magnetic_model == 3) THEN

      rr = sqrt((X7 - 1.0_wp)**2 + Y7**2)

      q17 = rr/a0


   !====================================================================
   ! Solovev equilibrium
   !
   ! q17 = rho_t = sqrt[T(rpsi)]
   !
   ! where
   !
   !     rpsi = sqrt(psi/psi_sep)
   !          = delta/(1-gamma)
   !
   ! and T(rpsi)=rho_t^2 is stored in Sol_data.
   !====================================================================

   ELSEIF (magnetic_model == 2) THEN

      delta_sep = 1.0_wp - gama

      DO i = 1, Np

         xi = X7(i)
         yi = Y7(i)

         !---------------------------------------------------------------
         ! Exact Solovev flux-surface label
         !
         ! delta^2 =
         !
         !   (R^2-1)^2
         ! + 4 Z^2 (R^2-gamma)/alpha^2
         !---------------------------------------------------------------

         delta2 = (xi*xi - 1.0_wp)**2 &
                + 4.0_wp*yi*yi*(xi*xi - gama)/(alfa*alfa)

         delta_geom = sqrt(max(delta2, 0.0_wp))

         rpsi = delta_geom/delta_sep


         !---------------------------------------------------------------
         ! Sol_data exists only up to SOL_RMAX.
         !
         ! This clipping should normally do nothing: particles should
         ! remain inside the tabulated Solovev domain.
         !---------------------------------------------------------------

         rpsi_tab = min(max(rpsi, 0.0_wp), SOL_RMAX)


         !---------------------------------------------------------------
         ! Uniform radial table lookup
         !
         ! SOL_INV_DR = (Nsol-1)/SOL_RMAX
         !---------------------------------------------------------------

         xsol = rpsi_tab*SOL_INV_DR

         isol = min(int(xsol), Nsol - 2)

         fsol = xsol - real(isol, wp)

         idxL = isol + 1
         idxR = idxL + 1


         !---------------------------------------------------------------
         ! Interpolate
         !
         !       T(rpsi) = rho_t^2
         !
         ! Near the magnetic axis T ~ rpsi^2, hence the quadratic
         ! interpolation in the first cell, matching gyrocenter_drifts.
         !---------------------------------------------------------------

         TL = Sol_data(SOL_T_OFF + idxL)
         TR = Sol_data(SOL_T_OFF + idxR)

         IF (isol == 0) THEN

            Ttor = TL + (TR - TL)*fsol*fsol

         ELSE

            Ttor = TL + (TR - TL)*fsol

         END IF


         q17(i) = sqrt(max(Ttor, 0.0_wp))

      END DO


   !====================================================================
   ! EFIT equilibria
   !====================================================================

   ELSEIF (magnetic_model == 1) THEN

      DO i = 1, Np

         xi = X7(i)
         yi = Y7(i)

         ! Grid indices
         poz1 = modulo(int((xi - minR)/stepR), NgridR)
         poz2 = modulo(int((yi - minZ)/stepZ), NgridZ)

         X1 = minR + poz1*stepR
         Y1 = minZ + poz2*stepZ

         Xef = (xi - X1)/stepR
         Yef = (yi - Y1)/stepZ

         !---------------------------------------------------------------
         ! EFIT quantity 11 = rho_t
         !---------------------------------------------------------------

         jax = 11

         F1 = Efit_data(NgridR*NgridZ*(jax - 1) + &
                        poz1*NgridR + poz2)

         F2 = Efit_data(NgridR*NgridZ*(jax - 1) + &
                        (poz1 + 1)*NgridR + poz2)

         F3 = Efit_data(NgridR*NgridZ*(jax - 1) + &
                        poz1*NgridR + (poz2 + 1))

         F4 = Efit_data(NgridR*NgridZ*(jax - 1) + &
                        (poz1 + 1)*NgridR + (poz2 + 1))

         efit_vals = F1 &
                   + (F2 - F1)*Xef &
                   + (F3 - F1)*Yef &
                   + Xef*Yef*(F1 + F4 - F2 - F3)

         q17(i) = efit_vals

      END DO


   !====================================================================
   ! Unknown model
   !====================================================================

   ELSE

      q17 = 0.0_wp

   END IF

END SUBROUTINE equilibrium_rhot_all
    !      ! Subroutine :: coordinates(q1,q2,q3,Qr,DQr,u1,u2,u3,hh3,Alfa,Beta)
    !      ! Purpose    :: Computes old (toroidal) coordinates (u1,u2,u3) from new (field-aligned) positions (q1,q2,q3)==(x,y,z)
    !                      Computes co/contravariant basis vectors (Alfa, Beta) of new basis (q1,q2,q3)
    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  01/01/2024   |     D. I. Palade       |   Remake
    !  ========================================================================================================================================================

    SUBROUTINE magnetic_norm(X, Y, Z, B)                                                                  ! F_ji = (dxj x b).dxi ;; matrix from the ExB term
       USE constants
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=wp), DIMENSION(Np), INTENT(IN)   :: X, Y, Z                                             ! new basis coordinates (input)
       REAL(KIND=wp), DIMENSION(Np), INTENT(OUT)  :: B                                                  ! Bx,  By, Bz, B components in the x,y,z system

       ! Local variables
       REAL(KIND=wp), DIMENSION(Np)                    :: hx, hy, hz, Bx, By, Bz, psir, psiz, Fpsi            ! lame coefficients
       REAL(KIND=wp), DIMENSION(Nqua, Np)               :: R                                                   ! G_ij = grad(qj).dr/dxi;

       hx = 1.0                                                       ! lame coefficient for x = R -> hx = 1
       hy = 1.0                                                       ! lame coefficient for y = Z -> hy = 1
       hz = X                                                         ! lame coefficient for z = zeta -> hz = R = X

       ! ------------------------------------------------------------------------------------------------------------------------
       ! R contains: d_r\psi, d_z\psi, F(psi), rot(b), grad(B)
       ! ------------------------------------------------------------------------------------------------------------------------

       CALL MagnModel(X, Y, R)

       psir = R(2, :)
       psiz = R(3, :)
       Fpsi = R(7, :)

       ! ------------------------------------------------------------------------------------------------------------------------
       ! Magnetic field components
       ! ------------------------------------------------------------------------------------------------------------------------

       Bx = -psiz/X                                                ! covariant component of B on x  (d_z\psi/R)
       By = psir/X                                                ! covariant component of B on y  (d_rpsi/R)
       Bz = Fpsi/X**2                                             ! covariant component of B on z  (F(psi)/R^2)

       B = sqrt(Bx**2*hx**2 + By**2*hy**2 + Bz**2*hz**2)             ! the norm of the magnetic field (with the covariant components)

    END SUBROUTINE magnetic_norm
