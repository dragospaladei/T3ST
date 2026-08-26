!========================================================================================================
! File: Gyro_center_kernels.f90
! SIMD-friendly gyro-center drift kernels used by the main particle pusher.
! Module name kept as drift_kernels for compatibility.
!========================================================================================================

MODULE drift_kernels
   USE constants
   IMPLICIT NONE

CONTAINS

   PURE SUBROUTINE gyrocenter_drifts(xi, yi, zi, vpi, mui, q1, q2, q3, time, &
                                     vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx, check_1, check_2, check_3, &
                                     vs_1, vs_2, vs_3, Qx1, Qy1, Qz1, Qw1, Qph1, Sol_data)
      !$omp declare simd(gyrocenter_drifts) uniform(time,Qx1,Qy1,Qz1,Qw1,Qph1, Sol_data) notinbranch

      !---------------------------------------------------------------------------------
      ! Arguments
      !---------------------------------------------------------------------------------
      REAL(wp), INTENT(in) :: xi, yi, zi, vpi, mui, time
      REAL(wp), INTENT(out) :: q1, q2, q3, vx, vy, vz, ap, vm, vbF, vbD, kin, Hi, FMi, Pc, B, Vtx, VFx
      REAL(wp), INTENT(out) :: check_1, check_2, check_3, vs_1, vs_2, vs_3
      REAL(wp), INTENT(in) :: Sol_data(SOL_NDATA)
      REAL(wp), INTENT(in), CONTIGUOUS :: Qx1(:), Qy1(:), Qz1(:), Qw1(:), Qph1(:)
      REAL(wp) :: Q0wrp1
      !---------------------------------------------------------------------------------
      ! Locals
      !---------------------------------------------------------------------------------

      ! Fields / drifts
      REAL(wp) :: Esx, Esy, Esz, Bsx, Bsy, Bsz, Bsp, Etx, Ety, Etz
      REAL(wp) :: Bx, By, Bz
      REAL(wp) :: gradBx, gradBy, gradBz
      REAL(wp) :: rotbx, rotby, rotbz
      REAL(wp) :: rotux, rotuy, rotuz
      REAL(wp) :: grad_u2_x, grad_u2_y, grad_u2_z, omega, R02avrg
      REAL(wp) :: phi0x, phi0y, phi0z, phi_ZF, phi_ZF_x

      ! Geometry
      REAL(wp) :: hx, hy, hz, hx2, hy2, hz2
      REAL(wp) :: G(3, 3), F(3, 3), M(3, 3), efit_vals(16)

      ! Coordinates / profiles
      REAL(wp) :: rr, rr2, theta, chi, chir, chiz, rhot, rhotr, rhotz
      REAL(wp) :: rhotrr, rhotrz, rhotzz, psiradius
      REAL(wp) :: qpsi, qprim, psiprim, psiprim2
      REAL(wp) :: psi, psir, psiz, psirr, psizz, psirz
      REAL(wp) :: Fpsi, Fprim, delta_q1

      ! Turbulence accumulators
      REAL(wp) :: phi0, phix, phiy, phiz, phixt, phiyt, phizt
      REAL(wp) :: Apar0, Aparx, Apary, Aparz
      REAL(wp) :: Tprofile, faza, zintc, zints, gpar, gprim
      ! --- hoisted invariants ---
      REAL(wp) :: env1
      REAL(wp) :: q3_over_C3, gamma_delta_q1, shear_factor, flr_pref
      REAL(wp) :: gxx, gyy, gxy
      ! --- SIMD-loop temporaries ---
      REAL(wp) :: amplu, amplu0, wfac, keff_x, keff_y, keff_z, frac, gg
      INTEGER :: ddm
      REAL(wp) :: Vx_1, Vy_1, Vz_1, Ap_1

      ! Numerics / helpers
      REAL(wp) :: xi_m1, xi2, yi2
      REAL(wp) :: one_m_rr2, root1, inv_root1
      REAL(wp) :: invxi, invxi2, invxi3, invhx, invhy, invhz
      REAL(wp) :: invB, invB3, invBsp
      REAL(wp) :: a02, tau, s_star
      REAL(wp) :: tsc, msc, temp, dens, Hi0, msg
      INTEGER :: poz1, poz2
      REAL(wp) :: X1, Y1, Xef, Yef
      REAL(wp) :: F1, F2, F3, F4
      REAL(wp) :: faza0, zintc0, zints0, frac0, cos_m, sin_m
      REAL(wp) :: cos_ddm(-dmmax:dmmax)
      REAL(wp) :: sin_ddm(-dmmax:dmmax)
      INTEGER :: n, jax, i

      ! Solovev flux-coordinate helpers
      REAL(wp) :: psi_sep, delta_sep
      REAL(wp) :: rpsi, rpsi_tab, rpsir, rpsiz
      REAL(wp) :: h2, hsqrt, delta2, delta_geom
      REAL(wp) :: eta, etar, etaz, ceta, seta

      ! Inline radial interpolation
      REAL(wp) :: xsol, fsol
      REAL(wp) :: w0sol, w0sol_r
      REAL(wp) :: Ttor, Jedge
      REAL(wp) :: vL, vR
      REAL(wp) :: bk, bkr

      ! q and straight-angle derivatives
      REAL(wp) :: Cqsol, Frpsi, qrpsi
      REAL(wp) :: chi_rpsi, chi_eta

      ! Fourier recurrence
      REAL(wp) :: sink, cosk, snext, cnext

      INTEGER :: isol, idxL, idxR
      INTEGER :: ibL, ibR, k

      !---------------------------------------------------------------------------------
      ! Quick aliases / hoists
      !---------------------------------------------------------------------------------
      xi_m1 = xi - 1.0_WP + eps_xi
      xi2 = xi*xi
      yi2 = yi*yi

      rr2 = xi_m1*xi_m1 + yi2
      rr = sqrt(rr2) + eps_rr

      one_m_rr2 = max(1.0_WP - rr2, eps_root)
      root1 = sqrt(one_m_rr2)
      inv_root1 = 1.0_WP/root1

      invxi = 1.0_WP/xi
      invxi2 = invxi*invxi
      invxi3 = invxi2*invxi

      a02 = a0*a0
      tau = Aeff*Te/(Ti + Te)
      s_star = (As/Zs)*(rhoi/R0)

      ! Metric (cylindrical-like)
      hx = 1.0_WP; hy = 1.0_WP; hz = xi
      hx2 = 1.0_WP; hy2 = 1.0_WP; hz2 = xi2
      invhx = 1.0_WP/hx
      invhy = 1.0_WP/hy
      invhz = 1.0_WP/hz

      !---------------------------------------------------------------------------------
      ! Geometry and magnetic equilibrium :: Straight-field-line angle, safety factor model, fluxes
      !---------------------------------------------------------------------------------
      IF (magnetic_model == 3) THEN   ! circular equilibrium

         theta = atan2(yi, xi_m1)

         ! chi
         chi = 2.0_WP*atan(sqrt((1.0_WP - rr)/(1.0_WP + rr))*tan(0.5_WP*theta))

         ! Derivatives of chi
         chir = sin(theta)*(rr*rr - xi)*(invxi/rr)*inv_root1
         chiz = -(rr - cos(theta))*(1.0_WP/rr)*inv_root1

         ! rho_t and derivatives
         rhot = rr/a0
         rhotr = xi_m1/(rr*a0)
         rhotz = yi/(rr*a0)
         rhotrr = (1.0_WP/a0)*yi2/(rr*rr*rr)
         rhotrz = -(1.0_WP/a0)*yi*xi_m1/(rr*rr*rr)
         rhotzz = (1.0_WP/a0)*(xi_m1*xi_m1)/(rr*rr*rr)

         ! q(r) via Horner
         qpsi = (s3*rhot + s2)*rhot + s1
         qprim = 2.0_WP*s3*rhot + s2

         ! psi'(r), psi''(r)
         psiprim = a0*rr/(qpsi*root1)
         psiprim2 = -a02*(rr/(qpsi*root1))*(qprim/(a0*qpsi) - 1.0_WP/(one_m_rr2*rr))

         ! Cartesian derivatives of psi
         psir = psiprim*rhotr
         psiz = psiprim*rhotz
         psirr = psiprim*rhotrr + psiprim2*rhotr*rhotr
         psirz = psiprim*rhotrz + psiprim2*rhotr*rhotz
         psizz = psiprim*rhotzz + psiprim2*rhotz*rhotz

         ! F(psi) profile
         Fpsi = 1.0_WP
         Fprim = 0.0_WP

         ! Circular poloidal flux
         psi = a02/s2*(rr/a0 - s1/s2*log(1.0_WP + s2/s1*rr/a0))

      ELSE IF (magnetic_model == 2) THEN

         ! psi and derivatives (analytical model)
         psi = (amp*((alfa**2*(-1.0 + xi**2)**2)/4.0_WP + (-gama + xi**2)*yi**2))/(2.0_WP*(1.0_WP + alfa**2))
         psir = (amp*(alfa**2*xi*(-1.0 + xi**2) + 2.0_WP*xi*yi**2))/(2.0_WP*(1.0_WP + alfa**2))
         psiz = (amp*(-gama + xi**2)*yi)/(1.0_WP + alfa**2)
         psirr = (amp*(2.0_WP*alfa**2*xi**2 + alfa**2*(-1.0 + xi**2) + 2.0_WP*yi**2))/(2.0_WP*(1.0_WP + alfa**2))
         psirz = (2.0_WP*amp*xi*yi)/(1.0_WP + alfa**2)
         psizz = (amp*(-gama + xi**2))/(1.0_WP + alfa**2)

         ! F(psi)
         Fpsi = sqrt(1.0_WP + 2.0_WP*amp*gama*psi/(1.0_WP + alfa**2))
         Fprim = (amp*gama)/(1.0_WP + alfa**2)/Fpsi

         !======================================================================
         ! Solovev flux coordinates
         ! rpsi  = sqrt(psi/psi_edge) : used ONLY as radial interpolation label
         ! rhot  = sqrt(Phi_t/Phi_t_edge) : true toroidal-flux radius
         ! chi   = straight-field-line poloidal angle in radians
         ! Precomputed Sol_data contains:
         !    w0(rpsi), T(rpsi)=rhot^2, b_k(rpsi), Jedge
         !======================================================================

         !----------------------------------------------------------------------
         ! Edge flux.
         !
         ! The chosen outer midplane boundary is R = 1+a0, Z=0.
         !
         ! delta = R^2 - 1 there, and psi = Psi0*delta^2
         !----------------------------------------------------------------------

         delta_sep = 1.0_WP - gama
         psi_sep = amp*alfa**2*delta_sep**2/(8.0_WP*(1.0_WP + alfa**2))

         !----------------------------------------------------------------------
         ! Exact Solovev geometrical quantities
         !
         ! delta^2 =  (R^2-1)^2 + 4 Z^2 (R^2-gamma)/alpha^2
         !   rpsi = delta/delta_edge = sqrt(psi/psi_edge).
         !----------------------------------------------------------------------

         h2 = max(xi2 - gama, eps_root)
         hsqrt = sqrt(h2)
         delta2 = (xi2 - 1.0_WP)**2 + 4.0_WP*yi2*h2/(alfa*alfa)
         delta_geom = sqrt(max(delta2, 0.0_WP))
         rpsi = delta_geom/delta_sep

         ! Clamp only the TABLE LOOKUP coordinate.
         ! Normally particles should satisfy 0 <= rpsi <= SOL_R_MAX.
         rpsi_tab = min(rpsi, SOL_RMAX)
         !======================================================================
         ! Inline interpolation on uniform rpsi grid
         ! xsol runs from 0 to Nsol-1.
         ! isol is a zero-based cell number: 0 ... Nsol-2.
         ! idxL,idxR are Fortran one-based table indices.
         !======================================================================
         xsol = rpsi_tab*SOL_INV_DR

         isol = min(int(xsol), Nsol - 2)
         fsol = xsol - REAL(isol, wp)
         idxL = isol + 1
         idxR = idxL + 1

         !----------------------------------------------------------------------
         ! w0(rpsi)
         ! w0 is even about the magnetic axis:      w0 = w00 + O(rpsi^2)
         ! Therefore use rpsi^2 interpolation in the first radial cell.
         ! Everywhere else use very cheap linear interpolation.
         !----------------------------------------------------------------------

         vL = Sol_data(SOL_W0_OFF + idxL)
         vR = Sol_data(SOL_W0_OFF + idxR)

         IF (isol == 0) THEN
            w0sol = vL + (vR - vL)*fsol*fsol
            w0sol_r = 2.0_WP*(vR - vL)*fsol*SOL_INV_DR
         ELSE
            w0sol = vL + (vR - vL)*fsol
            w0sol_r = (vR - vL)*SOL_INV_DR
         END IF

         !----------------------------------------------------------------------
         ! T(rpsi) = rhot^2
         !
         ! Near the axis T ~ const*rpsi^2, so use the same first-cell
         ! regularization.
         !----------------------------------------------------------------------

         vL = Sol_data(SOL_T_OFF + idxL)
         vR = Sol_data(SOL_T_OFF + idxR)

         IF (isol == 0) THEN
            Ttor = vL + (vR - vL)*fsol*fsol
         ELSE
            Ttor = vL + (vR - vL)*fsol
         END IF

         ! Jedge = integral_0^1 2*rpsi*q(rpsi) drpsi
         Jedge = Sol_data(SOL_JEDGE_IDX)
         !======================================================================
         ! Safety factor q(rpsi)
         !
         ! q = [(1+alpha^2)/(alpha*amp)] F w0
         !
         ! The code uses R_axis = 1 in these normalized coordinates.
         !======================================================================

         Cqsol = (1.0_WP + alfa**2)/(alfa*amp)
         qpsi = Cqsol*Fpsi*w0sol
         !----------------------------------------------------------------------
         ! dq/drpsi
         ! dF/drpsi = (dF/dpsi)*(dpsi/drpsi)  = Fprim * 2*psi_edge*rpsi
         !----------------------------------------------------------------------

         Frpsi = 2.0_WP*psi_sep*rpsi*Fprim
         qrpsi = Cqsol*(Frpsi*w0sol + Fpsi*w0sol_r)

         !======================================================================
         ! True toroidal-flux radius
         !======================================================================

         rhot = sqrt(max(Ttor, 0.0_WP))
         !----------------------------------------------------------------------
         ! rho_t spatial derivatives
         !
         ! d rho_t / d psi =   q /(2 psi_edge Jedge rho_t)
         !
         ! Thus:
         ! rho_t,R = q psi_R /(2 psi_edge Jedge rho_t)
         ! rho_t,Z = q psi_Z /(2 psi_edge Jedge rho_t)
         !
         ! rho_t is singular as a radial coordinate exactly on the magnetic
         ! axis, just like any polar radial/angle coordinate pair.
         !----------------------------------------------------------------------

         IF ((rhot > eps_rr) .AND. (abs(Jedge) > eps_B)) THEN
            rhotr = qpsi*psir/(2.0_WP*psi_sep*Jedge*rhot)
            rhotz = qpsi*psiz/(2.0_WP*psi_sep*Jedge*rhot)
         ELSE

            rhotr = 0.0_WP
            rhotz = 0.0_WP
         END IF
         !----------------------------------------------------------------------
         ! qprim in THIS kernel means dq/drhot, not dq/dpsi.
         ! drhot/drpsi = rpsi*q /(Jedge*rhot)
         ! dq/drhot = Jedge*rhot*(dq/drpsi)/(rpsi*q)
         !----------------------------------------------------------------------

         IF ((rpsi > eps_rr) .AND. (rhot > eps_rr) .AND. &
             (abs(qpsi) > eps_B)) THEN
            qprim = Jedge*rhot*qrpsi/(rpsi*qpsi)
         ELSE
            ! Regular magnetic-axis limit for an even q(rpsi) profile.
            qprim = 0.0_WP
         END IF
         !======================================================================
         ! Exact Solovev geometrical angle eta
         !
         ! cos eta = (R^2-1)/delta
         ! sin eta = 2 Z sqrt(R^2-gamma)/(alpha delta)
         !
         ! This orientation is chosen to match the B_R,B_Z convention
         ! used below in this routine.
         !======================================================================

         IF (delta_geom > eps_rr) THEN

            ceta = (xi2 - 1.0_WP)/delta_geom
            seta = 2.0_WP*yi*hsqrt/(alfa*delta_geom)
            eta = atan2(seta, ceta)
            !---------------------------------------------------------------
            ! Exact eta derivatives
            !---------------------------------------------------------------

            etar = -2.0_WP*xi*yi*(xi2 + 1.0_WP - 2.0_WP*gama)/(alfa*hsqrt*delta2)
            etaz = 2.0_WP*(xi2 - 1.0_WP)*hsqrt/(alfa*delta2)
         ELSE

            ! eta itself is undefined exactly at the magnetic axis.
            eta = 0.0_WP
            ceta = 1.0_WP
            seta = 0.0_WP
            etar = 0.0_WP
            etaz = 0.0_WP
         END IF
         !----------------------------------------------------------------------
         ! rpsi derivatives
         !
         ! psi = psi_edge*rpsi^2
         !----------------------------------------------------------------------

         IF (rpsi > eps_rr) THEN
            rpsir = psir/(2.0_WP*psi_sep*rpsi)
            rpsiz = psiz/(2.0_WP*psi_sep*rpsi)
         ELSE
            rpsir = 0.0_WP
            rpsiz = 0.0_WP
         END IF
         !======================================================================
         ! Straight-field-line angle
         !
         ! chi(rpsi,eta) =
         !
         !     eta + Sum_k [ b_k(rpsi)/k * sin(k eta) ]
         !
         ! chi_eta =
         !
         !     1 + Sum_k b_k cos(k eta)
         !
         ! chi_rpsi =
         !
         !     Sum_k b_k'(rpsi)/k * sin(k eta)
         !
         ! No sin(k*eta), cos(k*eta) calls:
         ! use recurrence from sin(eta),cos(eta).
         !======================================================================

         chi = eta
         chi_rpsi = 0.0_WP
         chi_eta = 1.0_WP

         sink = seta
         cosk = ceta
         DO k = 1, Ksol

            !---------------------------------------------------------------
            ! Inline interpolation of b_k(rpsi)
            !
            ! b_k ~ rpsi^k at the magnetic axis.
            ! Use that known regular behavior inside the first grid cell.
            !---------------------------------------------------------------

            ibL = SOL_B_OFF + (k - 1)*Nsol + idxL
            ibR = ibL + 1

            vL = Sol_data(ibL)
            vR = Sol_data(ibR)

            IF (isol == 0) THEN

               ! b_k(0)=0 by construction
               bk = vR*fsol**k
               bkr = REAL(k, wp)*vR*fsol**(k - 1)*SOL_INV_DR
            ELSE

               bk = vL + (vR - vL)*fsol
               bkr = (vR - vL)*SOL_INV_DR
            END IF
            !---------------------------------------------------------------
            ! Reconstruct chi and its two natural derivatives
            !---------------------------------------------------------------

            chi = chi + bk*sink/REAL(k, wp)
            chi_rpsi = chi_rpsi + bkr*sink/REAL(k, wp)
            chi_eta = chi_eta + bk*cosk

            !---------------------------------------------------------------
            ! sin[(k+1)eta], cos[(k+1)eta] recurrence
            !---------------------------------------------------------------

            IF (k < Ksol) THEN
               snext = sink*ceta + cosk*seta
               cnext = cosk*ceta - sink*seta
               sink = snext
               cosk = cnext
            END IF
         END DO
         !======================================================================
         ! Spatial derivatives of the straight-field-line angle
         !======================================================================

         chir = chi_rpsi*rpsir + chi_eta*etar
         chiz = chi_rpsi*rpsiz + chi_eta*etaz
      ELSE

         IF (magnetic_model == 1) THEN
            ! Grid indices
            poz1 = modulo(int((xi - minR)/stepR), NgridR)
            poz2 = modulo(int((yi - minZ)/stepZ), NgridZ)

            X1 = minR + poz1*stepR
            Y1 = minZ + poz2*stepZ

            Xef = (xi - X1)/stepR
            Yef = (yi - Y1)/stepZ

            !--------------------------------------------------------------
            ! Bilinear interpolation for EFIT quantities
            !
            !  1=psi, 2=psir, 3=psiz, 4=psirr, 5=psirz, 6=psizz,
            !  7=Fpsi, 8=Fprim,
            !  9=qpsi, 10=qprim=dq/drho,
            ! 11=rhot, 12=rhotr, 13=rhotz,
            ! 14=chi  (overwritten below), 15=chir, 16=chiz
            !
            ! Efit_data :: {psi,psir,psiz,psirr,psirz,psizz, F,Fprim,
            !               qpsi,qprim, rhot, rhotR, rhotZ, chi, chir, chiz}
            !--------------------------------------------------------------
            DO jax = 1, 16
               F1 = Efit_data(NgridR*NgridZ*(jax - 1) + poz1*NgridR + poz2)
               F2 = Efit_data(NgridR*NgridZ*(jax - 1) + (poz1 + 1)*NgridR + poz2)
               F3 = Efit_data(NgridR*NgridZ*(jax - 1) + poz1*NgridR + (poz2 + 1))
               F4 = Efit_data(NgridR*NgridZ*(jax - 1) + (poz1 + 1)*NgridR + (poz2 + 1))

               efit_vals(jax) = F1 + (F2 - F1)*Xef + (F3 - F1)*Yef + Xef*Yef*(F1 + F4 - F2 - F3)
            END DO

            psi = efit_vals(1)
            psir = efit_vals(2)
            psiz = efit_vals(3)
            psirr = efit_vals(4)
            psirz = efit_vals(5)
            psizz = efit_vals(6)

            Fpsi = efit_vals(7)
            Fprim = efit_vals(8)

            qpsi = efit_vals(9)
            qprim = efit_vals(10)

            rhot = efit_vals(11)
            rhotr = efit_vals(12)
            rhotz = efit_vals(13)
         END IF

         IF (magnetic_model == 1) THEN
            chi = efit_vals(14)
            chir = efit_vals(15)
            chiz = efit_vals(16)
         END IF

      END IF

      !---------------------------------------------------------------------------------
      ! Scaled field-aligned coords (GENE-like)
      !---------------------------------------------------------------------------------
      q1 = C1*rhot
      q2 = C2*(zi - qpsi*chi)
      q3 = C3*chi
      delta_q1 = q1 - C1*rhot0

      !---------------------------------------------------------------------------------
      ! Grad(qj) contravariant (G)
      !---------------------------------------------------------------------------------
      G(1, 1) = C1*rhotr; G(2, 1) = C1*rhotz; G(3, 1) = 0.0_WP
      G(1, 3) = C3*chir; G(2, 3) = C3*chiz; G(3, 3) = 0.0_WP

      G(1, 2) = -C2*qpsi*chir - C2*chi*qprim*rhotr
      G(2, 2) = -C2*qpsi*chiz - C2*chi*qprim*rhotz
      G(3, 2) = C2

      !---------------------------------------------------------------------------------
      ! Magnetic field (covariant comps) and |B|
      !---------------------------------------------------------------------------------
      Bx = -psiz*invxi
      By = psir*invxi
      Bz = Fpsi*invxi2

      B = sqrt(Bx*Bx*hx2 + By*By*hy2 + Bz*Bz*hz2)
      invB = 1.0_WP/max(B, eps_B)
      invB3 = invB*invB*invB

      !---------------------------------------------------------------------------------
      ! grad|B| and curls
      !---------------------------------------------------------------------------------
      gradBx = -((Fpsi*Fpsi + psir*psir + psiz*psiz) - Fpsi*Fprim*xi*psir - xi*psir*psirr - psiz*xi*psirz)*invB*invxi3
      gradBy = (Fpsi*Fprim*psiz + psir*psirz + psiz*psizz)*invB*invxi2
      gradBz = 0.0_WP

      rotbx = (Fprim*psiz*(psir*psir + psiz*psiz) - Fpsi*(psir*psirz + psiz*psizz))*invB3*invxi3
      rotby = -(Fpsi*Fpsi*Fpsi + Fpsi*(psir*psir + psiz*psiz) &
                - xi*Fpsi*(psirz*psiz + psir*psirr) + xi*Fprim*psir*(psir*psir + psiz*psiz))*invB3*invxi2*invxi2
      rotbz = -(Fpsi*Fprim*(psir*psir + psiz*psiz) - Fpsi*Fpsi*(psirr + psizz) &
                - psirr*psiz*psiz - psizz*psir*psir + 2.0_WP*psirz*psir*psiz)*invB3*invxi2*invxi2

      !---------------------------------------------------------------------------------
      ! u, rotu, gradu, phi0 (e pur si muove)
      !---------------------------------------------------------------------------------

      Omega = Omgt0*(1.0_WP + Omgtprim*(rhot - rhot0)) + eps_omega
      R02avrg = 1.0_WP + a0*rhot/2.0_WP
      rotux = xi*Omgt0*Omgtprim*rhotz
      rotuy = -xi*Omgt0*Omgtprim*rhotr - 2.0_WP*Omega
      rotuz = 0.0_WP

      ! u^2 gradient pieces and neoclassical Phi0 gradient
      grad_u2_x = xi2*Omega*Omega*(rhotr*Omgt0*Omgtprim/Omega + 1.0_WP/xi)
      grad_u2_y = xi2*Omega*Omega*(rhotz*Omgt0*Omgtprim/Omega)
      grad_u2_z = 0.0_WP

      phi0x = tau*(1.0_WP - R02avrg/xi2)*grad_u2_x + tau*Omega*Omega*R02avrg/xi
      phi0y = tau*(1.0_WP - R02avrg/xi2)*grad_u2_y
      phi0z = 0.0_WP

      !---------------------------------------------------------------------------------
      ! B* (covariant) and E*
      !---------------------------------------------------------------------------------
      Bsx = Bx + s_star*(vpi*rotbx + rotux)
      Bsy = By + s_star*(vpi*rotby + rotuy)
      Bsz = Bz + s_star*(vpi*rotbz + rotuz)

      Bsp = (Bsx*Bx*hx2 + Bsy*By*hy2 + Bsz*Bz*hz2)*invB
      invBsp = 1.0_WP/max(Bsp, eps_Bsp)

      Esx = -(phi0x + mui/Zs*gradBx - As/Zs*grad_u2_x)
      Esy = -(phi0y + mui/Zs*gradBy - As/Zs*grad_u2_y)
      Esz = -(phi0z + mui/Zs*gradBz - As/Zs*grad_u2_z)

      ! Baseline (neoclassical) copy for turbulent-only extraction later
      Etx = Esx
      Ety = Esy
      Etz = Esz

      !---------------------------------------------------------------------------------
      ! Cross-product matrix F = (grad(xj) x b) . grad(xi)
      !---------------------------------------------------------------------------------
      F = 0.0_WP
      F(1, 2) = hz*Bz*invB*invhx*invhy
      F(1, 3) = -hy*By*invB*invhx*invhz
      F(2, 3) = hx*Bx*invB*invhy*invhz
      F(2, 1) = -F(1, 2)
      F(3, 1) = -F(1, 3)
      F(3, 2) = -F(2, 3)

      !---------------------------------------------------------------------------------
      ! Turbulent contribution (SIMD)    !!! in the description of turbulence (x,y,z) (q1,q2,q3)= field-aligned; whereas in the remaining (X,Y,Z) = (R,Z,varphi)
      !---------------------------------------------------------------------------------
      phi0 = 0.0_WP; phix = 0.0_WP; phiy = 0.0_WP; phiz = 0.0_WP
      phixt = 0.0_WP; phiyt = 0.0_WP; phizt = 0.0_WP
      Apar0 = 0.0_WP; Aparx = 0.0_WP; Apary = 0.0_WP; Aparz = 0.0_WP
      msg = 0.0_WP

      ! turbulence must be ON and must have started
!      if ((USE_turb == ON) .and. ((time - tt)*sign(1.0_wp, tmax - t0) >= 0.0_wp)) then
      IF ((USE_turb == ON) .AND. ((time - tt) >= 0.0_WP)) THEN

         ! Defaults (kept as in your original snippet)
         q3_over_C3 = q3/C3
         shear_factor = C2*qprim/C1
         flr_pref = rhoi/R0*sqrt(2.0_WP*abs(mui)*As/(Zs**2*B))
         gxx = G(1, 1)*G(1, 1)*invhx**2 + G(2, 1)*G(2, 1)*invhy**2 + G(3, 1)*G(3, 1)*invhz**2
         gyy = G(1, 2)*G(1, 2)*invhx**2 + G(2, 2)*G(2, 2)*invhy**2 + G(3, 2)*G(3, 2)*invhz**2
         gxy = G(1, 1)*G(1, 2)*invhx**2 + G(2, 1)*G(2, 2)*invhy**2 + G(3, 1)*G(3, 2)*invhz**2

         gpar = exp((cos(q3_over_C3) - 1.0_WP)/lbalonz**2)*balloon   ! this must be periodic
         gprim = -sin(q3_over_C3)/lbalonz**2/C3*balloon ! note that this is g'/g
         env1 = noballoon + gpar             ! ballooning envelope multiplier
         IF (turb_model == 1) THEN
            IF (USE_polar == ON) THEN
               !$omp simd private(faza,zintc,zints,amplu,amplu0,wfac,keff_x,keff_y,keff_z,Q0wrp1) &
               !$omp& reduction(+:phi0,phix,phiy,phiz,phixt,phiyt,phizt)
               DO n = 1, Nc

                  wfac = Qw1(n)
                  Q0wrp1 = REAL(INT(C2*Qy1(n)*q00), wp)

                  keff_x = Qx1(n) + Qy1(n)*(q3_over_C3*usetilt*shear_factor)
                  keff_y = Qy1(n)
                  keff_z = Qz1(n) + usetilt*(C2*Qy1(n)*qpsi - Q0wrp1)/C3

                  faza = Qx1(n)*q1 + Qy1(n)*q2 + Qz1(n)*q3 - wfac*time + q3_over_C3*usetilt*(C2*Qy1(n)*qpsi - Q0wrp1) + Qph1(n)

                  amplu = Qx1(n)*Qx1(n)*gxx + Qy1(n)*Qy1(n)*gyy + 2.0_WP*Qx1(n)*Qy1(n)*gxy
                  amplu = flr_pref*sqrt(amplu)
                  ! Approximation for J0(z).
                  amplu = Cos(amplu - Pi/4.0*tanh(amplu**2)) &
                          /(1.0 + Tanh(amplu**2/5.0)*Sqrt(amplu))
                  amplu = USE_larmor*amplu + (1.0 - USE_larmor)
                  amplu = norm*amplu!*exp(-amplu/20_wp)!QL1(n)*exp(-amplu/20_wp)

                  zintc = cos(faza)
                  zints = sin(faza)

                  phi0 = phi0 + amplu*zints
                  phix = phix + amplu*keff_x*zintc
                  phiy = phiy + amplu*keff_y*zintc
                  phiz = phiz + amplu*keff_z*zintc + amplu*gprim*zints

                  phixt = phixt + amplu*keff_x*wfac*zints - amplu*(gamma_E*Qy1(n))*zintc
                  phiyt = phiyt + amplu*keff_y*wfac*zints
                  phizt = phizt + amplu*keff_z*wfac*zints - amplu*gprim*zintc*wfac

               END DO

               phixt = env1*phixt
               phiyt = env1*phiyt
               phizt = env1*phizt
            ELSE
               !$omp simd private(faza,zintc,zints,amplu,amplu0,wfac,keff_x,keff_y,keff_z,Q0wrp1) &
               !$omp& reduction(+:phi0,phix,phiy,phiz)
               DO n = 1, Nc

                  wfac = Qw1(n)
                  Q0wrp1 = REAL(INT(C2*Qy1(n)*q00), wp)

                  keff_x = Qx1(n) + Qy1(n)*(q3_over_C3*usetilt*shear_factor)
                  keff_y = Qy1(n)
                  keff_z = Qz1(n) + usetilt*(C2*Qy1(n)*qpsi - Q0wrp1)/C3

                  faza = Qx1(n)*q1 + Qy1(n)*q2 + Qz1(n)*q3 - wfac*time + q3_over_C3*usetilt*(C2*Qy1(n)*qpsi - Q0wrp1) + Qph1(n)

                  amplu = Qx1(n)*Qx1(n)*gxx + Qy1(n)*Qy1(n)*gyy + 2.0_WP*Qx1(n)*Qy1(n)*gxy
                  amplu = flr_pref*sqrt(amplu)
                  ! Approximation for J0(z).
                  amplu = Cos(amplu - Pi/4.0*tanh(amplu**2)) &
                          /(1.0 + Tanh(amplu**2/5.0)*Sqrt(amplu))
                  amplu = USE_larmor*amplu + (1.0 - USE_larmor)
                  amplu = norm*amplu!*exp(-amplu/20_wp)!QL1(n)*exp(-amplu/20_wp)

                  zintc = cos(faza)
                  zints = sin(faza)

                  phi0 = phi0 + amplu*zints
                  phix = phix + amplu*keff_x*zintc
                  phiy = phiy + amplu*keff_y*zintc
                  phiz = phiz + amplu*keff_z*zintc + amplu*gprim*zints

               END DO
            END IF

            phi0 = env1*phi0
            phix = env1*phix
            phiy = env1*phiy
            phiz = env1*phiz

         ELSEIF (turb_model == 2) THEN
! precompute once per call (or keep cached in a module if dmmax is fixed)
            DO i = -dmmax, dmmax
               cos_ddm(i) = cos(q3_over_C3*REAL(i, wp))
               sin_ddm(i) = sin(q3_over_C3*REAL(i, wp))
            END DO

            IF (USE_polar == ON) THEN
               !$omp simd private(faza0,zintc0,zints0,amplu,amplu0,wfac,gg,keff_x,keff_y,keff_z,frac0,frac,cos_m,sin_m,zintc,zints) &
               !$omp& reduction(+:phi0,phix,phiy,phiz,phixt,phiyt,phizt)
               DO n = 1, Nc
                  keff_x = Qx1(n) + Qy1(n)*(q3_over_C3*shear_factor)
                  keff_y = Qy1(n)
                  wfac = Qw1(n)

                  frac0 = C2*Qy1(n)*qpsi - int(C2*Qy1(n)*qpsi)    ! (same as your original)
                  faza0 = Qx1(n)*q1 + Qy1(n)*q2 + Qz1(n)*q3 - wfac*time + q3_over_C3*frac0 + Qph1(n)

                  zintc0 = cos(faza0)
                  zints0 = sin(faza0)

                  amplu0 = Qx1(n)*Qx1(n)*gxx + Qy1(n)*Qy1(n)*gyy + 2.0_WP*Qx1(n)*Qy1(n)*gxy
                  amplu0 = flr_pref*sqrt(amplu0)
                  ! Approximation for J0(z).
                  amplu0 = Cos(amplu0 - Pi/4.0*tanh(amplu0**2)) &
                           /(1.0 + Tanh(amplu0**2/5.0)*Sqrt(amplu0))
                  amplu0 = USE_larmor*amplu0 + (1.0 - USE_larmor)
                  amplu0 = norm*amplu0!*exp(-amplu/20_wp)!QL1(n)*exp(-amplu/20_wp)

                  DO ddm = -dmmax, dmmax
                     frac = frac0 + ddm
                     gg = exp(-frac**2*lbalonz**2/2.0_WP)*lbalonz/sqrt(2.0_WP*pi)
                     amplu = amplu0*gg

                     keff_z = Qz1(n) + frac/C3

                     cos_m = cos_ddm(ddm)
                     sin_m = sin_ddm(ddm)

                     zintc = zintc0*cos_m - zints0*sin_m
                     zints = zints0*cos_m + zintc0*sin_m

                     phi0 = phi0 + amplu*zints

                     phix = phix + amplu*keff_x*zintc + amplu*(-frac*lbalonz**2)*Qy1(n)*shear_factor*zints
                     phiy = phiy + amplu*keff_y*zintc
                     phiz = phiz + amplu*keff_z*zintc

                     phixt = phixt + amplu*keff_x*wfac*zints - amplu*(gamma_E*Qy1(n))*zintc &
                             - amplu*wfac*(-frac*lbalonz**2)*Qy1(n)*shear_factor*zintc
                     phiyt = phiyt + amplu*keff_y*wfac*zints
                     phizt = phizt + amplu*keff_z*wfac*zints
                  END DO
               END DO

            ELSE
               !$omp simd private(faza0,zintc0,zints0,amplu,amplu0,wfac,gg,keff_x,keff_y,keff_z,frac0,frac,cos_m,sin_m,zintc,zints) &
               !$omp& reduction(+:phi0,phix,phiy,phiz)
               DO n = 1, Nc
                  keff_x = Qx1(n) + Qy1(n)*(q3_over_C3*shear_factor)
                  keff_y = Qy1(n)
                  wfac = Qw1(n)

                  frac0 = C2*Qy1(n)*qpsi - int(C2*Qy1(n)*qpsi)    ! (same as your original)
                  faza0 = Qx1(n)*q1 + Qy1(n)*q2 + Qz1(n)*q3 - wfac*time + q3_over_C3*frac0 + Qph1(n)

                  zintc0 = cos(faza0)
                  zints0 = sin(faza0)

                  amplu0 = Qx1(n)*Qx1(n)*gxx + Qy1(n)*Qy1(n)*gyy + 2.0_WP*Qx1(n)*Qy1(n)*gxy
                  amplu0 = flr_pref*sqrt(amplu0)
                  ! Approximation for J0(z).
                  amplu0 = Cos(amplu0 - Pi/4.0*tanh(amplu0**2)) &
                           /(1.0 + Tanh(amplu0**2/5.0)*Sqrt(amplu0))
                  amplu0 = USE_larmor*amplu0 + (1.0 - USE_larmor)
                  amplu0 = norm*amplu0!exp(-amplu/20_wp)

                  DO ddm = -dmmax, dmmax
                     frac = frac0 + ddm
                     gg = exp(-frac**2*lbalonz**2/2.0_WP)*lbalonz/sqrt(2.0_WP*pi)
                     amplu = amplu0*gg

                     keff_z = Qz1(n) + frac/C3

                     cos_m = cos_ddm(ddm)
                     sin_m = sin_ddm(ddm)

                     zintc = zintc0*cos_m - zints0*sin_m
                     zints = zints0*cos_m + zintc0*sin_m

                     phi0 = phi0 + amplu*zints
                     phix = phix + amplu*keff_x*zintc + amplu*(-frac*lbalonz**2)*Qy1(n)*shear_factor*zints
                     phiy = phiy + amplu*keff_y*zintc
                     phiz = phiz + amplu*keff_z*zintc
                  END DO
               END DO
            END IF

         END IF ! turbulence+model
         ! Radial envelope profile: mainly to avoid numerical problems at r=0 or turbulent runaway particles at the edge
         Tprofile = tanh((turbprof*turbprof)*((rhot - 1.0_WP)**2)*(rhot**2))    ! radial profile of turbulence amplitude
         ! Scale to the mid-radius value. Its radial derivative is not included in the ExB drift.
         Tprofile = Tprofile/tanh((turbprof*turbprof + 0.000001_WP)/16.0_WP)

         ! Zonal flow shearing effects; the average zonal flow was absorbed via a Galilean referance frame change
         phi_ZF = -rhoi/R0/(0.000001_WP + Phi)*gamma_E/2.0_WP*delta_q1**2
         phi_ZF_x = -rhoi/R0/(0.000001_WP + Phi)*gamma_E*delta_q1

         phi0 = phi0 + phi_ZF
         phix = phix + phi_ZF_x
         ! Apply turbulent strength
         tsc = Tprofile*Phi
         Esx = Esx - tsc*(G(1, 1)*phix + G(1, 2)*phiy + G(1, 3)*phiz)
         Esy = Esy - tsc*(G(2, 1)*phix + G(2, 2)*phiy + G(2, 3)*phiz)
         Esz = Esz - tsc*(G(3, 1)*phix + G(3, 2)*phiy + G(3, 3)*phiz)

         IF (USE_polar == ON) THEN
            !---------------------------------------------------------------------------------
            ! Matrix M = (grad(qj) x b) . dr/dxi
            !---------------------------------------------------------------------------------
            M(1, 1) = hx*(Bz*hz*G(2, 1)*invhy - By*hy*G(3, 1)*invhz)*invB
            M(1, 2) = hx*(Bz*hz*G(2, 2)*invhy - By*hy*G(3, 2)*invhz)*invB
            M(1, 3) = hx*(Bz*hz*G(2, 3)*invhy - By*hy*G(3, 3)*invhz)*invB

            M(2, 1) = hy*(Bx*hx*G(3, 1)*invhz - Bz*hz*G(1, 1)*invhx)*invB
            M(2, 2) = hy*(Bx*hx*G(3, 2)*invhz - Bz*hz*G(1, 2)*invhx)*invB
            M(2, 3) = hy*(Bx*hx*G(3, 3)*invhz - Bz*hz*G(1, 3)*invhx)*invB

            M(3, 1) = hz*(By*hy*G(1, 1)*invhx - Bx*hx*G(2, 1)*invhy)*invB
            M(3, 2) = hz*(By*hy*G(1, 2)*invhx - Bx*hx*G(2, 2)*invhy)*invB
            M(3, 3) = hz*(By*hy*G(1, 3)*invhx - Bx*hx*G(2, 3)*invhy)*invB
            msc = tsc*(As/Zs)*(rhoi/R0)

            Esx = Esx + (msc*(M(1, 1)*phixt + M(1, 2)*phiyt + M(1, 3)*phizt))*invB
            Esy = Esy + (msc*(M(2, 1)*phixt + M(2, 2)*phiyt + M(2, 3)*phizt))*invB
            Esz = Esz + (msc*(M(3, 1)*phixt + M(3, 2)*phiyt + M(3, 3)*phizt))*invB
         END IF

         IF (USE_magnturb == ON) THEN
            msg = Tprofile
            Apar0 = phi0
            Aparx = phix
            Apary = phiy
            Aparz = phiz
            Esx = Esx + Beta*msg*vpi*(G(1, 1)*Aparx + G(1, 2)*Apary + G(1, 3)*Aparz)
            Esy = Esy + Beta*msg*vpi*(G(2, 1)*Aparx + G(2, 2)*Apary + G(2, 3)*Aparz)
            Esz = Esz + Beta*msg*vpi*(G(3, 1)*Aparx + G(3, 2)*Apary + G(3, 3)*Aparz)
         END IF

      END IF  ! use_turb condition

      ! Turbulent-only E (difference from baseline)
      Etx = Esx - Etx
      Ety = Esy - Ety
      Etz = Esz - Etz

      !---------------------------------------------------------------------------------
      ! Drifts & accelerations
      !---------------------------------------------------------------------------------
      vx = (vpi - Beta*Zs/As*Apar0)*Bsx*invBsp + (rhoi/R0)*invBsp*(F(1, 1)*Esx + F(1, 2)*Esy + F(1, 3)*Esz)
      vy = (vpi - Beta*Zs/As*Apar0)*Bsy*invBsp + (rhoi/R0)*invBsp*(F(2, 1)*Esx + F(2, 2)*Esy + F(2, 3)*Esz)
      vz = (vpi - Beta*Zs/As*Apar0)*Bsz*invBsp + (rhoi/R0)*invBsp*(F(3, 1)*Esx + F(3, 2)*Esy + F(3, 3)*Esz)

      ap = (Zs/As)*(Esx*Bsx + Esy*Bsy + Esz*Bsz)*invBsp

      !---------------------------------------------------------------------------------
      ! Drifts & accelerations::: the perturbative (turbulent) components
      !---------------------------------------------------------------------------------
      Vx_1 = (-Beta*Zs/As*Apar0)*Bsx*invBsp + (rhoi/R0)*invBsp*(F(1, 1)*Etx + F(1, 2)*Ety + F(1, 3)*Etz)
      Vy_1 = (-Beta*Zs/As*Apar0)*Bsy*invBsp + (rhoi/R0)*invBsp*(F(2, 1)*Etx + F(2, 2)*Ety + F(2, 3)*Etz)
      Vz_1 = (-Beta*Zs/As*Apar0)*Bsz*invBsp + (rhoi/R0)*invBsp*(F(3, 1)*Etx + F(3, 2)*Ety + F(3, 3)*Etz)
      ap_1 = (Zs/As)*(Etx*Bsx + Ety*Bsy + Etz*Bsz)*invBsp

      !---------------------------------------------------------------------------------
      ! purely turbulent drifts along rho for other purposes;
      !---------------------------------------------------------------------------------
      Vtx = vx_1*G(1, 1) + vy_1*G(2, 1) + vz_1*G(3, 1)
      VFx = vx*G(1, 1) + vy*G(2, 1) + vz*G(3, 1)

      ! collisions
      vm = 0.0_WP

      !---------------------------------------------------------------------------------
      ! Energy (Hamiltonian)
      !---------------------------------------------------------------------------------
      Hi = As*(vpi - Zs/As*Apar0)**2*0.5_WP + mui*B - As*xi2*Omega*Omega*0.5_WP &
           + Zs*xi2*Omega*Omega*0.5_WP*tau*(1.0_WP - R02avrg/xi2) &
           + Zs*Phi*phi0
      Hi0 = As*vpi*vpi*0.5_WP + mui*B - As*xi2*Omega*Omega*0.5_WP &
            + Zs*xi2*Omega*Omega*0.5_WP*tau*(1.0_WP - R02avrg/xi2)
      Pc = psi - As/Zs*rhoi/R0*Fpsi/B*(vpi + 1.0_WP*Fpsi/B*Omega)
      kin = As*0.5_WP*(vpi - Zs/As*Apar0)**2 + mui*B   ! pure kinetic energy, without rotational or potential contributions
      check_1 = phi0
      check_2 = phix
      check_3 = phiy

      check_1 = Bz/max(B, eps_B)

      check_2 = qpsi*(Bx*chir + By*chiz)/max(B, eps_B)

      check_3 = -chi*qprim*(Bx*rhotr + By*rhotz)/B

      !---------------------------------------------------------------------------------
      ! vW
      !---------------------------------------------------------------------------------

      temp = Ts*exp(delta_q1*a0/C1*Lts) !Ts*(1.0_wp + delta_q1*a0/C1*Lts)
      dens = 1.0*exp(delta_q1*a0/C1*Lns) !1.0*(1.0_wp + delta_q1*a0/C1*Lns)
      FMi = dens*exp(-Hi0/temp)/temp**1.5_WP

      ! The outer minus comes from rhs=-dfM/dt; the inner minus comes from the Maxwellian -E/T.
      vs_1 = -(-ap_1*As*vpi/temp &
               - mui/temp*(vx_1*gradBx + vy_1*gradBy + vz_1*gradBz))
      vs_2 = -(+Vtx)
      vs_3 = -(+Vtx*(Hi0/temp - 1.5_WP))

      vbD = -(-ap_1*As*vpi/temp &
              - mui/temp*(vx_1*gradBx + vy_1*gradBy + vz_1*gradBz) &
              + Vtx*a0/C1*(Lns + Lts*(Hi0/temp - 1.5_WP)))
      vbF = -(-ap*As*vpi/temp &
              - mui/temp*(vx*gradBx + vy*gradBy + vz*gradBz) &
              + VFx*a0/C1*(Lns + Lts*(Hi0/temp - 1.5_WP)))

   END SUBROUTINE gyrocenter_drifts

   PURE SUBROUTINE collision_kicks(sm1, sm2, dt_local, xi, yi, zi, vpi, mui, B, vcolx, vcoly, vcolz, vcolm, vcolp)
      !$omp declare simd(collision_kicks) uniform(dt_local) notinbranch

      !---------------------------------------------------------------------------------
      ! Arguments
      !---------------------------------------------------------------------------------
      REAL(wp), INTENT(in) :: sm1, sm2, dt_local, xi, yi, zi, vpi, mui, B
      REAL(wp), INTENT(out) :: vcolx, vcoly, vcolz, vcolm, vcolp

      !---------------------------------------------------------------------------------
      ! Locals
      !---------------------------------------------------------------------------------
      REAL(wp) :: v, zeta, ene, ene2, xx, Fphi, Fpsi, sm, ce
      REAL(wp) :: nub, nue, nueprim, nub1, nue1, nueprim1, gg, d

      !---------------------------------------------------------------------------------
      ! collision_operator
      !---------------------------------------------------------------------------------
      ene = As/2.0_WP*vpi**2 + B*abs(mui)
      v = sqrt(2.0_WP*abs(ene)/As)
      zeta = vpi/v

      xx = sqrt(abs(ene))*sqrt(Aeff)/sqrt(As)
      Fphi = erf(xx)
      Fpsi = (Fphi - 2.0_WP*xx*exp(-xx**2)/sqrt(pi))/(2.0_WP*xx**2)

      gg = (16.0_WP*R0*c0)/(vth**4*(2.0_WP/As)**1.5_WP)
      nub1 = (gg/4.0_WP)/ene**1.5_WP*(Fphi - Fpsi)
      nue1 = gg/ene**1.5_WP*Fpsi
      nueprim1 = gg*sqrt(Aeff/As/pi)*exp(-ene*Aeff/As)/ene**2 - 2.5_WP*nue1/ene

      d = delta
      ce = 1.0_WP/(1.0_WP + d*ene**(-1.0_WP)*(As/2.0_WP))

      nub = nub1*ce
      nue = nue1*ce
      nueprim = nueprim1*ce + nue1*2.0_WP*As*d/((As*d + 2.0_WP*ene)**2.0_WP)

      sm = 1.0_WP*sign(1.0_WP, sm1 - 0.5_WP)
      zeta = zeta*(1.0_WP - 1.0_WP*nub*dt_local) + sm*sqrt(dt_local*nub*(1.0_WP - zeta**2))

      sm = 1.0_WP*sign(1.0_WP, sm2 - 0.5_WP)
      ene2 = ene - 1.0_WP*nue*dt_local*(ene - (1.5_WP + ene/(nue + 0.000001_WP)*nueprim)) + 2.0_WP*sm*sqrt(ene*nue*dt_local)
      ene = abs(ene2)

      v = sqrt(2.0_WP*abs(ene)/As)

      !---------------------------------------------------------------------------------
      ! Defaults
      !---------------------------------------------------------------------------------
      vcolx = 0.0_WP; vcoly = 0.0_WP; vcolz = 0.0_WP
      vcolp = (v*zeta - vpi)/dt_local
      vcolm = (ene/B*(1.0_WP - zeta**2) - mui)/(dt_local + 0.00001_WP)

   END SUBROUTINE collision_kicks

END MODULE drift_kernels
