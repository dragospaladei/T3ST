module drift_kernels
   use constants
   implicit none

contains

   pure subroutine Drift2(dt, xi, yi, zi, vpi, mui, q1, q2, q3, time, &
                          vx, vy, vz, ap, vm, Hi,Pc, B, Vtx, Vty, check_1, check_2, check_3, &
                          Qx1, Qy1, Qz1, Qw1, Qph1, QL1, anormal)
      !$omp declare simd(Drift2) uniform(dt,anormal,time,Qx1,Qy1,Qz1,Qw1,Qph1,QL1) notinbranch

      !---------------------------------------------------------------------------------
      ! Arguments
      !---------------------------------------------------------------------------------
      real(rp), intent(in)  :: dt, xi, yi, zi, vpi, mui, time, anormal
!      real(rp), value, intent(in) :: 
      real(rp), intent(out) :: q1, q2, q3, vx, vy, vz, ap, vm, Hi, Pc, B, Vtx, Vty, check_1, check_2, check_3
      real(rp), intent(in), contiguous :: Qx1(:), Qy1(:), Qz1(:), Qw1(:), Qph1(:), QL1(:)
      real(rp) ::  Q0wrp1
      !---------------------------------------------------------------------------------
      ! Locals
      !---------------------------------------------------------------------------------

      ! Fields / drifts
      real(rp) :: Esx, Esy, Esz, Bsx, Bsy, Bsz, Bsp, Etx, Ety, Etz
      real(rp) :: Bx, By, Bz
      real(rp) :: gradBx, gradBy, gradBz
      real(rp) :: rotbx, rotby, rotbz
      real(rp) :: rotux, rotuy, rotuz
      real(rp) :: grad_u2_x, grad_u2_y, grad_u2_z, omega, R02avrg
      real(rp) :: phi0x, phi0y, phi0z

      ! Geometry
      real(rp) :: hx, hy, hz, hx2, hy2, hz2
      real(rp) :: G(3, 3), F(3, 3), M(3, 3), efit_vals(16)

      ! Coordinates / profiles
      real(rp) :: rr, rr2, theta, chi, chir, chiz, rhot, rhotr, rhotz
      real(rp) :: rhotrr, rhotrz, rhotzz, psiradius
      real(rp) :: qpsi, qprim, psiprim, psiprim2
      real(rp) :: psi, psir, psiz, psirr, psizz, psirz
      real(rp) :: Fpsi, Fprim, delta_q1, dBdt

      ! Turbulence accumulators
      real(rp) :: phi0, phix, phiy, phiz, phixt, phiyt, phizt
      real(rp) :: Tprofile, faza, zintc, zints, gpar, gprim
      ! --- hoisted invariants ---
	real(rp) :: env1, env2
	real(rp) :: ge_t
	! --- SIMD-loop temporaries ---
	real(rp) :: amplu, amplu0, wfac, keff_x,keff_y,keff_z, frac, gg
	integer  ::ddm

      ! Numerics / helpers
      real(rp) :: xi_m1, xi2, yi2
      real(rp) :: one_m_rr2, root1, inv_root1
      real(rp) :: invxi, invxi2, invxi3, invhx, invhy, invhz
      real(rp) :: invB, invB3, invBsp
      real(rp) :: Omgt02, a02, tau, s_star
      real(rp) :: t_tmp, sqrt_t1
      real(rp) :: tsc, msc
      integer  :: poz1, poz2
      real(rp) :: X1, Y1, Xef, Yef
      real(rp) :: F1, F2, F3, F4
      real(rp) :: faza0,zintc0,zints0, frac0, cos_m, sin_m
      real(rp) :: cos_ddm(-dmmax:dmmax)
      real(rp) :: sin_ddm(-dmmax:dmmax)
      integer :: n, jax, i

      !---------------------------------------------------------------------------------
      ! Quick aliases / hoists
      !---------------------------------------------------------------------------------
      xi_m1 = xi - 1.0_rp + eps_xi
      xi2   = xi*xi
      yi2   = yi*yi

      rr2 = xi_m1*xi_m1 + yi2
      rr  = sqrt(rr2) + eps_rr

      one_m_rr2  = max(1.0_rp - rr2, eps_root)
      root1      = sqrt(one_m_rr2)
      inv_root1  = 1.0_rp / root1

      invxi  = 1.0_rp / xi
      invxi2 = invxi*invxi
      invxi3 = invxi2*invxi

      a02    = a0*a0
      Omgt02 = Omgt0*Omgt0
      tau    = Aeff*Te / (Ti + Te)
      s_star = (As/Zs) * (rhoi/R0)

      ! Metric (cylindrical-like)
      hx = 1.0_rp;  hy = 1.0_rp;  hz = xi
      hx2 = 1.0_rp; hy2 = 1.0_rp; hz2 = xi2
      invhx = 1.0_rp / hx
      invhy = 1.0_rp / hy
      invhz = 1.0_rp / hz

      !---------------------------------------------------------------------------------
      ! Straight-field-line angle, safety factor model, fluxes
      !---------------------------------------------------------------------------------
      if (magnetic_model == 4) then

         theta = atan2(yi, xi_m1)

         ! chi
         chi = 2.0_rp * atan( sqrt((1.0_rp - rr)/(1.0_rp + rr)) * tan(0.5_rp*theta) )

         ! Derivatives of chi
         chir =  sin(theta) * (rr*rr - xi) * (invxi/rr) * inv_root1
         chiz = -(rr - cos(theta)) * (1.0_rp/rr) * inv_root1

         ! rho_t and derivatives
         rhot   = rr/a0
         rhotr  = xi_m1/(rr*a0)
         rhotz  = yi/(rr*a0)
         rhotrr = (1.0_rp/a0) * yi2/(rr*rr*rr)
         rhotrz = - (1.0_rp/a0) * yi*xi_m1/(rr*rr*rr)
         rhotzz = (1.0_rp/a0) * (xi_m1*xi_m1)/(rr*rr*rr)

         ! q(r) via Horner
         qpsi  = (s3*rhot + s2)*rhot + s1
         qprim = 2.0_rp*s3*rhot + s2

         ! psi'(r), psi''(r)
         psiprim  = a0*rr/(qpsi*root1)
         psiprim2 = -a02*(rr/(qpsi*root1)) * (qprim/(a0*qpsi) - 1.0_rp/(one_m_rr2*rr))

         ! Cartesian derivatives of psi
         psir  = psiprim*rhotr
         psiz  = psiprim*rhotz
         psirr = psiprim*rhotrr + psiprim2*rhotr*rhotr
         psirz = psiprim*rhotrz + psiprim2*rhotr*rhotz
         psizz = psiprim*rhotzz + psiprim2*rhotz*rhotz

         ! F(psi) profile
         Fpsi  = 1.0_rp
         Fprim = 0.0_rp

         ! psi (closed form) kept, algebraically simplified
         t_tmp   = a02*s1/s3
         sqrt_t1 = sqrt(t_tmp + 1.0_rp)
         psi     = -(a02/s3)/sqrt_t1 * (atanh(root1/sqrt_t1) - atanh(1.0_rp/sqrt_t1))

      else if (magnetic_model == 3) then

         ! psi and derivatives (analytical model)
         psi   = (amp*( (alfa**2*(-1.0 + xi**2)**2)/4.0_rp + (-gama + xi**2)*yi**2 )) / (2.0_rp*(1.0_rp + alfa**2))
         psir  = (amp*( alfa**2*xi*(-1.0 + xi**2) + 2.0_rp*xi*yi**2 )) / (2.0_rp*(1.0_rp + alfa**2))
         psiz  = (amp*(-gama + xi**2)*yi) / (1.0_rp + alfa**2)
         psirr = (amp*( 2.0_rp*alfa**2*xi**2 + alfa**2*(-1.0 + xi**2) + 2.0_rp*yi**2 )) / (2.0_rp*(1.0_rp + alfa**2))
         psirz = (2.0_rp*amp*xi*yi) / (1.0_rp + alfa**2)
         psizz = (amp*(-gama + xi**2)) / (1.0_rp + alfa**2)

         ! F(psi)
         Fpsi  = sqrt(1.0_rp + 2.0_rp*amp*gama*psi/(1.0_rp + alfa**2))
         Fprim = (amp*gama)/(1.0_rp + alfa**2) * Fpsi

         ! chi and its derivatives (analogous to circular model)
         psiradius = (amp*( (alfa**2*(-1.0 + (1.0_rp + a0)**2)**2)/4.0_rp )) / (2.0_rp*(1.0_rp + alfa**2))
         rr        = sqrt(psi/psiradius)*a0
         ! note how for Solovev we use the local definition of effective radius (also nondimensional sqrt(psi/psiedge)
         theta     = atan2(yi, xi_m1)

         chi  = 2.0_rp*atan( sqrt((1.0_rp - rr)/(1.0_rp + rr)) * tan(theta/2.0_rp) )
         chir = sin(theta) * (rr**2 - xi) / xi / rr / sqrt(1.0_rp - rr**2)
         chiz = -(rr - cos(theta)) / rr / sqrt(1.0_rp - rr**2)

         ! effective radius and derivatives
         rhot  = rr/a0
         rhotr = rr/2.0_rp/a0 * psir/psi
         rhotz = rr/2.0_rp/a0 * psiz/psi

         ! q(r) via Horner
         qpsi  = (s3*rhot + s2)*rhot + s1
         qprim = 2.0_rp*s3*rhot + s2
 !        note that this form of qpsi is an approximation and is not consistent with psi, thus, with field-alginement

      else

         if ((magnetic_model == 1) .or. (magnetic_model == 2)) then
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
            do jax = 1, 16
               F1 = Efit_data(NgridR*NgridZ*(jax - 1) +  poz1   *NgridR +  poz2)
               F2 = Efit_data(NgridR*NgridZ*(jax - 1) + (poz1+1)*NgridR +  poz2)
               F3 = Efit_data(NgridR*NgridZ*(jax - 1) +  poz1   *NgridR + (poz2+1))
               F4 = Efit_data(NgridR*NgridZ*(jax - 1) + (poz1+1)*NgridR + (poz2+1))

               efit_vals(jax) = F1 + (F2 - F1)*Xef + (F3 - F1)*Yef + Xef*Yef*(F1 + F4 - F2 - F3)
            end do

            psi   = efit_vals(1)
            psir  = efit_vals(2)
            psiz  = efit_vals(3)
            psirr = efit_vals(4)
            psirz = efit_vals(5)
            psizz = efit_vals(6)

            Fpsi  = efit_vals(7)
            Fprim = efit_vals(8)

            qpsi  = efit_vals(9)
            qprim = efit_vals(10)

            rhot  = efit_vals(11)
            rhotr = efit_vals(12)
            rhotz = efit_vals(13)
         end if

         if (magnetic_model == 1) then
            chi  = efit_vals(14)
            chir = efit_vals(15)
            chiz = efit_vals(16)
         end if

         if (magnetic_model == 2) then
            rr    = sqrt((xi - 1.0_rp)**2 + yi*yi) + 1.0e-7_rp
            theta = atan2(yi, xi - 1.0_rp)

            chi  = 2.0_rp * atan( sqrt((1.0_rp - rr)/(1.0_rp + rr)) * tan(0.5_rp*theta) )
            chir = sin(theta) * (rr*rr - xi) / (xi*rr*sqrt(1.0_rp - rr*rr))
            chiz = -(rr - cos(theta)) / (rr*sqrt(1.0_rp - rr*rr))
         end if

      end if

      !---------------------------------------------------------------------------------
      ! Scaled field-aligned coords (GENE-like)
      !---------------------------------------------------------------------------------
      q1 = C1*rhot
      q2 = C2*(zi - qpsi*chi)
      q3 = C3*chi
      delta_q1 = q1 - C1*q10

      !---------------------------------------------------------------------------------
      ! Grad(qj) contravariant (G)
      !---------------------------------------------------------------------------------
      G(1, 1) = C1*rhotr;  G(2, 1) = C1*rhotz;  G(3, 1) = 0.0_rp
      G(1, 3) = C3*chir;   G(2, 3) = C3*chiz;   G(3, 3) = 0.0_rp

      G(1, 2) = -C2*qpsi*chir - C2*chi*qprim*rhotr
      G(2, 2) = -C2*qpsi*chiz - C2*chi*qprim*rhotz
      G(3, 2) =  C2

      !---------------------------------------------------------------------------------
      ! Magnetic field (covariant comps) and |B|
      !---------------------------------------------------------------------------------
      Bx = -psiz*invxi
      By =  psir*invxi
      Bz =  Fpsi*invxi2

      B    = sqrt(Bx*Bx*hx2 + By*By*hy2 + Bz*Bz*hz2)
      invB = 1.0_rp / max(B, eps_B)
      invB3 = invB*invB*invB

      !---------------------------------------------------------------------------------
      ! grad|B| and curls
      !---------------------------------------------------------------------------------
      gradBx = -((Fpsi*Fpsi + psir*psir + psiz*psiz) - Fpsi*Fprim*xi*psir - xi*psir*psirr - psiz*xi*psirz)* invB * invxi3
      gradBy = (Fpsi*Fprim*psiz + psir*psirz + psiz*psizz) * invB * invxi2
      gradBz = 0.0_rp

      rotbx = (Fprim*psiz*(psir*psir + psiz*psiz) - Fpsi*(psir*psirz + psiz*psizz)) * invB3 * invxi3
      rotby = -(Fpsi*Fpsi*Fpsi + Fpsi*(psir*psir + psiz*psiz) &
                - xi*Fpsi*(psirz*psiz + psir*psirr) + xi*Fprim*psir*(psir*psir + psiz*psiz)) * invB3 * invxi2 * invxi2
      rotbz = -(Fpsi*Fprim*(psir*psir + psiz*psiz) - Fpsi*Fpsi*(psirr + psizz) &
                - psirr*psiz*psiz - psizz*psir*psir + 2.0_rp*psirz*psir*psiz) * invB3 * invxi2 * invxi2

      !---------------------------------------------------------------------------------
      ! u, rotu, gradu, phi0 (e pur si muove)
      !---------------------------------------------------------------------------------

      Omega = Omgt0*(1.0_rp + Omgtprim*(rhot - q10)) + eps_omega
      R02avrg = 1.0_rp + a0*rhot/2.0_rp
      rotux = xi*Omgt0*Omgtprim*rhotz
      rotuy = -xi*Omgt0*Omgtprim*rhotr- 2.0_rp*Omega
      rotuz = 0.0_rp
 
 ! u^2 gradient pieces and neoclassical Phi0 gradient
      grad_u2_x = xi2*Omega*Omega*(rhotr*Omgt0*Omgtprim/Omega + 1.0_rp/xi)!Omgt02*(phi0z_nc*phi0z_nc)*xi2*(rhotr*Omgtprim + invxi)
      grad_u2_y = xi2*Omega*Omega*(rhotz*Omgt0*Omgtprim/Omega)!Omgt02*(phi0z_nc*phi0z_nc)*xi2*(rhotz*Omgtprim)
      grad_u2_z = 0.0_rp
  
      phi0x = tau*(1.0_rp - R02avrg/xi2)*grad_u2_x + tau*Omega*Omega*R02avrg/xi
      phi0y = tau*(1.0_rp - R02avrg/xi2)*grad_u2_y
      phi0z = 0.0_rp
      
      !---------------------------------------------------------------------------------
      ! B* (covariant) and E*
      !---------------------------------------------------------------------------------
      Bsx = Bx + s_star*(vpi*rotbx + rotux)
      Bsy = By + s_star*(vpi*rotby + rotuy)
      Bsz = Bz + s_star*(vpi*rotbz + rotuz)

      Bsp    = (Bsx*Bx*hx2 + Bsy*By*hy2 + Bsz*Bz*hz2) * invB
      invBsp = 1.0_rp / max(Bsp, eps_Bsp)

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
      F = 0.0_rp
      F(1, 2) =  hz*Bz*invB*invhx*invhy
      F(1, 3) = -hy*By*invB*invhx*invhz
      F(2, 3) =  hx*Bx*invB*invhy*invhz
      F(2, 1) = -F(1, 2)
      F(3, 1) = -F(1, 3)
      F(3, 2) = -F(2, 3)
      
      !---------------------------------------------------------------------------------
      ! Zonal flow contribution ; note that by definition phi_ZF(x) is y,z independent;
      !---------------------------------------------------------------------------------
!     phi_ZF_x = 

      !---------------------------------------------------------------------------------
      ! Turbulent contribution (SIMD)    !!! everyone should know this: in the description of turbulence (x,y,z) = field-aligned; whereas in the remaining (X,Y,Z) = (R,Z,varphi)
      !---------------------------------------------------------------------------------
      phi0  = 0.0_rp; phix  = 0.0_rp; phiy  = 0.0_rp; phiz  = 0.0_rp
      phixt = 0.0_rp; phiyt = 0.0_rp; phizt = 0.0_rp
      
      ! turbulence must be ON and must have started
      if ( (USE_turb == ON) .and. (time > tt) ) then

        ! Defaults (kept as in your original snippet)
        gpar    = exp((cos(q3/C3)-1.0_rp)/lbalonz**2)*balloon   ! this must be periodic
        gprim   = -sin(q3/C3)/lbalonz**2/C3*balloon ! note that this is g'/g
        ge_t    = gamma_E * time
        env1    = noballoon + gpar             ! ballooning envelope multiplier
        env2    = anormal ! This is for normalization (outside the nc loop)
     
         if(turb_model == 1) then     
         !$omp simd private(faza,zintc,zints,amplu,amplu0,wfac,keff_x,keff_y,keff_z,Q0wrp1) &
         !$omp& reduction(+:phi0,phix,phiy,phiz,phixt,phiyt,phizt)
         do n = 1, Nc

            wfac  = Qw1(n) + Qy1(n)*gamma_E*delta_q1
            Q0wrp1  = REAL(INT(C2*Qy1(n)*q00), rp)

            ! Common x-derivative factor used in phix and phixt
            keff_x    = Qx1(n) + Qy1(n)*(q3/C3*usetilt*(C2*qprim/C1) - ge_t)
            keff_y    = Qy1(n)
            keff_z    = Qz1(n) + usetilt*(C2*Qy1(n)*qpsi - Q0wrp1)/C3
            
            faza  = Qx1(n)*q1 + Qy1(n)*q2 + Qz1(n)*q3 - wfac*time + q3/C3*usetilt*(C2*Qy1(n)*qpsi - Q0wrp1) + Qph1(n)

	    amplu = 2.0_rp*abs(mui)*As/Zs**2*(keff_x**2+keff_y**2)/B
	    amplu = norm*QL1(n)*exp(-amplu/20_rp)
	           
            ! ifx is typically good at fusing sin/cos, but paired calls are still fine
            zintc = cos(faza)
            zints = sin(faza)

            phi0  = phi0  + amplu*zints

            phix  = phix  + amplu*keff_x*zintc
            phiy  = phiy  + amplu*keff_y*zintc
            phiz  = phiz  + amplu*keff_z*zintc + amplu*gprim*zints

            ! Time-derivatives
            phixt = phixt + amplu*keff_x*wfac*zints - amplu*(gamma_E*Qy1(n))*zintc
            phiyt = phiyt + amplu*keff_y*wfac*zints
            phizt = phizt + amplu*keff_z*wfac*zints - amplu*gprim*zintc*wfac

         end do

         ! Apply the ballooning envelope once
            phi0  = env1*phi0
            phix  = env1*phix
            phiy  = env1*phiy
            phiz  = env1*phiz
            phixt = env1*phixt
            phiyt = env1*phiyt
            phizt = env1*phizt

         elseif(turb_model == 2) then     
! precompute once per call (or keep cached in a module if dmmax is fixed)
	    do i = -dmmax, dmmax
	      cos_ddm(i) = cos(q3 / C3 * real(i,rp))
	      sin_ddm(i) = sin(q3 / C3 * real(i,rp))
	    end do

	    !$omp simd private(faza0,zintc0,zints0,amplu,amplu0,wfac,gg,keff_x,keff_y,keff_z,frac0,frac,cos_m,sin_m,zintc,zints) &
	    !$omp& reduction(+:phi0,phix,phiy,phiz,phixt,phiyt,phizt)
	    do n = 1, Nc
	       keff_x = Qx1(n) + Qy1(n)*( q3/C3*(C2*qprim/C1) - ge_t )
	       keff_y = Qy1(n)
	       wfac   = Qw1(n) + Qy1(n)*gamma_E*delta_q1

		  frac0  = C2*Qy1(n)*qpsi - int(C2*Qy1(n)*qpsi)    ! (same as your original)
		  faza0  = Qx1(n)*q1 + Qy1(n)*q2 + Qz1(n)*q3 - wfac*time + (q3/C3)*frac0 + Qph1(n)

		  zintc0 = cos(faza0)
		  zints0 = sin(faza0)

!		    amplu = QL1(n)
	            amplu0 = 2.0_rp*abs(mui)*As/Zs**2*(keff_x**2+keff_y**2)/B
	            amplu0 = QL1(n)!*exp(-amplu0/20_rp)
!	            amplu0 = 1.0_rp
!	            amplu = Bessel_J0(amplu)
	            
		  do ddm = -dmmax, dmmax
		    frac  = frac0 + ddm
		    gg    = exp(-frac**2*lbalonz**2/2.0_rp)
		    amplu = amplu0*gg

		    keff_z = Qz1(n) + frac/C3

		    cos_m = cos_ddm(ddm)
		    sin_m = sin_ddm(ddm)

		    zintc = zintc0*cos_m - zints0*sin_m
		    zints = zints0*cos_m + zintc0*sin_m

		    phi0  = phi0  + amplu*zints

		    phix  = phix  + amplu*keff_x*zintc + amplu*(-frac*lbalonz**2)*C2*Qy1(n)*qprim/C1*zints
		    phiy  = phiy  + amplu*keff_y*zintc
		    phiz  = phiz  + amplu*keff_z*zintc

		    phixt = phixt + amplu*keff_x*wfac*zints - amplu*(gamma_E*Qy1(n))*zintc &
				  - amplu*wfac*(-frac*lbalonz**2)*C2*Qy1(n)*qprim/C1*zintc
		    phiyt = phiyt + amplu*keff_y*wfac*zints
		    phizt = phizt + amplu*keff_z*wfac*zints
		  end do
		end do
! Apply the normalization once
         phi0  = env2*phi0
         phix  = env2*phix
         phiy  = env2*phiy
         phiz  = env2*phiz
         phixt = env2*phixt
         phiyt = env2*phiyt
         phizt = env2*phizt

         endif ! turbulence+model
      Tprofile = tanh((turbprof*turbprof) * ((rhot - 1.0_rp)**2) * (rhot**2))    ! radial profile of turbulence amplitude
      Tprofile = Tprofile/tanh((turbprof*turbprof+0.000001_rp)/16.0_rp)                      ! scaled to mid-radius value; it's radial derivative is NOT taken into account in the ExB drift    

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
      ENDIF

      end if  ! use_turb condition
      
      ! Turbulent-only E (difference from baseline)
      Etx = Esx - Etx
      Ety = Esy - Ety
      Etz = Esz - Etz

      !---------------------------------------------------------------------------------
      ! Drifts & accelerations
      !---------------------------------------------------------------------------------
      vx = vpi*Bsx*invBsp + (rhoi/R0)*invBsp*(F(1, 1)*Esx + F(1, 2)*Esy + F(1, 3)*Esz)
      vy = vpi*Bsy*invBsp + (rhoi/R0)*invBsp*(F(2, 1)*Esx + F(2, 2)*Esy + F(2, 3)*Esz)
      vz = vpi*Bsz*invBsp + (rhoi/R0)*invBsp*(F(3, 1)*Esx + F(3, 2)*Esy + F(3, 3)*Esz)

      ap = (Zs/As) * (Esx*Bsx + Esy*Bsy + Esz*Bsz) * invBsp

      Vtx = (rhoi/R0)*invBsp*(F(1, 1)*Etx + F(1, 2)*Ety + F(1, 3)*Etz)
      Vty = (rhoi/R0)*invBsp*(F(2, 1)*Etx + F(2, 2)*Ety + F(2, 3)*Etz)

      Vtx = C1*(Vtx*rhotr + Vty*rhotz)
      Vty = C1*(Vx*rhotr  + Vy*rhotz)

      ! collisions
      vm = 0.0_rp

      !---------------------------------------------------------------------------------
      ! Energy (Hamiltonian)
      !---------------------------------------------------------------------------------
!     dBdt = gradBx*vx*hx2 + gradBy*vy*hy2 + gradBz*vz*hz2   ! if that’s your consistent “physical dot”
!     Hi = As*vpi*ap + mui*dBdt! - As u^2/2.0
     Hi = As*vpi*vpi*0.5_rp + mui*B - As*xi2*Omega*Omega*0.5_rp + Zs*xi2*Omega*Omega*0.5_rp*tau*(1.0_rp - R02avrg/xi2) + Zs*Phi*phi0!- (1.0_rp - tau)*As*(Omgt02*xi2)*0.5_rp + Zs*Phi*phi0
     Pc = psi - As/Zs*rhoi/R0*Fpsi/B*(vpi + 1.0_rp*Fpsi/B*Omega)
     check_1 = phi0
     check_2 = phix
     check_3 = phiy
      
! Hi = gradBy/B
! Use the following if you wanna test toroidal canonical momentum conservation Pctor
! there are several mistakes related to the inclusion of rotation

!q1 = gradBx/B
!q2 = gradBy/B
!q3 = rotbx
!Hi = rotby
!q1=psir
!q2=psiz
!q3=psirr
!Hi=psirz
end subroutine Drift2

   pure subroutine collisions_1_MP(sm1, sm2, dt_local, xi, yi, zi, vpi, mui, B, vcolx, vcoly, vcolz, vcolm, vcolp)
      !$omp declare simd(collisions_1_MP) uniform(dt_local) notinbranch

      !---------------------------------------------------------------------------------
      ! Arguments
      !---------------------------------------------------------------------------------
      real(rp), intent(in)  :: sm1, sm2, dt_local, xi, yi, zi, vpi, mui, B
      real(rp), intent(out) :: vcolx, vcoly, vcolz, vcolm, vcolp

      !---------------------------------------------------------------------------------
      ! Locals
      !---------------------------------------------------------------------------------
      real(rp) :: v, zeta, ene, ene2, xx, Fphi, Fpsi, sm, ce
      real(rp) :: nub, nue, nueprim, nub1, nue1, nueprim1, gg, d


      !---------------------------------------------------------------------------------
      ! Collisions
      !---------------------------------------------------------------------------------
      ene  = As/2.0_rp*vpi**2 + B*abs(mui)
      v    = sqrt(2.0_rp*abs(ene)/As)
      zeta = vpi / v

      xx   = sqrt(abs(ene)) * sqrt(Aeff) / sqrt(As)
      Fphi = erf(xx)
      Fpsi = (Fphi - 2.0_rp*xx*exp(-xx**2)/sqrt(pi)) / (2.0_rp*xx**2)

      gg       = (16.0_rp*R0*c0) / (vth**4 * (2.0_rp/As)**1.5_rp)
      nub1     = (gg/4.0_rp) / ene**1.5_rp * (Fphi - Fpsi)
      nue1     = gg / ene**1.5_rp * Fpsi
      nueprim1 = gg*sqrt(Aeff/As/pi) * exp(-ene*Aeff/As) / ene**2 - 2.5_rp*nue1/ene

      d  = delta
      ce = 1.0_rp / (1.0_rp + d*ene**(-1.0_rp)*(As/2.0_rp))

      nub     = nub1*ce
      nue     = nue1*ce
      nueprim = nueprim1*ce + nue1*2.0_rp*As*d / ((As*d + 2.0_rp*ene)**2.0_rp)

      sm   = 1.0_rp * sign(1.0_rp, sm1 - 0.5_rp)
      zeta = zeta*(1.0_rp - 1.0_rp*nub*dt_local) + sm*sqrt(dt_local*nub*(1.0_rp - zeta**2))

      sm   = 1.0_rp * sign(1.0_rp, sm2 - 0.5_rp)
      ene2 = ene - 1.0_rp*nue*dt_local*(ene - (1.5_rp + ene/(nue + 0.000001_rp)*nueprim)) + 2.0_rp*sm*sqrt(ene*nue*dt_local)
      ene  = abs(ene2)

      v     = sqrt(2.0_rp*abs(ene)/As)

      !---------------------------------------------------------------------------------
      ! Defaults
      !---------------------------------------------------------------------------------
      vcolx = 0.0_rp;  vcoly = 0.0_rp;  vcolz = 0.0_rp
      vcolp = (v*zeta - vpi)/dt_local
      vcolm = (ene/B*(1.0_rp - zeta**2) - mui)/(dt_local+0.00001_rp)

end subroutine collisions_1_MP

end module drift_kernels

