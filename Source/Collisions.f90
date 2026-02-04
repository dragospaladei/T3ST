    !      ! Subroutine :: Collisions
    !      ! Purpose    :: Collisions
    !      ! ??????????????????????????
    !      ! ??????
    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  01/01/2024   |     D. I. Palade       |   Remake
    !  ========================================================================================================================================================

    SUBROUTINE Collisions(dt, X, Y, Z, Vp, mut, vcolx, vcoly, vcolz, vcolm, vcolp, B, Bx, By, Bz)
       USE constants
       USE random_numbers
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN)        :: X, Y, Z, Vp, mut, B                              ! new basis coordinates (input)
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT)       :: vcolx, vcoly, vcolz, vcolm, vcolp

       ! Local variables
       REAL(KIND=dp), DIMENSION(Np)                       :: v, zeta, ene, alf, xx, Fphi, Fpsi     ! miscelaneous
       REAL(KIND=dp), DIMENSION(Np)                       :: Dpar, Dperp, hx, hy, hz, betx, bety, betz, Dparprim, Omega, nuf
       REAL(KIND=dp), DIMENSION(Np)                       :: Bx, By, Bz, Wx, Wy, Wz, Wp, Wm, DifR, Spp, Spm, Smp, Smm
       REAL(KIND=dp), DIMENSION(5, Np)                     :: grf!, grf1, grf2
       REAL(KIND=dp), DIMENSION(10*Np)                    :: grf1
       REAL(KIND=dp)                                       :: dt, ratio
       INTEGER                                               :: i

       logical :: contains_nan

       hx = 1.0
       hy = 1.0
       hz = X

       betx = Bx*hx/B
       bety = By*hy/B
       betz = Bz*hz/B
       Omega = B/Aw*wi*R0/vth

       ene = Aw/2.0*Vp**2 + B*abs(mut)
       v = sqrt(2.0*abs(ene)/Aw)
       zeta = Vp/v
       alf = ene/(v*B)
       xx = v/sqrt(2.0)/Sqrt(Aeff)
       Fphi = erf(xx)
       Fpsi = (Fphi - 2.0*xx*exp(-xx**2)/sqrt(pi))/2.0/xx**2

       Dpar = 4.0*R0*c0/vth**4*sqrt(Aeff/2.0)*Fpsi/xx
       Dperp = 2.0*R0*c0/vth**4*sqrt(Aeff/2.0)*(Fphi - Fpsi)
       Dparprim = -Dpar*((6.0*xx + 4.0*xx**3)*exp(-xx**2) - 3.0*Sqrt(Pi)*Erf(xx))/(xx*((2.0*xx)*exp(-xx**2) - Sqrt(Pi)*Erf(xx)))
       DifR = (Dpar*(1.0 - zeta**2) + Dperp*(1.0 + zeta**2))/2.0/Omega**2
       nuf = 4.0*R0*c0/vth**4*(sqrt(Aeff/2.0)**3)*Fpsi!!

       CALL PDF_G(5, Np, grf(1:5, :))           ! corr = 1 ~ (PDF_G, PDF_G, PDF_w_Bes)

       !        CALL RANDOM_NUMBER(grf1)
       !       grf(1,:) = sqrt(-2.0*log(grfaux(1,:)))*cos(2.0*pi*grfaux(2,:))
       !       grf(2,:) = sqrt(-2.0*log(grfaux(1,:)))*sin(2.0*pi*grfaux(2,:))
       !       grf(3,:) = sqrt(-2.0*log(grfaux(3,:)))*cos(2.0*pi*grfaux(4,:))
       !       grf(4,:) = sqrt(-2.0*log(grfaux(3,:)))*sin(2.0*pi*grfaux(4,:))
       !       grf(5,:) = sqrt(-2.0*log(grfaux(5,:)))*cos(2.0*pi*grfaux(6,:))

       Wm = grf(1, :)/sqrt(dt)!
       Wp = grf(2, :)/sqrt(dt)
       Wx = grf(3, :)/sqrt(dt)
       Wy = grf(4, :)/sqrt(dt)
       Wz = grf(5, :)/sqrt(dt)

       Spp = Sqrt(Dpar*zeta**2 + Dperp*(1.0 - zeta**2)) + 0.000001
       Spm = 0.0
       Smp = 2.0*Sqrt(2.0)*alf*zeta*(Dpar - Dperp)*(1.0 - zeta**2)/Spp
       Smm = 2.0*Sqrt(2.0)*alf*Sqrt(Dpar*Dperp*abs(1.0 - zeta**2))/Spp
       Spp = Sqrt(2.0)*Spp

       ratio = 1.0
       vcolp = ratio*(-nuf*Vp + zeta*(2.0*(Dpar - Dperp)/v + Dparprim))*1.0 + 1.0*(Spp*Wp + Spm*Wm)
       vcolm = ratio*(-2.0*nuf*mut + Aw*mut/ene*(v*Dparprim + 3.0*(Dpar - Dperp)) + 2.0*Aw*Dperp/B)*1.0 + 1.0*(Smp*Wp + Smm*Wm)
       vcolx = ratio*sqrt(2.0*DifR)*(Wx - Wx*betx*betx - Wy*bety*betx - Wz*betz*betx)/hx
       vcoly = ratio*sqrt(2.0*DifR)*(Wy - Wx*betx*bety - Wy*bety*bety - Wz*betz*bety)/hy
       vcolz = ratio*sqrt(2.0*DifR)*(Wz - Wx*betx*betz - Wy*bety*betz - Wz*betz*betz)/hz

       !       call print_error(X,Y,Z,Vp,mut,v,zeta,ene,B,xx,Fpsi,Fphi,Dpar,Dperp,Dparprim,DifR,vcolx,vcoly,vcolz,vcolm,vcolp,Spp,Spm,Smp,Smm)

    END SUBROUTINE Collisions

    SUBROUTINE Coll_Lorentz(dt, X, Y, Z, Vp, mut, vcolx, vcoly, vcolz, vcolm, vcolp, B)
       USE constants
       USE random_numbers
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN)        :: X, Y, Z, B                              ! new basis coordinates (input)
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT)       :: vcolx, vcoly, vcolz, vcolm, vcolp
       REAL(KIND=dp), DIMENSION(Np)                       :: Vp, mut, ene2

       ! Local variables
       REAL(KIND=dp), DIMENSION(Np)                       :: v, zeta, ene, xx, Fphi, Fpsi     ! miscelaneous
       REAL(KIND=dp), DIMENSION(Np)                       :: sm, ce, nub, nue, nueprim, nub1, nue1, nueprim1
       REAL(KIND=dp)                                       :: dt, ss, gg, d
       INTEGER                                               :: i, j, k

       logical :: contains_nan

       ene = Aw/2.0*Vp**2 + B*abs(mut) !good
       v = sqrt(2.0*abs(ene)/Aw)      !good
       zeta = Vp/v                       !good

       xx = sqrt(abs(ene))*sqrt(Aeff)/sqrt(Aw)                 ! v_24.03
       Fphi = erf(xx)                                           ! v_24.03
       Fpsi = (Fphi - 2.0*xx*exp(-xx**2)/sqrt(pi))/2.0/xx**2      ! v_24.03

       !       v_25.03
       gg = (16.0*R0*c0)/(vth**4*(2/Aw)**1.5)
       nub1 = (gg/4.0)/ene**1.5*(Fphi - Fpsi)
       nue1 = gg/ene**1.5*Fpsi
       nueprim1 = gg*sqrt(Aeff/Aw/pi)*exp(-ene*Aeff/Aw)/ene**2 - 2.5*nue1/ene

       d = delta! = 5.0*4./3./sqrt(pi)*(tmax-t0)/Nt*R0*2*logL*ndens*(Zeff*Zw*q0**2)**2/(pi*eps**2*mi**2*Aw**2*vth**4*sqrt(2.0/Aeff)**3)

       ce = 1.0/(1.0 + d*ene**(-1.0)*(Aw/2.0))
       nub = nub1*ce
       nue = nue1*ce
       nueprim = nueprim1*ce + nue1*2.0*Aw*d/((Aw*d + 2.0*ene)**2.0)

       call random_number(sm)
       sm = 1.0*sign(1.0, sm - 0.5)
       zeta = zeta*(1.0 - nub*dt) + sm*sqrt(dt*nub*(1.0 - zeta**2))    ! this line is crazy hard
       call random_number(sm)
       sm = 1.0*sign(1.0, sm - 0.5)
       ene2 = ene - 2.0*nue*dt*(ene - (1.5 + ene/(nue + 0.000001)*nueprim)) + 2.0*sm*sqrt(ene*nue*dt)
       ene = abs(ene2)

       v = sqrt(2.0*abs(ene)/Aw)
       Vp = v*zeta
       mut = ene/B*(1.0 - zeta**2)

       vcolx = 0.0
       vcoly = 0.0
       vcolz = 0.0
       vcolm = 0.0
       vcolp = 0.0

    END SUBROUTINE Coll_Lorentz

    SUBROUTINE Coll_Lorentz_drag(dt, X, Y, Z, Vp, mut, vcolx, vcoly, vcolz, vcolm, vcolp, B, Bx, By, Bz)
       USE constants
       USE random_numbers
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN)        :: X, Y, Z, B                              ! new basis coordinates (input)
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT)       :: vcolx, vcoly, vcolz, vcolm, vcolp
       REAL(KIND=dp), DIMENSION(Np)                       :: Vp, mut, ene2

       ! Local variables
       REAL(KIND=dp), DIMENSION(Np)                       :: v, zeta, ene, alf, xx, Fphi, Fpsi, se, vene, vzeta     ! miscelaneous
    REAL(KIND = dp),  DIMENSION(Np)                       :: hx, hy, hz, betx, bety, betz, Omega,sm, ce, nub, nue, nueprim, nub1, nue1, nueprim1
       REAL(KIND=dp), DIMENSION(Np)                       :: Bx, By, Bz
       REAL(KIND=dp)                                       :: dt, ratio, ss, gg, d
       INTEGER                                               :: i, j, k

       logical :: contains_nan

       ene = Aw/2.0*Vp**2 + B*abs(mut) !good
       v = sqrt(2.0*abs(ene)/Aw)      !good
       zeta = Vp/v                       !good

       xx = sqrt(abs(ene))*sqrt(Aeff)/sqrt(Aw)                 ! v_24.03
       Fphi = erf(xx)                                           ! v_24.03
       Fpsi = (Fphi - 2.0*xx*exp(-xx**2)/sqrt(pi))/2.0/xx**2      ! v_24.03

    !! v_25.03
       gg = (16.0*R0*c0)/(vth**4*(2/Aw)**1.5)
       nub1 = (gg/4.0)/ene**1.5*(Fphi - Fpsi)
       nue1 = gg/ene**1.5*Fpsi
       nueprim1 = gg*sqrt(Aeff/Aw/pi)*exp(-ene*Aeff/Aw)/ene**2 - 2.5*nue1/ene

       d = delta! = 5.0*4./3./sqrt(pi)*(tmax-t0)/Nt*R0*2*logL*ndens*(Zeff*Zw*q0**2)**2/(pi*eps**2*mi**2*Aw**2*vth**4*sqrt(2.0/Aeff)**3)

       ce = 1.0/(1.0 + d*ene**(-1.0)*(Aw/2.0))
       nub = nub1*ce
       nue = nue1*ce
       nueprim = nueprim1*ce + nue1*2.0*Aw*d/((Aw*d + 2.0*ene)**2.0)

       CALL random_number(sm)
       sm = 1.0*sign(1.0, sm - 0.5)
       zeta = zeta*(1.0 - nub*dt) + sm*sqrt(dt*nub*(1.0 - zeta**2))
       CALL random_number(sm)
       sm = 1.0*sign(1.0, sm - 0.5)
       ene2 = ene - 2.0*nue*dt*(ene - (1.5 + ene/(nue + 0.000001)*nueprim)) + 2.0*sm*sqrt(ene*nue*dt)
       ene = abs(ene2)
       v = sqrt(2.0*abs(ene)/Aw)

       vcolx = 0.0
       vcoly = 0.0
       vcolz = 0.0
       vcolm = (ene/B*(1.0 - zeta**2) - mut)/dt
       vcolp = (v*zeta - Vp)/dt

    END SUBROUTINE Coll_Lorentz_drag
