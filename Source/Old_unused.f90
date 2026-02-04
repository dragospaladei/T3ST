    !  ========================================================================================================================================================
    !      ! Subroutine :: drift
    !      ! Purpose    :: Computes drift velocities of gyro-centers projected onto the x,y,z,vp coordinate system
    !      ! x,y,z = R,R0\zeta,Z ; the ortogonal right-handed cylindrical coordinates
    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  --/--/2021   |     D. I. Palade       |   Initial version
    !      !  --/--/2023   |     L. M. Pomarjanschi |   Corections and development
    !      !  01/01/2024   |     D. I. Palade       |   Remake
    !  ========================================================================================================================================================

    SUBROUTINE Drift(dt, X, Y, Z, Vp, mu, q1, q2, q3, time, Vx, Vy, Vz, Ap, Vm, Ham, B, Vtx, Vty)
       USE constants
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN)    :: X, Y, Z, Vp, mu                               ! gyro-center coordinates (cylindrical)
       REAL(KIND=dp), INTENT(IN)    :: time, dt                                          ! the moment in time
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT)   :: Vx, Vy, Vz, Vm, Ap, Ham, q1, q2, q3, B, Vtx, Vty                           ! GC velocities

       ! Local variables
       REAL(KIND=dp), DIMENSION(Np)                    :: Esx, Esy, Esz, Bsx, Bsy, Bsz, Bsp, Etx, Ety, Etz             ! Effective (star, "s") electric/magnetic fields projected on (x,y,z)
       REAL(KIND=dp), DIMENSION(3, 3, Np)                :: F                                             ! The F_(ij)=<(dx_j x b).dx_i> matrix from the ExB term
       REAL(KIND=dp), DIMENSION(Np)                    :: vcolx, vcoly, vcolz, vcolm, vcolp

       ! ------------------------------------------------------------------------------------------------------------------------
       ! Estar/Bstar components in field-aligned coordinate system
       ! ------------------------------------------------------------------------------------------------------------------------

    CALL EBstar(dt,X,Y,Z,Vp,mu,q1,q2,q3,time,Esx,Esy,Esz,Bsx,Bsy,Bsz,Bsp,B,F,Ham,vcolx,vcoly,vcolz,vcolp,vcolm,Etx,Ety,Etz)                          ! The components of the Estar field
       ! ------------------------------------------------------------------------------------------------------------------------
       ! (Vx,Vy,Vz) = covariant components of the velocity :: Vx = dx/dt = dr/dt*grad(x)
       ! ------------------------------------------------------------------------------------------------------------------------
       Vx = Vp*Bsx/Bsp + rhoi/R0*(F(1, 1, :)*Esx + F(1, 2, :)*Esy + F(1, 3, :)*Esz)/Bsp             ! vp*Bs/Bsp + Es x b/Bs
       Vy = Vp*Bsy/Bsp + rhoi/R0*(F(2, 1, :)*Esx + F(2, 2, :)*Esy + F(2, 3, :)*Esz)/Bsp             ! vp*Bs/Bsp + Es x b/Bs
       Vz = Vp*Bsz/Bsp + rhoi/R0*(F(3, 1, :)*Esx + F(3, 2, :)*Esy + F(3, 3, :)*Esz)/Bsp             ! vp*Bs/Bsp + Es x b/Bs
       Ap = Zw/Aw*(Esx*Bsx + Esy*Bsy + Esz*Bsz)/Bsp                                           ! q/m*(Es*Bs)/Bsp

       Vtx = rhoi/R0*(F(1, 1, :)*Etx + F(1, 2, :)*Ety + F(1, 3, :)*Etz)/Bsp             ! vp*Bs/Bsp + Es x b/Bs
       Vty = rhoi/R0*(F(2, 1, :)*Etx + F(2, 2, :)*Ety + F(2, 3, :)*Etz)/Bsp             ! vp*Bs/Bsp + Es x b/Bs

!           Vtx = Etx
!          Vty = Ety

       Vx = Vx + vcolx
       Vy = Vy + vcoly
       Vz = Vz + vcolz
       Ap = Ap + vcolp
       Vm = vcolm

    END subroutine Drift
    !  ========================================================================================================================================================
    
        !      ! Subroutine :: EBstar
    !      ! Purpose    :: Computes the components of E,Bstar fields projected onto the (x,y,z) coordinates
    !      ! x,y,z = R,Z,\zeta ; the ortogonal right-handed cylindrical coordinates
    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  --/--/2021   |     D. I. Palade       |   Initial version
    !      !  --/--/2023   |     L. M. Pomarjanschi |   Corections and development
    !      !  01/01/2024   |     D. I. Palade       |   Remake
    !  ========================================================================================================================================================

    SUBROUTINE EBstar(dt,X,Y,Z,Vp,mu,qx,qy,qz,time,Esx,Esy,Esz,Bsx,Bsy,Bsz,Bsp,B,F,Ham,vcolx,vcoly,vcolz,vcolp,vcolm,Etx,Ety,Etz)
       USE constants
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN)     :: X, Y, Z, Vp, mu                                     ! GY coordinates (input)
       REAL(KIND=dp), INTENT(IN)     :: time                                            ! moment in time
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT)    :: Esx, Esy, Esz, Bsx, Bsy, Bsz, Bsp, Ham, Etx, Ety, Etz              ! Effective fields components (output)
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT)    :: vcolx, vcoly, vcolz, vcolm, vcolp, qx, qy, qz
       REAL(KIND=dp), DIMENSION(3, 3, Np), INTENT(OUT)    :: F                                                   ! F (exb) matrix

       ! Local variables
       REAL(KIND=dp), DIMENSION(Np)                    :: phi0, phix, phiy, phiz, phixt, phiyt, phizt         ! Lagrangian derivatives of the turbulent potential
       REAL(KIND=dp), DIMENSION(Np)                    :: phi0x, phi0y, phi0z, Tprofile, Armp                                ! The gradient of the neoclassical field Phi_0 (for quasineutrality in rotating plasmas)
       REAL(KIND=dp), DIMENSION(Np)                    :: gradBx, gradBy, gradBz, gradu2x, gradu2y, gradu2z   ! Grad B, grad u**2
       REAL(KIND=dp), DIMENSION(Np)                    :: rotbx, rotby, rotbz, rotux, rotuy, rotuz            ! curl(b), curl(u) components in the x,y,z system
       REAL(KIND=dp), DIMENSION(Np)                    :: Bx, By, Bz, B                                       ! Bx,  By, Bz, B, Bsp components in the x,y,z system
       REAL(KIND=dp), DIMENSION(3, 3, Np)                :: G, M                                                ! Matrices for fieldaligned-cylindrical (gradqi*dr/dxj), polarization
       REAL(KIND=dp), DIMENSION(Np)                    :: hx, hy, hz                                          ! lame coefficients
       REAL(KIND=dp)                                    :: dt                                             ! moment in time

       ! ------------------------------------------------------------------------------------------------------------------------
       ! Magnetic equilibrium quantities and turbulent field
       ! ------------------------------------------------------------------------------------------------------------------------

       CALL Equil(X, Y, Z, Bx, By, Bz, gradBx, gradBy, gradBz, gradu2x, gradu2y, gradu2z,&
           &rotbx, rotby, rotbz, rotux, rotuy, rotuz, phi0x, phi0y, phi0z, B, G, qx, qy, qz)
       CALL matrices2(X, Bx, By, Bz, G, F, M)

       hx = 1.0                                                       ! lame coefficient for x = R -> hx = 1
       hy = 1.0                                                       ! lame coefficient for y = Z -> hy = 1
       hz = X                                                         ! lame coefficient for z = zeta -> hz = R = X

       ! ------------------------------------------------------------------------------------------------------------------------
       ! Bstar covariant components in (x,y,z) coordinate system
       ! ------------------------------------------------------------------------------------------------------------------------
       Bsx = Bx + Aw/Zw*rhoi/R0*(Vp*rotbx + rotux)                    ! Bstar = B + m/q*[vp*rot(b)+ rot(u)]
       Bsy = By + Aw/Zw*rhoi/R0*(Vp*rotby + rotuy)
       Bsz = Bz + Aw/Zw*rhoi/R0*(Vp*rotbz + rotuz)

       ! ------------------------------------------------------------------------------------------------------------------------
       ! Bstar covariant components in (x,y,z) coordinate system from the RMP
       ! ------------------------------------------------------------------------------------------------------------------------
       IF (USE_magnturb == ON) THEN
          call RMPs(X, Y, Z, Armp)
          Bsx = Bsx + rhoi/R0*(Armp*rotbx)                    ! Bstar = B + m/q*[vp*rot(b)+ rot(u)]
          Bsy = Bsy + rhoi/R0*(Armp*rotby)
          Bsz = Bsz + rhoi/R0*(Armp*rotbz)
       END IF

       Bsp = (Bsx*Bx*hx**2 + Bsy*By*hy**2 + Bsz*Bz*hz**2)/B           ! Bsp = Bstar*b

       ! ------------------------------------------------------------------------------------------------------------------------
       ! Estar contravariant components in (x,y,z) coordinate system// the neoclassical part
       ! ------------------------------------------------------------------------------------------------------------------------
       Esx = -(phi0x + mu/Zw*gradBx - Aw/Zw*gradu2x)                 ! Estar = -del(phi0 + mu*B/q - m/q*u^2/2)
       Esy = -(phi0y + mu/Zw*gradBy - Aw/Zw*gradu2y)
       Esz = -(phi0z + mu/Zw*gradBz - Aw/Zw*gradu2z)

       ! ------------------------------------------------------------------------------------------------------------------------
       ! neoclassical electric field components  Etx, ...
       ! ------------------------------------------------------------------------------------------------------------------------
       Etx = Esx
       Ety = Esy
       Etz = Esz
       ! ------------------------------------------------------------------------------------------------------------------------
       ! Estar contravariant components in (x,y,z) coordinate system// the turbulent part
       ! ------------------------------------------------------------------------------------------------------------------------
       IF ((USE_turb == ON) .and. (time .gt. tt)) THEN

          IF (USE_real == OFF) THEN
             CALL fields_3D(qx, qy, qz, time, USE_larmor, mu, B, phi0, phix, phiy, phiz, phixt, phiyt, phizt)
          ELSEIF (USE_real == ON) THEN
             CALL fields_3D_single(qx, qy, qz, time, USE_larmor, mu, B, phi0, phix, phiy, phiz, phixt, phiyt, phizt)
          END IF

          Tprofile = 1.0*tanh(turbprof**2*(Y**2 + (X - 1.0)**2)/a0**2)     ! Phi = tanh(r/a)

          Esx = Esx - Tprofile*Phi*(G(1, 1, :)*phix + G(1, 2, :)*phiy + G(1, 3, :)*phiz)
          Esy = Esy - Tprofile*Phi*(G(2, 1, :)*phix + G(2, 2, :)*phiy + G(2, 3, :)*phiz)
          Esz = Esz - Tprofile*Phi*(G(3, 1, :)*phix + G(3, 2, :)*phiy + G(3, 3, :)*phiz)

          IF (USE_polar == ON) THEN
             Esx = Esx + Tprofile*Phi*Aw/Zw*rhoi/R0*(M(1, 1, :)*phixt + M(1, 2, :)*phiyt + M(1, 3, :)*phizt)/B
             Esy = Esy + Tprofile*Phi*Aw/Zw*rhoi/R0*(M(2, 1, :)*phixt + M(2, 2, :)*phiyt + M(2, 3, :)*phizt)/B
             Esz = Esz + Tprofile*Phi*Aw/Zw*rhoi/R0*(M(3, 1, :)*phixt + M(3, 2, :)*phiyt + M(3, 3, :)*phizt)/B
          END IF
       ELSE
          phi0 = 0.0
       END IF

       ! ------------------------------------------------------------------------------------------------------------------------
       ! Energy, in fact, Hamiltonian
       ! ------------------------------------------------------------------------------------------------------------------------

       Ham = Aw*Vp**2/2.0 + mu*B - (1.0 - Te/(Ti + Te))*Aw*(Omgt0**2*X**2)/2.0 + Zw*Phi*phi0

       ! ------------------------------------------------------------------------------------------------------------------------
       ! Collisions
       ! ------------------------------------------------------------------------------------------------------------------------
       IF ((USE_coll == ON) .and. (time .gt. tc)) THEN
          !         CALL Collisions(dt,X,Y,Z,Vp,mu,vcolx,vcoly,vcolz,vcolm,vcolp,B,Bx,By,Bz)
          CALL Coll_Lorentz(dt, X, Y, Z, Vp, mu, vcolx, vcoly, vcolz, vcolm, vcolp, B)
          !         CALL Coll_Lorentz_drag(dt,X,Y,Z,Vp,mu,vcolx,vcoly,vcolz,vcolm,vcolp,B,Bx,By,Bz)
       ELSE
          vcolx = 0.0; vcoly = 0.0; vcolz = 0.0; vcolp = 0.0; vcolm = 0.0
       END IF

       ! ------------------------------------------------------------------------------------------------------------------------
       ! Turbulent electric field components  Etx, ...
       ! ------------------------------------------------------------------------------------------------------------------------
       Etx = Esx - Etx
       Ety = Esy - Ety
       Etz = Esz - Etz

!          Etx = phix
!          Ety = phiy

    END SUBROUTINE EBstar
    ! =======================================================================================================================================
    
        !      ! Subroutine :: matrices
    !      ! Purpose    :: Computes the matrix elements for the ExB drift and for polarization components
    !      ! x,y,z = R,R0\zeta,Z ; the ortogonal right-handed cylindrical coordinates
    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  --/--/2021   |     D. I. Palade       |   Initial version
    !      !  --/--/2023   |     L. M. Pomarjanschi |   Corections and development
    !      !  01/01/2024   |     D. I. Palade       |   Remake
    !  ========================================================================================================================================================

    SUBROUTINE matrices(X, B, g, F, M)
       USE constants
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN)   :: X, B                                 ! positions, magnetic field
       REAL(KIND=dp), DIMENSION(3, 3, Np), INTENT(IN)   :: g                                    ! G_ij = grad(qj).dr/dxi;
       REAL(KIND=dp), DIMENSION(3, 3, Np), INTENT(OUT)  :: F, M                                 ! F_ij = (grad(xj)x b).grad(xi) ; M = (grad(xj)x b).dr/dxi;

       ! Local variables
       REAL(KIND=dp), DIMENSION(3, 3, Np)                  :: eta                                  ! auxiliary matrix
       REAL(KIND=dp), DIMENSION(Np)                      :: hx, hy, hz                           ! lame coefficients

       hx = 1.0                                                       ! lame coefficient for x = R -> hx = 1
       hy = 1.0                                                       ! lame coefficient for y = Z -> hy = 1
       hz = X                                                         ! lame coefficient for z = zeta -> hz = R = X

       !      B = sqrt(Bx**2*hx**2 + By**2*hy**2 + Bz**2*hz**2)                                 ! the norm of the magnetic field (with the covariant components)

       ! ------------------------------------------------------------------------------------------------------------------------
       ! The matrix elements F_ij = (grad(xj)x b).grad(xi);
       ! ------------------------------------------------------------------------------------------------------------------------
       F(1, 1, :) = 0.0
       F(2, 2, :) = 0.0
       F(3, 3, :) = 0.0
       F(1, 2, :) = (g(1, 2, :)*g(2, 1, :) - g(1, 1, :)*g(2, 2, :))/B/C1/C2/hx**2/hy**2
       F(1, 3, :) = (g(1, 2, :)*g(3, 1, :) - g(1, 1, :)*g(3, 2, :))/B/C1/C2/hx**2/hz**2
       F(2, 3, :) = (g(2, 2, :)*g(3, 1, :) - g(2, 1, :)*g(3, 2, :))/B/C1/C2/hy**2/hz**2
       F(2, 1, :) = -F(1, 2, :)
       F(3, 1, :) = -F(1, 3, :)
       F(3, 2, :) = -F(2, 3, :)
       !      F(1,2,:) = Bz/B*hz/(hx*hy)
       !      F(1,3,:) = - By/B*hy/(hx*hz)
       !      F(2,3,:) = Bx/B*hx/(hy*hz)

       ! ------------------------------------------------------------------------------------------------------------------------
       ! The matrix elements eta_ij = grad(qi).grad(qj);
       ! ------------------------------------------------------------------------------------------------------------------------
       eta(1, 1, :) = g(1, 1, :)*g(1, 1, :)/hx**2 + g(1, 2, :)*g(1, 2, :)/hy**2 + g(1, 3, :)*g(1, 3, :)/hz**2
       eta(1, 2, :) = g(1, 1, :)*g(2, 1, :)/hx**2 + g(1, 2, :)*g(2, 2, :)/hy**2 + g(1, 3, :)*g(2, 3, :)/hz**2
       eta(1, 3, :) = g(1, 1, :)*g(3, 1, :)/hx**2 + g(1, 2, :)*g(3, 2, :)/hy**2 + g(1, 3, :)*g(3, 3, :)/hz**2
       eta(2, 2, :) = g(2, 1, :)*g(2, 1, :)/hx**2 + g(2, 2, :)*g(2, 2, :)/hy**2 + g(2, 3, :)*g(2, 3, :)/hz**2
       eta(2, 3, :) = g(2, 1, :)*g(3, 1, :)/hx**2 + g(2, 2, :)*g(3, 2, :)/hy**2 + g(2, 3, :)*g(3, 3, :)/hz**2
       eta(3, 3, :) = g(3, 1, :)*g(3, 1, :)/hx**2 + g(3, 2, :)*g(3, 2, :)/hy**2 + g(3, 3, :)*g(3, 3, :)/hz**2
       eta(2, 1, :) = eta(1, 2, :)
       eta(3, 1, :) = eta(1, 3, :)
       eta(3, 2, :) = eta(2, 3, :)

       ! ------------------------------------------------------------------------------------------------------------------------
       ! The matrix elements M_ij = (grad(qj)x b).dr/dxi;
       ! ------------------------------------------------------------------------------------------------------------------------

       M(1, 1, :) = (g(1, 2, :)*eta(1, 1, :) - g(1, 1, :)*eta(1, 2, :))/B/C1/C2
       M(1, 2, :) = (g(1, 2, :)*eta(2, 1, :) - g(1, 1, :)*eta(2, 2, :))/B/C1/C2
       M(1, 3, :) = (g(1, 2, :)*eta(3, 1, :) - g(1, 1, :)*eta(3, 2, :))/B/C1/C2
       M(2, 1, :) = (g(2, 2, :)*eta(1, 1, :) - g(2, 1, :)*eta(1, 2, :))/B/C1/C2
       M(2, 2, :) = (g(2, 2, :)*eta(2, 1, :) - g(2, 1, :)*eta(2, 2, :))/B/C1/C2
       M(2, 3, :) = (g(2, 2, :)*eta(3, 1, :) - g(2, 1, :)*eta(3, 2, :))/B/C1/C2
       M(3, 1, :) = (g(3, 2, :)*eta(1, 1, :) - g(3, 1, :)*eta(1, 2, :))/B/C1/C2
       M(3, 2, :) = (g(3, 2, :)*eta(2, 1, :) - g(3, 1, :)*eta(2, 2, :))/B/C1/C2
       M(3, 3, :) = (g(3, 2, :)*eta(3, 1, :) - g(3, 1, :)*eta(3, 2, :))/B/C1/C2

    END SUBROUTINE matrices
    !  ========================================================================================================================================================
    
        !      ! Subroutine :: matrices
    !      ! Purpose    :: Computes the matrix elements for the ExB drift and for polarization components
    !      ! x,y,z = R,R0\zeta,Z ; the ortogonal right-handed cylindrical coordinates
    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  --/--/2021   |     D. I. Palade       |   Initial version
    !      !  --/--/2023   |     L. M. Pomarjanschi |   Corections and development
    !      !  01/01/2024   |     D. I. Palade       |   Remake
    !  ========================================================================================================================================================

    SUBROUTINE matrices2(X, Bx, By, Bz, G, F, M)
       USE constants
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN)   :: X, Bx, By, Bz                        ! positions, magnetic field
       REAL(KIND=dp), DIMENSION(3, 3, Np), INTENT(IN)   :: G                                    ! G_ij = grad(qj).dr/dxi;
       REAL(KIND=dp), DIMENSION(3, 3, Np), INTENT(OUT)  :: F, M                                 ! F_ij = (grad(xj)x b).grad(xi) ; M = (grad(xj)x b).dr/dxi;

       ! Local variables
       REAL(KIND=dp), DIMENSION(Np)                      :: hx, hy, hz, B                        ! lame coefficients

       hx = 1.0                                                       ! lame coefficient for x = R -> hx = 1
       hy = 1.0                                                       ! lame coefficient for y = Z -> hy = 1
       hz = X                                                         ! lame coefficient for z = zeta -> hz = R = X

       B = sqrt(Bx**2*hx**2 + By**2*hy**2 + Bz**2*hz**2)                                 ! the norm of the magnetic field (with the covariant components)

       ! ------------------------------------------------------------------------------------------------------------------------
       ! The matrix elements F_ij = (grad(xj)x b).grad(xi);
       ! ------------------------------------------------------------------------------------------------------------------------
       F(1, 1, :) = 0.0
       F(2, 2, :) = 0.0
       F(3, 3, :) = 0.0
       F(1, 2, :) = hz*Bz/(B*hx*hy)
       F(1, 3, :) = -hy*By/(B*hx*hz)
       F(2, 3, :) = hx*Bx/(B*hy*hz)
       F(2, 1, :) = -F(1, 2, :)
       F(3, 1, :) = -F(1, 3, :)
       F(3, 2, :) = -F(2, 3, :)

       ! ------------------------------------------------------------------------------------------------------------------------
       ! The matrix elements M_ij = (grad(qj)x b).dr/dxi;
       ! ------------------------------------------------------------------------------------------------------------------------
       if (USE_polar == ON) then
          M(1, 1, :) = hx*(Bz*hz*G(2, 1, :)/hy - By*hy*G(3, 1, :)/hz)/B
          M(1, 2, :) = hx*(Bz*hz*G(2, 2, :)/hy - By*hy*G(3, 2, :)/hz)/B
          M(1, 3, :) = hx*(Bz*hz*G(2, 3, :)/hy - By*hy*G(3, 3, :)/hz)/B

          M(2, 1, :) = hy*(Bx*hx*G(3, 1, :)/hz - Bz*hz*G(1, 1, :)/hx)/B
          M(2, 2, :) = hy*(Bx*hx*G(3, 2, :)/hz - Bz*hz*G(1, 2, :)/hx)/B
          M(2, 3, :) = hy*(Bx*hx*G(3, 3, :)/hz - Bz*hz*G(1, 3, :)/hx)/B

          M(3, 1, :) = hz*(By*hy*G(1, 1, :)/hx - Bx*hx*G(2, 1, :)/hy)/B
          M(3, 2, :) = hz*(By*hy*G(1, 2, :)/hx - Bx*hx*G(2, 2, :)/hy)/B
          M(3, 3, :) = hz*(By*hy*G(1, 3, :)/hx - Bx*hx*G(2, 3, :)/hy)/B

       else
          M = 0.0
       end if

    END SUBROUTINE matrices2
    !  ========================================================================================================================================================

    !  ========================================================================================================================================================
    !      ! Subroutine :: RK4
    !      ! Purpose    :: Runge-Kutta-4 method for time evolution of x'(t)=f(x(t),t)
    !      ! The case of 4 variables : X,Y,Z,Vp is considered with velocities (Vx,Vy,Vz,Ap)
    !      ! Record of revisions:
    !      !    Date       |    Programmer       |   Description of change
    !      !  =======      |   =============     |   =========================
    !      !  --/--/2021   |     D. I. Palade    |   Initial version
    !      !  05/09/2022   |     L.Pomarjanschi  |   Corrections and development
    !      !  03/12/2023   |     D.Palade        |   Corrections and development
    !  ========================================================================================================================================================

    SUBROUTINE RK4(X, Y, Z, Vp, mu, pb, t01, tmax1, Xtot,Ytot,Ztot,Vptot,mutot,Htot,q1tot,q2tot,q3tot, Xtr,Ytr,Ztr,Vptr,mutr,Htr,q1tr,q2tr,q3tr,vit,dif,Vcorff,VcorTT,VcorTN,VcorNT) ! time-propagation of test-particles
       USE constants
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), INTENT(IN)     :: t01, tmax1                  ! initial time, final time
       REAL(KIND=dp), DIMENSION(Np, 2), INTENT(OUT)    :: Xtot, Ytot, Ztot, Vptot, mutot, Htot, q1tot, q2tot, q3tot! store gyro-center coordinates at each moment in time
       REAL(KIND=dp), DIMENSION(ntraj, Nt + 1), INTENT(OUT)    :: Xtr, Ytr, Ztr, Vptr, mutr, Htr, q1tr, q2tr, q3tr! store gyro-center coordinates at each moment in time
       REAL(KIND=dp), DIMENSION(4, Nt + 1), INTENT(OUT)    :: vit                         ! transport coefficients, pinch
       REAL(KIND=dp), DIMENSION(4, Nt + 1), INTENT(OUT)    :: dif                         ! transport coefficients, diffusion

       REAL(KIND=dp), DIMENSION(Np, Nt + 1)                 :: Vrr, VrT                    ! transport coefficients, diffusion
       REAL(KIND=dp), DIMENSION(Nt + 1, Nt + 1)               :: Vcorff, VcorTT, VcorTN, VcorNT !                     ! transport coefficients, diffusion

       ! Local variables
       REAL(KIND=dp), DIMENSION(Np)                        :: X, Y, Z, Vp, mu             ! gyro-center coordinates
       REAL(KIND=dp), DIMENSION(Np)                        :: q1, q2, q3, BB, pb          ! gyro-center coordinates
       REAL(KIND=dp), DIMENSION(Np)                        :: Vx, Vy, Vz, Ap, Vm          ! gyro-center velocities
       REAL(KIND=dp), DIMENSION(Np)                        :: Wx, Wy, Wz, Wp, Wm          ! auxiliary gyro-center velocities
       REAL(KIND=dp), DIMENSION(Np)                        :: Ham, Vrinit, VrinitT, Vtx, Vty ! store gyro-center coordinates at each moment in time
       REAL(KIND=dp), DIMENSION(Nt + 1)                      :: t                           ! discrete time vector
       REAL(KIND=dp)                                                    :: dt                          ! time step
       INTEGER                                                               :: k, i                        ! auxiliary integer for time steps
       REAL(KIND=dp), DIMENSION(Np)                      :: mask, maskR, maskZ, rrrr, ttt   !, q2tot, q3tot,Btot  ! store GC coordinates

       dt = (tmax1 - t01)/REAL(Nt)                           ! time step
       t = [(i*dt + t01, i=0, Nt)]                             ! discrete time vector

       !     dt = -dt

       Xtot(:, 1) = X                                         ! store GC (X) at t=0
       Ytot(:, 1) = Y                                         ! store GC (Y) at t=0
       Ztot(:, 1) = Z                                         ! store GC (Z) at t=0
       Vptot(:, 1) = Vp                                       ! store GC (Vp) at t=0
       mutot(:, 1) = mu                                       ! store GC (mu) at t=0
    CALL Drift(dt, X, Y, Z, Vp, mu, q1, q2, q3, t(1), Vx, Vy, Vz, Ap, Vm, Ham,BB, Vtx, Vty)                                                   ! velocities
       Htot(:, 1) = Ham                                       ! store GC (E) at t=0
       q1tot(:, 1) = q1                                       ! store GC (x) at t=0
       q2tot(:, 1) = q2                                       ! store GC (y) at t=0
       q3tot(:, 1) = q3                                       ! store GC (z) at t=0
!     q1tot(:,1) = (Vx*(X-1.0) + Vy*Y)/sqrt((X-1.0)**2 + Y**2)
       !   q1tot(:,1) = Vty       !*                                  ! store GC (x) at t=0

       Vrinit = (Vx*(X - 1.0) + Vy*Y)/sqrt((X - 1.0)**2 + Y**2)  ! V_x(t=0)
       VrinitT = (Vtx*(X - 1.0) + Vty*Y)/sqrt((X - 1.0)**2 + Y**2)! V_x^turb(t=0)

       DO k = 1, Nt + 1                                          ! start of the time loop

          CALL Drift(dt, X, Y, Z, Vp, mu, q1, q2, q3, t(k), Vx, Vy, Vz, Ap, Vm, Ham, BB, Vtx, Vty)     ! velocities
          ! because this call of "drift" offers us the evaluation of quantities at non-intermediate times, we proceed to save/compute Lagrangian quant
          !      CALL Diagnostic(k, X,Y,Z,Vp,mu,Ham,q1,q2,q3,Vx,Vy, Vtx,Vty,Vrinit,VrinitT,Xtr(:,k),Ytr(:,k),Ztr(:,k),Vptr(:,k),mutr(:,k),Htr(:,k),q1tr(:,k),q2tr(:,k),q3tr(:,k),Vrr(:,k),VrT(:,k),vit,dif)

          ! = = = = = = = = = = = = = = =
          Xtr(:, k) = X(1:ntraj)                               ! store GC (X) at each "k"-th moment
          Ytr(:, k) = Y(1:ntraj)                               ! store GC (Y) at each "k"-th moment
          Ztr(:, k) = Z(1:ntraj)                               ! store GC (Z) at each "k"-th moment
          Vptr(:, k) = Vp(1:ntraj)                             ! store GC (Vp)at each "k"-th moment
          mutr(:, k) = mu(1:ntraj)                             ! store GC (Vp)at each "k"-th moment
          Htr(:, k) = Ham(1:ntraj)                             ! store GC (X) at each "k"-th moment
          q1tr(:, k) = q1(1:ntraj)                             ! store GC (Vp)at each "k"-th moment
          q2tr(:, k) = q2(1:ntraj)                             ! store GC (Vp)at each "k"-th moment
          q3tr(:, k) = q3(1:ntraj)                             ! store GC (X) at each "k"-th moment

          Vrr(:, k) = (Vx*(X - 1.0) + Vy*Y)/sqrt((X - 1.0)**2 + Y**2)
          VrT(:, k) = (Vtx*(X - 1.0) + Vty*Y)/sqrt((X - 1.0)**2 + Y**2)

          !       Vrr(:,k) = Vtx    !*
          !       VrT(:,k) = Vty  !*
!          q1 = Vtx

          if ((magnetic_model .eq. 1) .or. (magnetic_model .eq. 2)) then
             maskR = (X - minR)*(maxR - X)
             maskZ = (Y - minZ)*(maxZ - Y)
             mask = merge(1.0, 0.0, maskR >= 0.0)*merge(1.0, 0.0, maskZ >= 0.0)
             mask = mask*real(Np)/sum(mask)
          else
             mask = 1.0
          end if

          vit(1, k) = sum(pb*mask*q1)/Np
          dif(1, k) = sum(pb*mask*q1**2)/Np
          vit(2, k) = sum(pb*mask*VrT(:, k))/Np                   ! turbulent drift
          dif(2, k) = sum(pb*mask*VrT(:, k)*VrinitT)/Np           ! turbulent L correlation
          vit(3, k) = sum(pb*mask*(Vrr(:, k) - VrT(:, k)))/Np
          dif(3, k) = sum(pb*mask*VrT(:, k)*(Vrinit - VrinitT))/Np
          vit(4, k) = sum(pb*mask*(Vrr(:, k) - VrT(:, k))*VrinitT)/Np
          dif(4, k) = sum(pb*mask*(Vrr(:, k) - VrT(:, k))*(Vrinit - VrinitT))/Np
          ! = = = = = = = = = = = = = = =

          Wx = Vx/6.0                                                             ! intermediate velocities between RK4 steps
          Wy = Vy/6.0                                                             ! intermediate velocities between RK4 steps
          Wz = Vz/6.0                                                             ! intermediate velocities between RK4 steps
          Wm = Vm/6.0                                                             ! intermediate velocities between RK4 steps
          Wp = Ap/6.0                                                             ! intermediate velocities between RK4 steps

        CALL Drift(dt, X+dt*Vx/2.0,Y+dt*Vy/2.0,Z+dt*Vz/2.0,Vp+dt*Ap/2.0, mu+dt*Vm/2.0, q1, q2, q3,t(k)+dt/2.0,Vx, Vy, Vz, Ap, Vm, Ham, BB, Vtx, Vty)        ! velocities
          Wx = Wx + Vx/3.0                                                        ! intermediate velocities between RK4 steps
          Wy = Wy + Vy/3.0                                                        ! intermediate velocities between RK4 steps
          Wz = Wz + Vz/3.0                                                        ! intermediate velocities between RK4 steps
          Wm = Wm + Vm/3.0                                                        ! intermediate velocities between RK4 steps
          Wp = Wp + Ap/3.0                                                        ! intermediate velocities between RK4 steps

        CALL Drift(dt, X+dt*Vx/2.0,Y+dt*Vy/2.0,Z+dt*Vz/2.0,Vp+dt*Ap/2.0, mu+dt*Vm/2.0, q1, q2, q3,t(k)+dt/2.0,Vx, Vy, Vz, Ap, Vm, Ham, BB, Vtx, Vty)        ! velocities
          Wx = Wx + Vx/3.0                                                        ! intermediate velocities between RK4 steps
          Wy = Wy + Vy/3.0                                                        ! intermediate velocities between RK4 steps
          Wz = Wz + Vz/3.0                                                        ! intermediate velocities between RK4 steps
          Wm = Wm + Vm/3.0                                                        ! intermediate velocities between RK4 steps
          Wp = Wp + Ap/3.0                                                        ! intermediate velocities between RK4 steps

        CALL Drift(dt, X+dt*Vx,Y+dt*Vy,Z+dt*Vz,Vp+dt*Ap, mu+dt*Vm, q1, q2, q3,t(k)+dt,Vx, Vy, Vz, Ap, Vm, Ham, BB, Vtx, Vty)                            ! velocities
          Wx = Wx + Vx/6.0                                                        ! intermediate velocities between RK4 steps
          Wy = Wy + Vy/6.0                                                        ! intermediate velocities between RK4 steps
          Wz = Wz + Vz/6.0                                                        ! intermediate velocities between RK4 steps
          Wm = Wm + Vm/6.0                                                        ! intermediate velocities between RK4 steps
          Wp = Wp + Ap/6.0                                                        ! intermediate velocities between RK4 steps

          X = X + dt*Wx                                                          ! (k+1)-th  positions on X
          Y = Y + dt*Wy                                                          ! (k+1)-th  positions on Y
          Z = Z + dt*Wz                                                          ! (k+1)-th  positions on Z
          Vp = Vp + dt*Wp                                                         ! (k+1)-th  positions on Vp
          mu = mu + dt*Wm                                                         ! (k+1)-th  positions on Vp
       END DO

       Xtot(:, 2) = X                                                           ! store GC (X) at each "k"-th moment
       Ytot(:, 2) = Y                                                           ! store GC (Y) at each "k"-th moment
       Ztot(:, 2) = Z                                                           ! store GC (Z) at each "k"-th moment
       Vptot(:, 2) = Vp                                                         ! store GC (Vp)at each "k"-th moment
       mutot(:, 2) = mu                                                         ! store GC (Vp)at each "k"-th moment
       Htot(:, 2) = Ham                                                         ! store GC (X) at each "k"-th moment
       q1tot(:, 2) = q1                                       ! store GC (x) at t=0
       q2tot(:, 2) = q2                                       ! store GC (y) at t=0
       q3tot(:, 2) = q3                                       ! store GC (z) at t=0
!     q1tot(:,2) = (Vx*(X-1.0) + Vy*Y)/sqrt((X-1.0)**2 + Y**2)
       !   q1tot(:,2) = Vty           !*                              ! store GC (x) at t=0

       If (USE_corr .eq. 1) then

          !$omp parallel do private(i,k) shared(Vcorff,Vrr,VrT,Np,pb)! default(none)
          do i = 1, Nt + 1
             do k = 1, Nt + 1
                Vcorff(i, k) = sum(pb*Vrr(:, i)*Vrr(:, k))/Np - sum(pb*Vrr(:, i))*sum(pb*Vrr(:, k))/(Np*Np)
                VcorTT(i, k) = sum(pb*VrT(:, i)*VrT(:, k))/Np - sum(pb*VrT(:, i))*sum(pb*VrT(:, k))/(Np*Np)
             VcorNT(i, k) = sum(pb*(Vrr(:, i) - VrT(:, i))*VrT(:, k))/Np - sum(pb*(Vrr(:, i) - VrT(:, i)))*sum(pb*VrT(:, k))/(Np*Np)
             VcorTN(i, k) = sum(pb*VrT(:, i)*(Vrr(:, k) - VrT(:, k)))/Np - sum(pb*VrT(:, i))*sum(pb*(Vrr(:, k) - VrT(:, k)))/(Np*Np)
             end do
          end do
          !$omp end parallel do

       end if
    END SUBROUTINE RK4
    
    !  ========================================================================================================================================================
    !   This file contains subroutines dedicated to the evaluation of GRF at multiple space& time locations
    !            First and second order derivatives of fields are also implemented
    !  1.  fields_3D :: compute GRF values (and derivatives) on a set of space-time locations in 3D
    !  2.  fields_2D :: compute GRF values (and derivatives) on a set of space-time locations in 2D
    !  3.  potent_3D :: compute GRF values on a set of space-time locations in 3D
    !  4.  potent_2D :: compute GRF values on a set of space-time locations in 2D
    !  5.  colliz    :: compute GRF values for two only-time dependent GRFs;
    !========================================================================================================
    !========================================================================================================
    ! Subroutine :: fields_3D
    ! Purpose    :: compute GRF values (and derivatives) on a set of space-time locations
    ! Definition of input : X,Y,Z = particle space locations ; time = time of evaluation
    !                     : Nc = the no of partial waves used ;; Np = the no. of particles(space locations)
    !                     : kx,ky,kz = wavenumbers ;; ph = phases ;; w = frequencies
    ! Definition of output: phi0 = the values of the field ;; phix,phiy,phiz = the values of space derivatives ;; phixt,phiyt,phizt=the space-time derivatives
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  01/10/2022   |     L.Pomarjanschi  |   Testing

    SUBROUTINE RMPs(X, Y, Z, Armp)
       USE constants
       USE omp_lib
       USE random_numbers
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN) :: X, Y, Z                   ! positions
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT):: Armp!, rotbx, rotby, rotbz

       ! Local variables
       INTEGER                         :: i
       REAL(KIND=dp), DIMENSION(Nc)                :: zintc, zints               ! auxiliary
       REAL(KIND=dp), DIMENSION(Np)                :: rr, theta, zzeta

       !      rr = sqrt((X - 1.0)**2 + Y**2) + 0.00001
       !      theta = atan2(Y, X-1.0)      ! atan(z,rr)
       !      zzeta = Z

       !      !$OMP PARALLEL DO &
       !      !$OMP SHARED (rr, theta, zzeta) &
       !      !$OMP PRIVATE (i,zintc)
       !      DO i=1,Np
       !       zintc = exp(-((rr - rrmp(:,i))/sigm)**4)*cos(Mrmp(:,i)*theta(i) + Nrmp(:,i)*zzeta + phrmp(:,i))     ! sin(k_i x+ph_i) terms for all realizations (i) and all wavenumbers
       !       phi0(i)  = sum(zints)
       !      END DO
       !      !$OMP END parallel DO

       Armp = 0.0
       !      rotbx = 0.0
       !      rotby = 0.0
       !      rotbz = 0.0
    END SUBROUTINE RMPs

    !========================================================================================================
    ! Subroutine :: fields_3D
    ! Purpose    :: compute GRF values (and derivatives) on a set of space-time locations
    ! Definition of input : X,Y,Z = particle space locations ; time = time of evaluation
    !                     : Nc = the no of partial waves used ;; Np = the no. of particles(space locations)
    !                     : kx,ky,kz = wavenumbers ;; ph = phases ;; w = frequencies
    ! Definition of output: phi0 = the values of the field ;; phix,phiy,phiz = the values of space derivatives ;; phixt,phiyt,phizt=the space-time derivatives
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  01/10/2022   |     L.Pomarjanschi  |   Testing

    SUBROUTINE fields_3D(X, Y, Z, time, ind, mut, B, phi0, phix, phiy, phiz, phixt, phiyt, phizt)
       USE constants
       USE omp_lib
       USE random_numbers
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN) :: X, Y, Z                   ! positions
       REAL(KIND=dp), INTENT(IN) :: time                    ! time
       REAL(KIND=dp), INTENT(IN) :: ind                     ! index for Larmor radius;; 1. -> use Larmor radius;; anything else -> don't use Larmor radius
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN) :: mut, B                  ! magnetic moment and |B|

       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT):: phi0, phix, phiy, phiz, phixt, phiyt, phizt

       ! Local variables
       INTEGER                         :: i
       REAL(KIND=dp), DIMENSION(Nc)             :: zintc, zints, faza              ! auxiliary
       !      REAL(KIND=dp), DIMENSION(Nc,Np)          :: L, L1                          ! auxiliary for Larmor radius effects

       !   ky = BESSEL_J0(kperp*SPREAD(SQRT(2.0*Aw*abs(mu)/Zw**2/B), DIM = 1, NCOPIES = Nc))
       !   not fast enough; a little bit faster than series do, but sloeer than parallel
       ! i have inverted the indexes to get it 3 times faster (larmor)
       !IF(ind==1.) THEN
       !        !$OMP PARALLEL DO &
       !        !$OMP SHARED (mut,B) &
       !        !$OMP PRIVATE (i)
       !        DO i=1, Np
       !          L(:,i) = BESSEL_J0(kperp(:,i)*SQRT(2.0*Aw*abs(mut(i))/Zw**2/B(i)))
       !        END DO
       !        !$OMP END parallel DO
       !      ELSE
       !        L = 1.
       !      END IF

       !      call omp_set_num_threads(30)

       !$OMP PARALLEL DO &
       !$OMP SHARED (X,Y,Z,time) &
       !$OMP PRIVATE (i,zintc,zints,faza) schedule(dynamic)
       DO i = 1, Np
          faza = kx(:, i)*X(i) + ky(:, i)*Y(i) + kz(:, i)*Z(i) - w(:, i)*time + ph(:, i)
          zintc = cos(faza)     ! sin(k_i x+ph_i) terms for all realizations (i) and all wavenumbers
          zints = sin(faza)
          !zintc = cos(kx(:,i)*X(i)+ky(:,i)*Y(i)+kz(:,i)*Z(i)-w(:,i)*time+ph(:,i))
          !zints = sin(kx(:,i)*X(i)+ky(:,i)*Y(i)+kz(:,i)*Z(i)-w(:,i)*time+ph(:,i))
          phi0(i) = DOT_PRODUCT(L(:, i), zints)                                        ! summing over all waves in a realization
          phix(i) = DOT_PRODUCT(L(:, i)*kx(:, i), zintc)
          phiy(i) = DOT_PRODUCT(L(:, i)*ky(:, i), zintc)
          phiz(i) = DOT_PRODUCT(L(:, i)*kz(:, i), zintc)
          phixt(i) = DOT_PRODUCT(L(:, i)*kx(:, i)*w(:, i), zints)
          phiyt(i) = DOT_PRODUCT(L(:, i)*ky(:, i)*w(:, i), zints)
          phizt(i) = DOT_PRODUCT(L(:, i)*kz(:, i)*w(:, i), zints)
       END DO
       !$OMP END parallel DO

       phi0 = phi0/sqrt(Nc/2.0)                            ! normalization
       phix = phix/sqrt(Nc/2.0)                            ! normalization
       phiy = phiy/sqrt(Nc/2.0)
       phiz = phiz/sqrt(Nc/2.0)
       phixt = phixt/sqrt(Nc/2.0)
       phiyt = phiyt/sqrt(Nc/2.0)
       phizt = phizt/sqrt(Nc/2.0)

    END SUBROUTINE fields_3D

    !========================================================================================================
    ! Subroutine :: fields_3D
    ! Purpose    :: compute GRF values (and derivatives) on a set of space-time locations
    ! Definition of input : X,Y,Z = particle space locations ; time = time of evaluation
    !                     : Nc = the no of partial waves used ;; Np = the no. of particles(space locations)
    !                     : kx,ky,kz = wavenumbers ;; ph = phases ;; w = frequencies
    ! Definition of output: phi0 = the values of the field ;; phix,phiy,phiz = the values of space derivatives ;; phixt,phiyt,phizt=the space-time derivatives
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  01/10/2022   |     L.Pomarjanschi  |   Testing

    SUBROUTINE fields_3D_single(X, Y, Z, time, ind, mut, B, phi0, phix, phiy, phiz, phixt, phiyt, phizt)
       USE constants
       USE omp_lib
       USE random_numbers
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN) :: X, Y, Z                   ! positions
       REAL(KIND=dp), INTENT(IN) :: time                    ! time
       REAL(KIND=dp), INTENT(IN) :: ind                     ! index for Larmor radius;; 1. -> use Larmor radius;; anything else -> don't use Larmor radius
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN) :: mut, B                  ! magnetic moment and |B|

       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT):: phi0, phix, phiy, phiz, phixt, phiyt, phizt

       ! Local variables
       INTEGER                         :: i
       REAL(KIND=dp), DIMENSION(Nc)             :: zintc, zints, Ls2, faza!, kxs,kys,kzs,ws,phs,Ls,faza              ! auxiliary

       !$OMP PARALLEL DO &
       !$OMP SHARED (X,Y,Z,time) &
       !$OMP PRIVATE (i,zintc,zints,faza,Ls2) schedule(dynamic)
       DO i = 1, Np
          Ls2 = Ls(:, i)
          faza = kxs*X(i) + kys*Y(i) + kzs*Z(i) - ws*time + phs
          zintc = cos(faza)     ! sin(k_i x+ph_i) terms for all realizations (i) and all wavenumbers
          zints = sin(faza)
          phi0(i) = DOT_PRODUCT(Ls2, zints)                                        ! summing over all waves in a realization
          phix(i) = DOT_PRODUCT(Ls2*kxs, zintc)
          phiy(i) = DOT_PRODUCT(Ls2*kys, zintc)
          phiz(i) = DOT_PRODUCT(Ls2*kzs, zintc)
          phixt(i) = DOT_PRODUCT(Ls2*kxs*ws, zints)
          phiyt(i) = DOT_PRODUCT(Ls2*kys*ws, zints)
          phizt(i) = DOT_PRODUCT(Ls2*kzs*ws, zints)
       END DO
       !$OMP END parallel DO

       phi0 = phi0/sqrt(Nc/2.0)                            ! normalization
       phix = phix/sqrt(Nc/2.0)                            ! normalization
       phiy = phiy/sqrt(Nc/2.0)
       phiz = phiz/sqrt(Nc/2.0)
       phixt = phixt/sqrt(Nc/2.0)
       phiyt = phiyt/sqrt(Nc/2.0)
       phizt = phizt/sqrt(Nc/2.0)

    END SUBROUTINE fields_3D_single

    !========================================================================================================
    ! Subroutine :: fields_2D
    ! Purpose    :: compute GRF values (and derivatives) on a set of space-time locations
    ! Definition of input : X,Y = particle space locations ; time = time of evaluation
    !                     : Nc = the no of partial waves used ;; Np = the no. of particles(space locations)
    !                     : kx,ky = wavenumbers ;; ph = phases ;; w = frequencies
    ! Definition of output: phi0 = the values of the field ;; phix,phiy = the values of space derivatives ;; phixt,phiyt=the space-time derivatives
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  01/10/2022   |     L.Pomarjanschi  |   Testing

    SUBROUTINE fields_2D(X, Y, time, phi0, phix, phiy, phixt, phiyt)
       USE constants
       USE random_numbers
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN) :: X, Y
       REAL(KIND=dp), INTENT(IN) :: time

       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT):: phi0, phix, phiy, phixt, phiyt

       ! Local variables
       INTEGER                         :: i
       REAL(KIND=dp), DIMENSION(Nc)             :: zintc, zints

       !$OMP PARALLEL DO &
       !$OMP SHARED (X,Y,time) &
       !$OMP PRIVATE (i,zintc,zints)
       DO i = 1, Np
          zintc = cos(kx(:, i)*X(i) + ky(:, i)*Y(i) - w(:, i)*time + ph(:, i))     ! sin(k_i x+ph_i) terms for all realizations (i) and all wavenumbers
          zints = sin(kx(:, i)*X(i) + ky(:, i)*Y(i) - w(:, i)*time + ph(:, i))
          phi0(i) = sum(zints)                                                        ! summing over all waves in a realization
          phix(i) = DOT_PRODUCT(kx(:, i), zintc)
          phiy(i) = DOT_PRODUCT(ky(:, i), zintc)
          phixt(i) = DOT_PRODUCT(kx(:, i)*w(:, i), zints)
          phiyt(i) = DOT_PRODUCT(ky(:, i)*w(:, i), zints)
       END DO
       !$OMP END parallel DO

       phi0 = phi0/sqrt(Nc/2.0)                            ! normalization
       phix = phix/sqrt(Nc/2.0)                            ! normalization
       phiy = phiy/sqrt(Nc/2.0)
       phixt = phixt/sqrt(Nc/2.0)
       phiyt = phiyt/sqrt(Nc/2.0)

    END SUBROUTINE fields_2D

    !========================================================================================================
    !========================================================================================================
    ! Subroutine :: potent_3D
    ! Purpose    :: compute only the GRF values on a set of space-time locations
    ! Definition of input : X,Y = particle space locations ; time = time of evaluation
    !                     : Nc = the no of partial waves used ;; Np = the no. of particles(space locations)
    !                     : kx,ky = wavenumbers ;; ph = phases ;; w = frequencies
    ! Definition of output: phi0 = the values of the field
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  01/10/2022   |     L.Pomarjanschi  |   Testing

    SUBROUTINE potent_3D(X, Y, Z, time, phi0)
       USE constants
       USE random_numbers
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN) :: X, Y, Z
       REAL(KIND=dp), INTENT(IN) :: time

       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT):: phi0

       ! Local variables
       INTEGER                         :: i

       !$OMP PARALLEL DO &
       !$OMP SHARED (X,Y,Z,time) &
       !$OMP PRIVATE (i)
       DO i = 1, Np
          phi0(i) = sum(sin(kx(:, i)*X(i) + ky(:, i)*Y(i) + kz(:, i)*Z(i) - w(:, i)*time + ph(:, i))) ! sin(k_i x+ph_i) terms for all realizations (i) and all wavenumbers; summing over all waves in a realization
       END DO
       !$OMP END parallel DO

       phi0 = phi0/sqrt(Nc/2.0)                            ! normalization

    END SUBROUTINE potent_3D

    !========================================================================================================
    ! Subroutine :: potent_2D
    ! Purpose    :: compute only the GRF values on a set of space-time locations
    ! Definition of input : X,Y = particle space locations ; time = time of evaluation
    !                     : Nc = the no of partial waves used ;; Np = the no. of particles(space locations)
    !                     : kx,ky = wavenumbers ;; ph = phases ;; w = frequencies
    ! Definition of output: phi0 = the values of the field
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  01/10/2022   |     L.Pomarjanschi  |   Testing

    SUBROUTINE potent_2D(X, Y, time, phi0)
       USE constants
       USE random_numbers
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN) :: X, Y
       REAL(KIND=dp), INTENT(IN) :: time

       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT):: phi0

       ! Local variables
       INTEGER                         :: i

       !$OMP PARALLEL DO &
       !$OMP SHARED (X,Y,time) &
       !$OMP PRIVATE (i)
       DO i = 1, Np
          phi0(i) = sum(sin(kx(:, i)*X(i) + ky(:, i)*Y(i) - w(:, i)*time + ph(:, i)))     ! sin(k_i x+ph_i) terms for all realizations (i) and all wavenumbers, summing over all waves in a realization
       END DO
       !$OMP END parallel DO

       phi0 = phi0/sqrt(Nc/2.0)                            ! normalization

    END SUBROUTINE potent_2D

    !========================================================================================================
    ! Subroutine :: colliz
    ! Purpose    :: computes the values of two GRFs at the same locations; the fields are uncorrelated
    !               The name stems from its utility as generator of collisional terms
    ! Definition of input : X,Y = particle space locations ; time = time of evaluation
    !                     : Nc = the no of partial waves used ;; Np = the no. of particles(space locations)
    !                     : kx,ky = wavenumbers ;; ph = phases ;; w = frequencies
    ! Definition of output: phi0 = the values of the field
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  01/10/2022   |     L.Pomarjanschi  |   Testing

    SUBROUTINE colliz(time, Cx, Cy)
       USE constants
       USE random_numbers
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), INTENT(IN) :: time
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT):: Cx, Cy

       ! Local variables
       INTEGER                         :: i

       !$OMP PARALLEL DO &
       !$OMP SHARED (time) &
       !$OMP PRIVATE (i)
       DO i = 1, Np
          Cx(i) = sum(sin(wc(:, i)*time + phx(:, i)))     ! sin(w*t+ph) terms for all realizations (i) and all wavenumbers
          Cy(i) = sum(sin(wc(:, i)*time + phy(:, i)))
       END DO
       !$OMP END parallel DO

       Cx = Cx/sqrt(Nc/2.0)                         ! normalization ;;  the Pc factor is related to the strength of such fields;
       Cy = Cy/sqrt(Nc/2.0)

    END SUBROUTINE colliz
