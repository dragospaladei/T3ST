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
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN)   :: X, Y                  ! new basis coordinates (input)
       REAL(KIND=dp), DIMENSION(Nqua, Np), INTENT(OUT)  :: R                     ! new basis coordinates (input)

       IF (magnetic_model .eq. 1) THEN
          CALL EFIT2(X, Y, R)                                               ! realistic equilibrium:: interpolation from Efit files
       ELSEIF (magnetic_model .eq. 2) THEN
          CALL EFITsa(X, Y, R)                                              ! Realistic with chi circular
       ELSEIF (magnetic_model .eq. 3) THEN
          CALL Solovev(X, Y, R)                                             ! Solovev equilibrium ....
       ELSEIF (magnetic_model .eq. 4) THEN
          CALL Circular(X, Y, R)                                            ! Circular equilibria ....
       ELSE
          WRITE (*, *) 'Error: The index for the magnetic model does not correspond to any choice'
       END IF

    END SUBROUTINE MagnModel

    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  01/02/2024   |     D. I. Palade       |   Initial version
    !  ========================================================================================================================================================

    SUBROUTINE OnlyQ(X7, Y7, Q7, Q8)
       USE constants
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), INTENT(IN)   :: X7, Y7                  ! new basis coordinates (input)
       REAL(KIND=dp), INTENT(OUT)  :: Q7, Q8                  ! new basis coordinates (input)

       ! Local variables
       REAL(KIND=dp), DIMENSION(Np)                       :: X8, Y8    ! Effective (star, "s") electric/magnetic fields
       REAL(KIND=dp), DIMENSION(Nqua, Np)                :: R8       ! Effective (star, "s") electric/magnetic fields

       X8 = X7
       Y8 = Y7

       IF (magnetic_model .eq. 1) THEN
          CALL EFIT2(X8, Y8, R8)                                               ! realistic equilibrium:: interpolation from Efit files
       ELSEIF (magnetic_model .eq. 2) THEN
          CALL EFITsa(X8, Y8, R8)                                              ! Realistic with chi circulari
       ELSEIF (magnetic_model .eq. 3) THEN
          CALL Solovev(X8, Y8, R8)                                             ! Solovev equilibrium ....
       ELSEIF (magnetic_model .eq. 4) THEN
          CALL Circular(X8, Y8, R8)                                            ! Circular equilibria ....
       ELSE
          WRITE (*, *) 'Error: The index for the magnetic model does not correspond to any choice'
       END IF

       Q7 = R8(9, 1)
       Q8 = R8(1, 1)

    END SUBROUTINE OnlyQ

    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  01/02/2024   |     D. I. Palade       |   Initial version
    !  ========================================================================================================================================================

    SUBROUTINE Onlycoord(X7, Y7, Z7, q17, q27, q37, B7, Vp)
       USE constants
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN)   :: X7, Y7, Z7, Vp                 ! new basis coordinates (input)
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT)  :: q17, q27, q37, B7                ! new basis coordinates (input)

       ! Local variables
       REAL(KIND=dp), DIMENSION(Np)                       :: X8, Y8, Q7, hx, hy, hz    ! Effective (star, "s") electric/magnetic fields
       REAL(KIND=dp), DIMENSION(Nqua,  Np)                :: R8          ! Effective (star, "s") electric/magnetic fields

       X8 = X7
       Y8 = Y7
       hx = 1.0                                                      ! lame coefficient for x = R -> hx = 1
       hy = 1.0                                                      ! lame coefficient for y = Z -> hy = 1
       hz = X8                                                        ! lame coefficient for z = zeta -> hz = R = X

       IF (magnetic_model .eq. 1) THEN
          CALL EFIT2(X8, Y8, R8)                                               ! realistic equilibrium:: interpolation from Efit files
       ELSEIF (magnetic_model .eq. 2) THEN
          CALL EFITsa(X8, Y8, R8)                                              ! Realistic with chi circulari
       ELSEIF (magnetic_model .eq. 3) THEN
          CALL Solovev(X8, Y8, R8)                                             ! Solovev equilibrium ....
       ELSEIF (magnetic_model .eq. 4) THEN
          CALL Circular(X8, Y8, R8)                                            ! Circular equilibria ....
       ELSE
          WRITE (*, *) 'Error: The index for the magnetic model does not correspond to any choice'
       END IF

       B7 = sqrt((R8(3, :)/X8)**2*hx**2 + (R8(2, :)/X8)**2*hy**2 + (R8(7, :)/X8**2)**2*hz**2)             ! the norm of the magnetic field (with the covariant components)

       q17 = C1*R8(11, :)
       q27 = C2*(Z7 - R8(9, :)*R8(14, :))
       q37 = C3*R8(14, :)

    END SUBROUTINE Onlycoord


    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  15/01/2026   |     D. I. Palade       |   Initial version
    !  ========================================================================================================================================================

    SUBROUTINE Onlycoord_new(X7, Y7, Z7, q17)!, q27, q37, B7, Vp)
       USE constants
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), INTENT(IN)   :: X7, Y7, Z7!, Vp                 ! new basis coordinates (input)
       REAL(KIND=dp), INTENT(OUT)  :: q17!, q27, q37, B7                ! new basis coordinates (input)

       ! Local variables
       REAL(KIND=dp)         :: xi,yi,zi,rhot,rr, psi, psiradius, X1, Y1, Xef,Yef, F1, F2, F3, F4, efit_vals
      INTEGER :: poz1,poz2, jax
   
        xi = X0
        yi = Y0
        zi = Z0
	rr  = sqrt((xi - 1.0)**2 + yi**2)

if (magnetic_model == 4) then

         ! rho_t and derivatives
         rhot   = rr/a0
      else if (magnetic_model == 3) then

         ! psi and derivatives (analytical model)
         psi   = (amp*( (alfa**2*(-1.0 + xi**2)**2)/4.0_dp + (-gama + xi**2)*yi**2 )) / (2.0_dp*(1.0_dp + alfa**2))
         psiradius = (amp*( (alfa**2*(-1.0 + (1.0_dp + a0)**2)**2)/4.0_dp )) / (2.0_dp*(1.0_dp + alfa**2))
         rr        = sqrt(psi/psiradius)*a0
         rhot  = rr/a0
      else
         if ((magnetic_model == 1) .or. (magnetic_model == 2)) then
            ! Grid indices
            poz1 = modulo(int((xi - minR)/stepR), NgridR)
            poz2 = modulo(int((yi - minZ)/stepZ), NgridZ)

            X1 = minR + poz1*stepR
            Y1 = minZ + poz2*stepZ

            Xef = (xi - X1)/stepR
            Yef = (yi - Y1)/stepZ

            do jax = 11, 11
               F1 = Efit_data(NgridR*NgridZ*(jax - 1) +  poz1   *NgridR +  poz2)
               F2 = Efit_data(NgridR*NgridZ*(jax - 1) + (poz1+1)*NgridR +  poz2)
               F3 = Efit_data(NgridR*NgridZ*(jax - 1) +  poz1   *NgridR + (poz2+1))
               F4 = Efit_data(NgridR*NgridZ*(jax - 1) + (poz1+1)*NgridR + (poz2+1))

               efit_vals = F1 + (F2 - F1)*Xef + (F3 - F1)*Yef + Xef*Yef*(F1 + F4 - F2 - F3)
            end do

            rhot  = efit_vals
         end if
      end if
      
       q17 = rhot
      
    END SUBROUTINE Onlycoord_new


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
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN)   :: X, Y, Z                                             ! new basis coordinates (input)
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT)  :: B                                                  ! Bx,  By, Bz, B components in the x,y,z system

       ! Local variables
       REAL(KIND=dp), DIMENSION(Np)                    :: hx, hy, hz, Bx, By, Bz, psir, psiz, Fpsi            ! lame coefficients
       REAL(KIND=dp), DIMENSION(Nqua, Np)               :: R                                                   ! G_ij = grad(qj).dr/dxi;

       hx = 1.0                                                       ! lame coefficient for x = R -> hx = 1
       hy = 1.0                                                       ! lame coefficient for y = Z -> hy = 1
       hz = X                                                         ! lame coefficient for z = zeta -> hz = R = X

       ! ------------------------------------------------------------------------------------------------------------------------
       ! R contains: d_r\psi, d_z\psi, F(psi), rot(b), grad(B)
       ! ------------------------------------------------------------------------------------------------------------------------

       CALL MagnModel(X, Z, R)

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
