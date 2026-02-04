    !      ! Subroutine :: initial_conditions
    !      ! Purpose    :: fixes the initial GY coordinates of particles. Consequently their energetic distribution
    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  --/--/2021   |     D. I. Palade       |   Initial version
    !      !  --/--/2023   |     L. M. Pomarjanschi |   Corections and development
    !      !  01/01/2024   |     D. I. Palade       |   Remake
    !  ========================================================================================================================================================

    SUBROUTINE initial_conditions(X, Y, Z, mu, Vp, pb)
       USE constants
       USE random_numbers
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT) :: X, Y, Z, mu, vp, pb

       ! Local variables
       REAL(KIND=dp), DIMENSION(Np)                 :: En, PA, aux
       REAL(KIND=dp), DIMENSION(Np)                 :: phi0, phi1, phi2, B, Bx, By, Bz, q1, q2, q3
       REAL(KIND=dp)                                :: m1, m2
       ! ------------------------------------------------------------------------------------------------------------------------
       ! Pitch angle PA = (-pi, pi) or fixed
       ! ------------------------------------------------------------------------------------------------------------------------

       IF (pitch_type == 1) THEN
          PA = pitch
       ELSEIF (pitch_type == 2) THEN
          CALL random_number(PA)
          PA = (2.*PA - 1.)
!          m1 = sum(PA)/Np
!          m2 = sum(PA**2)/Np       ! this was an old attempt to create a distribution with the correct first 2 moments; i found out that it induces correlation and leads to |PA|>1
!          PA = (PA - m1)*sqrt(1.0/3.0/(m2 - m1**2))
       ELSE
          WRITE (*, *) 'Error [initial_conditions]: unknown type for pitch angle distribution specification'
       END IF

       ! ------------------------------------------------------------------------------------------------------------------------
       ! Energy En distribution :: note, these are kinetic energies only.
       ! ------------------------------------------------------------------------------------------------------------------------

       IF (energy_type == 1) THEN
          En = Ew
       ELSEIF (energy_type == 2) THEN
          CALL PDF_Boltz(Np, Tw, En)                                ! the distribution of real kinetic energies mv^2/2
       ELSE
          WRITE (*, *) 'Error [initial_conditions]: unknown type for energy distribution specification'
       END IF

       Vp = SQRT(2.0*En/Aw)*PA                      !SQRT(2*(En-Phi*Zw*phi0)/Aw)*COS(PA)

       ! ------------------------------------------------------------------------------------------------------------------------
       ! Initial positions
       ! ------------------------------------------------------------------------------------------------------------------------
       CALL initial_position(X, Y, Z, Vp)
       ! ------------------------------------------------------------------------------------------------------------------------
       ! |B| and phi (used for mu and vp)
       ! ------------------------------------------------------------------------------------------------------------------------

       CALL magnetic_norm(X, Y, Z, B)                              ! the norm of the magnetic field evaluated at initial positions
       CALL Onlycoord(X, Y, Z, q1, q2, q3, B, Vp)
       CALL potent_3D(q1, q2, q3, t0, phi1)                            ! the turbulent field values evaluated at initial positions
       !       CALL phi_zero(X,Y,Z,t0,phi2)                             ! the macroscopic electric potential field values evaluated at initial positions

       phi2 = Omgt0**2*Te/(Ti + Te)*(X - 1.0)*(3.0 + X)/2.0

       ! there is a question weather it is ok to evaluate the potential without gyro-average...
       ! ------------------------------------------------------------------------------------------------------------------------
       ! Note:: the energies are shifted with the potential value directly in "mu, vp" for simplicity

       IF (USE_turb == ON) THEN
          pb = E**(-Zw*(Phi*phi1 + phi2)/Tw - 0.5*(Zw*Phi/Tw)**2)
       ELSE
          pb = E**(-Zw*phi2/Tw)
       END IF

       pb = pb*Np/sum(pb)

       mu = En/B*(1.0 - PA**2.)                               !(En-Phi*Zw*phi0)/normB*SIN(PA)**2.

    END SUBROUTINE


    !      ! Subroutine :: initial_positions
    !      ! Purpose    :: fixes the initial GY coordinates of particles (only spatial)
    !  ========================================================================================================================================================
    !      ! Record of revisions:
    !      !    Date       |    Programmer          |   Description of change
    !      !  =======      |   =============        |   =========================
    !      !  --/--/2021   |     D. I. Palade       |   Initial version
    !      !  --/--/2023   |     L. M. Pomarjanschi |   Corections and development
    !      !  01/01/2024   |     D. I. Palade       |   Remake
    !  ========================================================================================================================================================

    SUBROUTINE initial_position(X, Y, Z, Vp)
       USE constants
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT) :: X, Y, Z
       REAL(KIND=dp), DIMENSION(Np)                 :: Vp
       REAL(KIND=dp), DIMENSION(Np)                 :: aux, aux2, aux1

       IF (position_type .eq. 1) THEN
          X = X0
          Y = Y0
          Z = Z0
          !         X = 1.
          !         Y = r00
       ELSEIF (position_type .eq. 2) THEN
          CALL random_number(aux)
          aux = pi*(2.0*aux - 1.0)
          if ((magnetic_model .eq. 1) .or. (magnetic_model .eq. 2)) then
             call psisurf_efit3(X, Y, Vp)
             Z = aux
          elseif (magnetic_model .eq. 3) then
             call psisurf_solov2(X, Y)
             Z = aux
          elseif (magnetic_model .eq. 4) then
             call psisurf_circ(X, Y, Vp)
             Z = aux
          end if

       ELSEIF (position_type .eq. 3) THEN
          CALL random_number(aux)
          CALL random_number(aux2)
          aux = pi*(2.0*aux - 1.0)
          aux2 = a0*sqrt(aux2*(0.55**2 - 0.45**2) + 0.45**2)
          Y = aux2*sin(aux)
          X = 1.0 + aux2*cos(aux)
          Z = 0.0
       ELSE
          WRITE (*, *) 'Error [initial_conditions]: unknown tyoe for GC position distribution specification'
       END IF
       

    END SUBROUTINE
    
    
! a new version (2026) made because routine restructuring
SUBROUTINE potent_3D(q_1, q_2, q_3, time, phi0)
       USE constants
       USE random_numbers
       IMPLICIT NONE

       ! I/O variables
       REAL(KIND=dp), DIMENSION(Np), INTENT(IN) :: q_1, q_2, q_3
       REAL(KIND=dp), INTENT(IN) :: time

       REAL(KIND=dp), DIMENSION(Np), INTENT(OUT):: phi0

       REAL(KIND=dp), DIMENSION(Np) :: gpar, gprim, env
       ! Local variables
       INTEGER                         :: i


      ! Defaults (kept as in your original snippet)
      gpar    = exp((cos(q_3/C3)-1.0_dp)/lbalonz**2)*balloon   ! this must be periodic
      gprim   = -sin(q_3/C3)/lbalonz**2/C3*balloon ! note that this is g'/g
      env     = noballoon + gpar             ! ballooning envelope multiplier

       if(USE_real .eq. OFF) then 
	       !$OMP PARALLEL DO &
	       !$OMP SHARED (q_1,q_2,q_3,time) &
	       !$OMP PRIVATE (i)
	       DO i = 1, Np
		  phi0(i) = sum(L(:,i)*sin(kx(:, i)*q_1(i) + ky(:, i)*q_2(i) + kz(:, i)*q_3(i) - w(:, i)*time + ph(:, i))) ! sin(k_i x+ph_i) terms for all realizations (i) and all wavenumbers; summing over all waves in a realization
	       END DO
	       !$OMP END parallel DO
       elseif(USE_real .eq. ON) then 

	       !$OMP PARALLEL DO &
	       !$OMP SHARED (q_1,q_2,q_3,time) &
	       !$OMP PRIVATE (i)
	       DO i = 1, Np
		  phi0(i) = sum(Ls(:,i)*sin(kxs(:)*q_1(i) + kys(:)*q_2(i) + kzs(:)*q_3(i) - ws(:)*time + phs(:))) 
	       END DO
	       !$OMP END parallel DO
	endif
       phi0 = phi0/sqrt(Nc/2.0)                            ! normalization

    END SUBROUTINE potent_3D

