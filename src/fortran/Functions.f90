    !   This file contains subroutines for various useful functions
    !  1. mybessK0 :: integral calculus of the modified bessel function of the second kind K(n,z))
    !  2. mybessK1 :: integral calculus of K(n,z)/2^(n-1)/Gamma(n)
    !  3. mybessK2 :: integral calculus of z^n*K(n,z)/2^(n-1)/Gamma(n)
    !  4. mybessI0 :: approximation of I(0,z))
    !  5. mybessJ0 :: approximation of J(0,z))

    !  ========================================================

    ! Subroutine :: mybessK0
    ! Purpose    :: computes the modified Bessel function K_znu(z)
    !            :: It uses the integral definition of K_znu(z) and an equidistant grid for Riemman integration
    !            :: K_znu(z) = int_0^\infty exp(-z*cosh(t))*cosh(znu*t) dt
    !            :: other integral formulations could be explored @Ligia
    !            :: improvements might be possible with the correct evaluation of the end points
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  06/10/2022   |     L.Pomarjanschi  |   Testing
    !  01/12/2022   |     L.Pomarjanschi  |   Final

    SUBROUTINE mybessK0(n, z, znu, rez)

       ! I/O variables
       INTEGER                        :: n
       REAL                           :: znu
       REAL, DIMENSION(n), INTENT(IN)  :: z
       REAL, DIMENSION(n), INTENT(OUT) :: rez

       ! Local variables
       PARAMETER(nx=100)
       INTEGER                        :: i
       REAL                           :: zmax
       REAL*16, DIMENSION(nx)         :: xs, csh0, csh
       REAL*16, DIMENSION(n)          :: zk, rezi
       REAL*16, DIMENSION(n, nx)       :: zintz

       zmax = 10.0

       DO i = 1, nx
          xs(i) = (i - 1)*zmax/nx                     ! equidistant grid on the (0,10) interval;; nx points
       END DO

       csh0 = cosh(xs)
       csh = cosh(znu*xs)

       DO j = 1, n
          zintz(j, :) = exp(-z(j)*csh0)
       END DO

       rezi = matmul(zintz, csh)*zmax/nx
       rez = rezi                                    ! double precision goes single
    END SUBROUTINE mybessK0

    !  ========================================================

    ! Subroutine :: mybessK1
    ! Purpose    :: computes a normalized version of the modified Bessel function K(znu,z)/2^(znu-1)/Gamma(znu)
    !            :: It uses the integral definition of K_znu(z) and an equidistant grid for Riemman integration
    !            :: K_znu(z) = int_0^\infty exp(-z*cosh(t))*cosh(znu*t) dt
    !            :: other integral formulations could be explored @Ligia
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  06/10/2022   |     L.Pomarjanschi  |   Testing
    !  01/12/2022   |     L.Pomarjanschi  |   Final

    SUBROUTINE mybessK1(n, z, znu, rez)

       ! I/O variables
       INTEGER                        :: n
       REAL                           :: znu
       REAL, DIMENSION(n), INTENT(IN)  :: z
       REAL, DIMENSION(n), INTENT(OUT) :: rez

       ! Local variables
       INTEGER                        :: i
       PARAMETER(nx=100)
       REAL                           :: zmax
       REAL*16, DIMENSION(nx)         :: xs, csh0, csh
       REAL*16, DIMENSION(n)          :: zk, rezi
       REAL*16, DIMENSION(n, nx)       :: zintz

       zmax = 10.0

       DO i = 1, nx
          xs(i) = (i - 1)*zmax/nx                                  !equidistant grid on the (0,10) interval;; nx points
       END DO

       csh0 = cosh(xs)
       csh = cosh(znu*xs)/gamma(znu)

       DO j = 1, n
          zintz(j, :) = exp(-z(j)*csh0)
       END DO

       rezi = matmul(zintz, csh)*zmax/nx/2.0**(znu - 1)
       rez = rezi                                               ! double precision goes single
    END SUBROUTINE mybessK1

    !  ========================================================

    ! Subroutine :: mybessK2
    ! Purpose    :: computes a normalized version of the modified Bessel function z^nu*K(nu,z)/2^(nu-1)/Gamma(nu)
    !            :: It uses the integral definition of K_znu(z) and an equidistant grid for Riemmann integration
    !            :: K_znu(z) = int_0^\infty exp(-z*cosh(t))*cosh(znu*t) dt
    !            :: other integral formulations could be explored @Ligia
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  06/10/2022   |     L.Pomarjanschi  |   Testing
    !  01/12/2022   |     L.Pomarjanschi  |   Final

    SUBROUTINE mybessK2(n, z, znu, rez)
       ! I/O variables
       INTEGER                        :: n
       REAL                           :: znu
       REAL, DIMENSION(n), INTENT(IN)  :: z
       REAL, DIMENSION(n), INTENT(OUT) :: rez

       ! Local variables
       INTEGER                        :: i
       PARAMETER(nx=10000)
       REAL                           :: zmax
       REAL*16, DIMENSION(nx)         :: xs, csh0, csh
       REAL*16, DIMENSION(n)          :: zk, rezi
       REAL*16, DIMENSION(n, nx)       :: zintz

       zmax = 10.0

       DO i = 1, nx
          xs(i) = (i - 1)*zmax/nx
       END DO

       csh0 = cosh(xs)
       csh = cosh(znu*xs)/gamma(znu)

       DO j = 1, n
          zintz(j, :) = (z(j)**znu)*exp(-z(j)*csh0)
       END DO

       rezi = matmul(zintz, csh)*zmax/nx/2.0**(znu - 1)
       rez = rezi                                                            ! double precision goes single
    END SUBROUTINE mybessK2

    !  ========================================================

    ! Subroutine :: mybessI0
    ! Purpose    :: compute the modified Bessel function I(0,z)
    !            :: It uses a very good approximation ~ 0.1% global error when weightening with exponential
    !            :: The approximation is due to : J. Olivares et al 2018 J. Phys.: Conf. Ser. 1043 012003
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  06/10/2022   |     L.Pomarjanschi  |   Testing
    !  01/12/2022   |     L.Pomarjanschi  |   Final

    SUBROUTINE mybessI0(n, z, rez)
       ! I/O variables
       INTEGER                        :: n
       REAL, DIMENSION(n), INTENT(IN)  :: z
       REAL, DIMENSION(n), INTENT(OUT) :: rez

       rez = (1.0 + 0.25*z**2)**(-1.0/4.0)*cosh(z)*(1.0 + 0.24273*z**2)/(1.0 + 0.43023*z**2)

    END SUBROUTINE mybessI0

    !  ========================================================

    ! Subroutine :: mybessJ0
    ! Purpose    :: compute the modified Bessel function J(0,z)
    !            :: It uses a very good approximation ~ 0.1% global error when weightening with exponential
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  06/10/2022   |     L.Pomarjanschi  |   Testing
    !  01/12/2022   |     L.Pomarjanschi  |   Final

    SUBROUTINE mybessJ0(n, z, rez)
       USE constants
       ! I/O variables
       INTEGER                        :: n
       REAL                           :: znu
       REAL, DIMENSION(n), INTENT(IN)  :: z
       REAL, DIMENSION(n), INTENT(OUT) :: rez
       ! Local variables
       REAL, DIMENSION(n)             :: cc1, cc2, cc3, cc4, zla, q, d1, d2

       zla = 0.8152
       q = 0.570
       d1 = 0.59
       d2 = 0.454

       cc1 = zla*q/Sqrt(pi)
       cc2 = -q/zla/8.0/Sqrt(pi)
       cc3 = q/zla/Sqrt(pi)
       cc4 = zla/8.0*q/(pi)

    rez =  ((d1 + cc1*z**2 + cc2*Sqrt(1.0 + z**2*zla**4))*Cos(z) + ((cc4*z**2 + (d2 + cc3*z**2)*Sqrt(1.0 + z**2*zla**4))*Sin(z))/z)/((1.0 + q*z**2)*(1.0 + z**2*zla**4)**0.25)

    END SUBROUTINE mybessJ0

    !  ========================================================

    ! Subroutine :: mybessJ0
    ! Purpose    :: compute the modified Bessel function J(0,z)
    !            :: It uses a very good approximation ~ 0.1% global error when weightening with exponential
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  06/10/2022   |     L.Pomarjanschi  |   Testing
    !  01/12/2022   |     L.Pomarjanschi  |   Final

    SUBROUTINE mybessJ0simple(n, z, rez)
       USE constants
       ! I/O variables
       INTEGER                        :: n
       REAL                           :: znu
       REAL, DIMENSION(n), INTENT(IN)  :: z
       REAL, DIMENSION(n), INTENT(OUT) :: rez

    rez =  ((0.59 + 0.26215*z**2 - 0.04931*Sqrt(1. + 0.44162*z**2))*Cos(z) + 0.018488*z*Sin(z) + ((0.454 + 0.39448*z**2)*Sqrt(1. + 0.441628*z**2)*Sin(z))/z)/((1. + 0.441628*z**2)**0.25*(1. + 0.57*z**2))

    END SUBROUTINE mybessJ0simple

    SUBROUTINE mybessJ0simplest(n, z, rez)
       USE constants
       ! I/O variables
       INTEGER                        :: n
       REAL                           :: znu
       REAL, DIMENSION(n), INTENT(IN)  :: z
       REAL, DIMENSION(n), INTENT(OUT) :: rez

       rez = Cos(z - Pi/4.0*tanh(z**2))/(1.0 + Tanh(z**2/5.0)*Sqrt(z))
       !     rez =  ((0.59 + 0.26215*z**2 - 0.04931*Sqrt(1. + 0.44162*z**2))*Cos(z) + 0.018488*z*Sin(z) + ((0.454 + 0.39448*z**2)*Sqrt(1. + 0.441628*z**2)*Sin(z))/z)/((1. + 0.441628*z**2)**0.25*(1. + 0.57*z**2))

    END SUBROUTINE mybessJ0simplest
    !  ========================================================

    ! Subroutine :: Aphi
    ! Purpose    :: simple function used in non-Gaussian transformations
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  06/10/2022   |     L.Pomarjanschi  |   Testing
    !  01/12/2022   |     L.Pomarjanschi  |   Final

    Subroutine Aphi(size, index, a, b, phiz, Afi)
       INTEGER, INTENT(IN)   :: size, index
       REAL, DIMENSION(size) :: phiz, Afi
       REAL                  :: a, b

       IF (index .eq. 1) THEN
          Afi = (1.0 + 2.0*a*phiz + 3.0*b*phiz**2)/sqrt(1.0 + 2.0*a**2 + 3.0*b*(2.0 + 5.0*b))
       ELSEIF (index .eq. 2) THEN
          Afi = sin(a*phiz)
       ELSEIF (index .eq. 3) THEN
          Afi = phiz
       END IF

    END subroutine Aphi
