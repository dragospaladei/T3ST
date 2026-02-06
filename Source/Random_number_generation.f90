    !   This file contains subroutines for random generators for several PDF's
    !  1.  rand_p_sphere(n,d,x) generates <n> uniform points in a <d>-dimensional sphere of unit radius
    !  2.  PDF_G(m,n,g,mean,wavelength,variance) generates <m> arrays of <n> random normal numbers, with <mean> average
    !                                            and <wavelength> or <variance> variance;; SUBROUTINE FOUND IN "CONSTANTS"
    !  2.  PDF_G_D(m,n,g,mean,wavelength,variance) generates <m> arrays of <n> random normal numbers, with <mean> average
    !                                            and <wavelength> or <variance> variance;; no need to move subroutine in constnats, since is not "present" anymore
    !  3.  PDF_ky(n,zlam,k0,ky) generates <n> random numbers computed for
    !            PDF(ky) = ky/k0 (Exp[-(ky - k0)^2zlam^2/2] - Exp[-(ky + k0)^2 zlam^2/2])
    !  4.  PDF_Cauchy(n,tau,w00,xf) generates <n> random numbers computed for
    !             PDF(xf) = 1/(1+\tau^2(xf-w00)^2)
    !  6.  PDF_w_Bes(n,tau,znu,xf)  generates <n> random numbers computed for
    !             PDF(xf) ~ xf^(nu-1/2)*K(nu-1/2,xf*sqrt(2nu)) !
    !  7.  PDF_k_iso(n,zlam,znu,xf,yf)  generates <n> random numbers computed for
    !             PDF(xf) ~ xf^(nu-1/2)*K(nu-1/2,xf*sqrt(2nu)) !
    !  8.  PDF_k_inert(n,zlam,znu,xf,yf)  generates <n> random numbers computed for
    !             PDF(kx,ky) ~ sqrt(kx^2+ky^2)^znu*exp(-(kx^2+ky^2)/4*(2+znu))  !
    !  9.  PDF_Boltz(n,tau,xf) Boltzmann distribution
    !========================================================================================================
    !  10. wavenum - generates wavenumbers kx,ky,kz; phases ph,phx,phy; frequencies w,wc
    !========================================================================================================
    MODULE random_numbers
    USE constants
    IMPLICIT NONE
    SAVE

    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:),target    :: ph,kx,ky,kperp,kz,w,phx,phy,wc,dm, L, phrmp, Ls, q00wrap
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:),target      :: phs,kxs,kys,kperps,kzs,ws,phxs,phys,wcs,dms, phrmps, q00wraps
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:)  :: Wien
    CONTAINS

    ! Subroutine :: rand_p_sphere(n,x,d)
    ! Purpose    :: uniform generation of points in a dim-dimensional sphere
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  01/09/2022   |     L.Pomarjanschi  |   Testing
    !  01/12/2022   |     L.Pomarjanschi  |   Final

    SUBROUTINE rand_p_sphere(n,d,x)
    USE constants, ONLY: pi
    IMPLICIT NONE

    ! I/O variables
    INTEGER,                       INTENT(IN)  :: n, d
    REAL(KIND=dp),DIMENSION(d, n), INTENT(OUT) :: x

    ! Local variables
    INTEGER :: i
    REAL(KIND=dp), DIMENSION(n)    :: r1, r2, r3

    CALL RANDOM_NUMBER(r1)                                  ! generate random vectors uniformly in (0,1) ;; r1 = radius

    r1=r1**(1.0/d)

    IF (d.EQ.1) THEN                                        ! 1D case
        x(1,:) = 2.0*r1 - 1.0
    ELSE IF (d.EQ.2) THEN                                   ! 2D case
        CALL random_number (r2)                             ! generate random vectors uniformly in (0,1) ;; r2 = poloidal angle
        r2 = 2.0*pi*r2
        x(1,:) = r1*sin(r2)
        x(2,:) = r1*cos(r2)
    ELSE IF (d.EQ.3) THEN                                   ! 3D case
        CALL random_number (r2)                             ! generate random vectors uniformly in (0,1) ;; r2 = poloidal angle
        CALL random_number (r3)                             ! generate random vectors uniformly in (0,1) ;; r3 = azimuthal angle
        r2 = 2.0*pi*r2
        r3 = acos(1.0 - 2.0*r3)
        x(1,:) = r1*sin(r2)*sin(r3)
        x(2,:) = r1*cos(r2)*sin(r3)
        x(3,:) = r1*cos(r3)
    ELSE
        WRITE(*,*) 'Error: keep it in maximum 3D! ~ rand_p_sphere'            ! don't go beyond 3D space
        PAUSE
        STOP
    END IF

    END SUBROUTINE rand_p_sphere


    ! Subroutine :: PDF_Boltz(n,tau,xf)
    ! Purpose    :: generate Boltzmann distributed random numbers ~ exp(-x/T)sqrt(x/T)
    ! Note       :: this is related to the distribution of kinetic energies in the Maxwell-Boltzmann distribution and not the exponetnial bolzmann
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  01/09/2022   |     L.Pomarjanschi  |   Testing
    !  01/12/2022   |     L.Pomarjanschi  |   Final

    SUBROUTINE PDF_Boltz2(n,tau,xf)
    USE constants
    IMPLICIT NONE

    ! I/O variables
    INTEGER,                        INTENT(IN)  :: n
    REAL(KIND=dp),                  INTENT(IN)  :: tau
    REAL(KIND=dp), DIMENSION(n),    INTENT(OUT) :: xf

    ! Local variables
    INTEGER                                     :: i, q
    REAL(KIND=dp)                               :: a, b, c
    REAL(KIND=dp), DIMENSION(MAX(n/500,10))     :: x1, x2, y1, y2
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE    :: x3

    IF(n.lt.10) THEN
        PRINT*, '--------------------------------------------------------------------------------------'
        PRINT*, 'ERROR'
        PRINT*, "Please generate at least 10 numbers! ~ PDF_Boltz"
        PAUSE
        STOP
    END IF

    a   = 0.0
    b   = 20.0
    c   = 1.0/SQRT(2.0*E)
    xf  = 0.*xf
    i   = 0

    DO WHILE (i<n)

        CALL RANDOM_NUMBER(x1)
        CALL RANDOM_NUMBER(x2)

        x1 = a + (b-a)*x1
        y1 = c*x2
        y2 = SQRT(x1)*EXP(-x1)
        x3 = PACK(x1,y2-y1>=0)

        q = MIN(n, i + SIZE(x3))

        xf(i+1:q) = x3(1:q-i)

        i = i + SIZE(x3)
        DEALLOCATE(x3)
    END DO

    xf = xf*tau

    END SUBROUTINE PDF_Boltz2
    
    SUBROUTINE PDF_Boltz(n, tau, xf)
  USE constants, ONLY: dp, E
  IMPLICIT NONE
  INTEGER,               INTENT(IN)  :: n
  REAL(KIND=dp),         INTENT(IN)  :: tau
  REAL(KIND=dp),         INTENT(OUT) :: xf(n)

  INTEGER, PARAMETER :: CHUNK = 65536
  INTEGER :: out, m, j
  REAL(KIND=dp) :: a, b, c
  REAL(KIND=dp), ALLOCATABLE :: x1(:), x2(:), y1(:), y2(:)
  LOGICAL,        ALLOCATABLE :: keep(:)

  IF (n < 10) THEN
     PRINT *, 'ERROR: generate at least 10 numbers ~ PDF_Boltz'
     STOP
  END IF

  a = 0.0_dp;  b = 20.0_dp
  c = 1.0_dp/SQRT(2.0_dp*E)

  xf  = 0.0_dp
  out = 0
  ALLOCATE(x1(CHUNK), x2(CHUNK), y1(CHUNK), y2(CHUNK), keep(CHUNK))

  DO WHILE (out < n)
     m = MIN(CHUNK, n - out)
     CALL RANDOM_NUMBER(x1(1:m))
     CALL RANDOM_NUMBER(x2(1:m))

     x1(1:m) = a + (b - a) * x1(1:m)
     y1(1:m) = c * x2(1:m)
     y2(1:m) = SQRT(x1(1:m)) * EXP(-x1(1:m))

     keep(1:m) = (y2(1:m) >= y1(1:m))
     DO j = 1, m
        IF (keep(j)) THEN
           out = out + 1
           xf(out) = x1(j)
           IF (out == n) EXIT
        END IF
     END DO
  END DO

  xf = tau * xf
  DEALLOCATE(x1, x2, y1, y2, keep)
END SUBROUTINE PDF_Boltz


    !==========================================================================================
    ! Subroutine :: PDF_Cauchy
    ! Purpose    :: Inversion Method for the following PDF (the Cauchy distribution)
    !               PDF(xf) =  1/(pi*\tau*(1+(xf-w)^2/\tau^2)) + 1/(pi*\tau*(1+(xf+w)^2/\tau^2))
    !               Which has the inverse CDF:
    !               CDF^(-1)(y) = Tan(pi(y-1/2))
    !
    ! Note       :: The distribution of the generated numbers has considerable errors in comparison to the real PDF
    !
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  01/09/2022   |     L.Pomarjanschi  |   Testing
    !  01/12/2022   |     L.Pomarjanschi  |   Final

    SUBROUTINE PDF_Cauchy(n,w00,tau,xf)
    USE constants, ONLY: pi
    IMPLICIT NONE

    ! I/O variables
    REAL(KIND=dp),               INTENT(IN)  :: tau, w00
    INTEGER,                     INTENT(IN)  :: n
    REAL(KIND=dp), DIMENSION(n), INTENT(OUT) :: xf

    ! Local variables
    REAL(KIND=dp), DIMENSION(n)              :: x1, x2, aux

    CALL RANDOM_NUMBER(x1)                                  ! generate random vectors uniformly in (0,1)
    CALL RANDOM_NUMBER(x2)                                  ! generate random vectors uniformly in (0,1)

    x2 = TAN((x2-0.5d0)*pi)                                 ! Inverse of CDF computation
    aux = SIGN(1.d0,x1-0.5d0)                                ! this is an auxiliary vector with random +1,-1 values ;; it is used further to make a double-symetric PDF
    xf = x2/tau + w00*aux                                 ! Scaling Cauchy no. with \tau and center them at +w00 and -w00

    END SUBROUTINE PDF_Cauchy

    !  ========================================================

    ! Subroutine :: PDF_ky
    ! Purpose    :: Acceptance-Rejection Method for the following PDF
    !               PDF(x) = x/zk0 (Exp[-(x - zk0)^2/2] - Exp[-(x + zk0)^2 zlam^2/2])
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  01/09/2022   |     L.Pomarjanschi  |   Testing
    !  01/12/2022   |     L.Pomarjanschi  |   Final


    SUBROUTINE PDF_ky(n, zlam, zk0, xf)
  USE constants, ONLY: dp, pi, E
  IMPLICIT NONE
  INTEGER,               INTENT(IN)  :: n
  REAL(KIND=dp),         INTENT(IN)  :: zlam, zk0
  REAL(KIND=dp),         INTENT(OUT) :: xf(n)

  INTEGER, PARAMETER :: CHUNK = 65536   ! size of reusable buffers
  REAL(KIND=dp) :: a, b, c
  INTEGER :: out, m, j
  REAL(KIND=dp), ALLOCATABLE :: x1(:), x2(:), y1(:), y2(:)
  LOGICAL,        ALLOCATABLE :: keep(:)

  IF (n < 10) THEN
     PRINT *, 'ERROR: generate at least 10 numbers ~ PDF_ky'
     STOP
  END IF

  a = -5.0_dp/zlam - zk0
  b =  5.0_dp/zlam + zk0
  c = SQRT(2.0_dp/pi)/E*zlam   ! an upper bound; same as your code

  ALLOCATE(x1(CHUNK), x2(CHUNK), y1(CHUNK), y2(CHUNK), keep(CHUNK))
  out = 0

  DO WHILE (out < n)
     m = MIN(CHUNK, n - out)

     CALL RANDOM_NUMBER(x1(1:m))
     CALL RANDOM_NUMBER(x2(1:m))

     x1(1:m) = a + (b - a) * x1(1:m)
     y1(1:m) = c * x2(1:m)
     y2(1:m) = ( x1(1:m)/zk0 * ( EXP(-zlam*zlam*(x1(1:m)-zk0)**2/2.0_dp) - &
                                 EXP(-zlam*zlam*(x1(1:m)+zk0)**2/2.0_dp) ) ) * zlam/SQRT(2.0_dp*pi)/2.0_dp

     keep(1:m) = (y2(1:m) >= y1(1:m))
     DO j = 1, m
        IF (keep(j)) THEN
           out = out + 1
           xf(out) = x1(j)
           IF (out == n) EXIT
        END IF
     END DO
  END DO

  DEALLOCATE(x1, x2, y1, y2, keep)
END SUBROUTINE PDF_ky

    
    SUBROUTINE PDF_ky2(n,zlam,zk0,xf)
    USE constants, ONLY: pi
    IMPLICIT NONE

    ! I/O variables
    INTEGER,                        INTENT(IN)  :: n                 ! number of random numbers
    REAL(KIND=dp),                  INTENT(IN)  :: zlam, zk0         ! parameters of the distribution
    REAL(KIND=dp), DIMENSION(n),    INTENT(OUT) :: xf                ! final array that stores the numbers

    ! Local variables
    INTEGER                            :: i, q
    REAL(KIND=dp)                               :: a, b, c
    REAL(KIND=dp), DIMENSION(MAX(n/500,10))     :: x1, x2, y1, y2    ! fragment the auxiliary arrays into smaller ones;; speeds up computation
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE    :: x3

    IF(n.lt.10) THEN                                        ! take into account the minimum size(x1,2 etc.)
        PRINT*, '--------------------------------------------------------------------------------------'
        PRINT*, 'ERROR'
        Print*, "Please generate at least 10 numbers! ~ PDF_ky"
        PAUSE
        STOP
    END IF

    a = -5.0/zlam - zk0                                     ! the random no. "x" are generated in the (a,b) interval:: PDF(x)->0, for x<a && x>b
    b =  5.0/zlam + zk0
    c = SQRT(2.0/pi)/2.71828*zlam                           ! c = maxval(PDF(x))
    xf = 0.*xf
    i = 0                                                   ! stores the number of "good" numbers found (that satisfy the condition y2>=y1); initially 0

    DO WHILE (i<n)

        CALL RANDOM_NUMBER(x1)
        CALL RANDOM_NUMBER(x2)

        x1 = a + (b-a)*x1
        y1 = c*x2                                             ! y1 ~ uniform distribution of height c, always greater than or equal to PDF(x2)
        y2 = x1/zk0*(EXP(-zlam**2*(x1 - zk0)**2/2.0)-EXP(-zlam**2*(x1 + zk0)**2/2.0))*zlam/SQRT(2.0*pi)/2.0   ! y2 ~ PDF(x1)
        x3 = PACK(x1,y2-y1>=0)                                ! store in x3 the values of x1 that satisfy the Acceptance-Rejection condition y2>=y1

        q = MIN(n, i + SIZE(x3))                              ! use q to decide where the values of x3 are stored in xf; use MIN in order to avoid overwriting xf

        xf(i+1:q) = x3(1:q-i)                                 ! fill the array xf with the "good" values found

        i = i + SIZE(x3)                                      ! with each iteration, increase i with the size of x3 (variable for each iteration)
        DEALLOCATE(x3)
    END DO

    END SUBROUTINE PDF_ky2

    !  ========================================================

    ! Subroutine :: PDF_G
    ! Purpose    :: Box-Muller transform ;;
    !               Generate 3 vectors of normal no. with 0 mean and l1, l2, l3 square-variance
    !               included in 0constants.mod
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  01/09/2022   |     L.Pomarjanschi  |   Testing
    !  01/12/2022   |     L.Pomarjanschi  |   Final

    SUBROUTINE PDF_G2(m,n,g,mean,wavelength,variance)
    USE constants, ONLY: pi
    IMPLICIT NONE

    ! I/O variables
    INTEGER,               INTENT(IN)  :: m                 ! number of Gaussians
    INTEGER,               INTENT(IN)  :: n                 ! size of each Gaussian
    REAL(KIND=dp), DIMENSION(m, n), INTENT(OUT) :: g                 ! vector that stores the values

    ! Local variables
    INTEGER :: i, k
    REAL(KIND=dp), DIMENSION(m+1, n)            :: r
    REAL(KIND=dp), DIMENSION(m)                 :: avg, var
    REAL(KIND=dp), DIMENSION(m), OPTIONAL       :: mean, variance, wavelength  ! format :: [:[,:,...]]
    ! optional argument; if not provided, will use default;
    ! should only use variance <var0> OR wavelength <l0>
    ! <l0> overwrites <var0>

    IF(PRESENT(mean)) THEN                                  ! default average set to 0.
        avg = mean
    ELSE
        avg = 0.
    END IF

    IF(PRESENT(wavelength)) THEN                            ! default variance set to 1.; wavelength overwrites variance
        var = 1./wavelength
    ELSE IF(PRESENT(variance)) THEN
        var = variance
    ELSE
        var = 1.
    END IF


    IF(MOD(m,2)/=0) THEN                                    ! we need an even iterator <k> for the number of arrays generated with Box-Muller
        k = m + 1
    ELSE
        k = m
    END IF


    DO i = 1, k

        CALL RANDOM_NUMBER(r(i,:))

        WHERE(r(i,:)==0.)                                     ! prevent NaN results (which come from trying to compute log(0.) )
            r(i,:)=10.**(-8)
        END WHERE

        IF(MOD(i,2)==0) THEN                                  ! we always need the index to be an even number in order to use the Box-Muller transform
            IF ((m-i)>=0) THEN                                ! if the initial number of arrays <m> is even, then generate the arrays in pairs of two

                g(i-1,:) = var(i-1)*SQRT(-2.0*LOG(r(i-1,:)))*SIN(2.0*pi*r(i,:))+avg(i-1)
                g(i,:) = var(i)*SQRT(-2.0*LOG(r(i-1,:)))*COS(2.0*pi*r(i,:))+avg(i)

            ELSE                                              ! if the initial number of arrays <m> is odd, then only compute one array for the last index

                g(i-1,:) = var(i-1)*SQRT(-2.0*LOG(r(i-1,:)))*SIN(2.0*pi*r(i,:))+avg(i-1)

            END IF
        END IF
    END DO


    END SUBROUTINE PDF_G2
    
    SUBROUTINE PDF_G(m, n, g, mean, wavelength, variance)
  USE constants, ONLY: dp, pi
  IMPLICIT NONE
  INTEGER,                    INTENT(IN)  :: m, n
  REAL(KIND=dp),              INTENT(OUT) :: g(m, n)
  REAL(KIND=dp), OPTIONAL,    INTENT(IN)  :: mean(m), wavelength(m), variance(m)

  INTEGER :: i
  REAL(KIND=dp) :: avg(m), sig(m)
  REAL(KIND=dp), ALLOCATABLE :: u1(:), u2(:), r(:), s(:), c(:)

  ! Defaults: avg=0, sig=1
  avg = 0.0_dp;  sig = 1.0_dp
  IF (PRESENT(mean))       avg = mean
  IF (PRESENT(wavelength)) sig = 1.0_dp / wavelength
  IF (PRESENT(variance))   sig = variance   ! wavelength wins only if provided; same as your logic

  ALLOCATE(u1(n), u2(n), r(n), s(n), c(n))

  DO i = 1, m, 2
     CALL RANDOM_NUMBER(u1);  CALL RANDOM_NUMBER(u2)
     WHERE (u1 < 1.0e-12_dp) u1 = 1.0e-12_dp
     r = SQRT(-2.0_dp * LOG(u1))
     s = SIN(2.0_dp*pi*u2)
     c = COS(2.0_dp*pi*u2)

     g(i, :) = sig(i) * r * s + avg(i)
     IF (i+1 <= m) g(i+1, :) = sig(i+1) * r * c + avg(i+1)
  END DO

  DEALLOCATE(u1, u2, r, s, c)
END SUBROUTINE PDF_G

    !  ========================================================

    ! Subroutine :: PDF_G_D
    ! Purpose    :: Box-Muller transform ;;
    !               Generate 3 vectors of normal no. with 0 mean and l1, l2, l3 square-variance
    !               included in 0constants.mod
    ! Record of revisions:
    !    Date       |    Programmer       |   Description of change
    !  =======      |   =============     |   ===========================
    !  01/07/2021   |     D.Palade        |   Conversion to Fortran 95
    !  01/09/2022   |     L.Pomarjanschi  |   Testing
    !  01/12/2022   |     L.Pomarjanschi  |   Final
    !  27/12/2023   |     D.Palade        |   without optional arguments

    SUBROUTINE PDF_G_D(m,n,g,mean,wavelength,variance)
    USE constants, ONLY: pi
    IMPLICIT NONE

    ! I/O variables
    INTEGER,               INTENT(IN)  :: m                 ! number of Gaussians
    INTEGER,               INTENT(IN)  :: n                 ! size of each Gaussian
    REAL(KIND=dp), DIMENSION(m, n), INTENT(OUT) :: g                 ! vector that stores the values

    ! Local variables
    INTEGER :: i, k
    REAL(KIND=dp), DIMENSION(m+1, n)            :: r
    REAL(KIND=dp), DIMENSION(m)                 :: avg, var
    REAL(KIND=dp), DIMENSION(m)                 :: mean, variance, wavelength  ! format :: [:[,:,...]]
    ! optional argument; if not provided, will use default;
    ! should only use variance <var0> OR wavelength <l0>
    ! <l0> overwrites <var0>

    avg = mean

    var = 1./wavelength
    var = variance

    IF(MOD(m,2)/=0) THEN                                    ! we need an even iterator <k> for the number of arrays generated with Box-Muller
        k = m + 1
    ELSE
        k = m
    END IF


    DO i = 1, k

        CALL RANDOM_NUMBER(r(i,:))

        WHERE(r(i,:)==0.)                                     ! prevent NaN results (which come from trying to compute log(0.) )
            r(i,:)=10.**(-8)
        END WHERE

        IF(MOD(i,2)==0) THEN                                  ! we always need the index to be an even number in order to use the Box-Muller transform
            IF ((m-i)>=0) THEN                                ! if the initial number of arrays <m> is even, then generate the arrays in pairs of two

                g(i-1,:) = var(i-1)*SQRT(-2.0*LOG(r(i-1,:)))*SIN(2.0*pi*r(i,:))+avg(i-1)
                g(i,:) = var(i)*SQRT(-2.0*LOG(r(i-1,:)))*COS(2.0*pi*r(i,:))+avg(i)

            ELSE                                              ! if the initial number of arrays <m> is odd, then only compute one array for the last index

                g(i-1,:) = var(i-1)*SQRT(-2.0*LOG(r(i-1,:)))*SIN(2.0*pi*r(i,:))+avg(i-1)

            END IF
        END IF
    END DO


    END SUBROUTINE PDF_G_D
    !  ========================================================
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SUBROUTINE wavenum
    !!! Generates wavenumbers, phases and frequencies;; all constant for each 'run'
    !!! Should be called in 'main' at the start of every 'run', after parameters(run)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE wavenum_old(normB, Vstar1, Vstar2, Vstar3, gees)
    IMPLICIT NONE
    REAL(KIND=dp), DIMENSION(Np)    :: w1,w2
    REAL(KIND=dp), DIMENSION(5,Np)  :: g
    INTEGER                :: i
    REAL(KIND=dp)                   :: normB, Vstar1, Vstar2, Vstar3
    REAL(KIND=dp),  DIMENSION(3,3)  :: gees            ! matrix of kperp^2
    real(kind=dp) :: w000, kx000

    IF(ALLOCATED(ph))     THEN
        IF(SIZE(ph).ne.Np)  THEN
            DEALLOCATE(ph,kx,ky,kperp,kz,w, phrmp, phx,phy,wc,dm)

            ALLOCATE(ph(Nc,Np))
            ALLOCATE(kx,ky,kperp, kz, w, phrmp, phx, phy, wc,dm, mold=ph)
        END IF
    ELSE
        ALLOCATE(ph(Nc,Np))
        ALLOCATE(kx,ky,kperp, kz, w, phrmp, phx, phy, wc,dm, mold=ph)
    END IF

    ! ------------------------------------------------------------------------------------------------------------------------
    IF(USE_turb == ON) THEN

        CALL random_number(ph)
        CALL random_number(phrmp)
        CALL random_number(phx)
        CALL random_number(phy)
        CALL random_number(dm)            ! this is for further parrallel fluctuations

        ph  = pi*(2.0*ph-1.0)
        phrmp = pi*(2.0*phrmp-1.0)
        phx = pi*(2.0*phx-1.0)
        phy = pi*(2.0*phy-1.0)
        dm = 2.0*(dm - 0.5)               ! dm = [R0q0/lambdaz*k_paralel]
        dm = int(q00/lambdaz*dm)

        ! ------------------------------------------------------------------------------------------------------------------------
        ! Wavenumbers Kx, Ky, Kz, frequency w (depending on 'corr'; used in Fields); phx, phy, w (used in colliz)
        ! ------------------------------------------------------------------------------------------------------------------------
        w000 = 0.0
        kx000 = 0.0    ! central radial wavenumber in the case of cauchy like kx spectrum

        DO i = 1, Nci     !ITG
            CALL PDF_G(4,Np,g(1:4,:),wavelength=[lambdax,lambdaz,tauc,taucc])           ! corr = 1 ~ (PDF_G, PDF_G, PDF_w_Bes)
            CALL PDF_ky(Np,lambday,k0i,g(5,:))

            kz(i,:) = g(2,:)
            wc(i,:) = g(4,:)
            ky(i,:) = g(5,:)
            if(x_corr .eq.1) then          
               kx(i,:) = g(1,:)
            elseif(x_corr.eq.2) then
                CALL PDF_Cauchy(Np,kx000,lambdax,kx(i,:))
            endif
            if(t_corr .eq.1) then          
                w1      = g(3,:)
            elseif(t_corr.eq.2) then
                CALL PDF_Cauchy(Np,w000,lambdax,w1)
            endif

            ky(i,:) = 1.0*int(ky(i,:)*C2)/C2
            kz(i,:) = 0.0*(C2*ky(i,:)*q00- int(C2*ky(i,:)*q00))/C3 ! ATENTIE, am pus un 0.0 aici fiindca am schimbat reprezentarea
            kz(i,:) = kz(i,:) + dm(i,:)/C3

            kperp(i,:) = (rhoi/R0)*SQRT(kx(i,:)*gees(1,1)*kx(i,:) + 2.0*kx(i,:)*gees(1,2)*ky(i,:) + &
                & 2.0*kx(i,:)*gees(1,3)*kz(i,:) + ky(i,:)*gees(2,2)*ky(i,:)  + &
                & 2.0*ky(i,:)*gees(2,3)*kz(i,:) + kz(i,:)*gees(3,3)*kz(i,:))
            w2 = + Li*(Ti/Ti)*((rhoi/R0/normB)**2)*(kx(i,:)*Vstar1 + ky(i,:)*Vstar2 + kz(i,:)*Vstar3)
            w2 = w2/(1.+Aeff*(Te/Ti)*kperp(i,:)**2.)
            w(i,:) = w1 + 1.0*w2
        END DO

        DO i = Nci+1, Nc    ! TEM

            CALL PDF_G(4,Np,g(1:4,:),wavelength=[lambdax,lambdaz,tauc,taucc])           ! corr = 1 ~ (PDF_G, PDF_G, PDF_w_Bes)
            CALL PDF_ky(Np,lambday,k0e,g(5,:))


            kz(i,:) = g(2,:)

            wc(i,:) = g(4,:)
            ky(i,:) = g(5,:)
            if(x_corr .eq.1) then          
               kx(i,:) = g(1,:)
            elseif(x_corr.eq.2) then
                CALL PDF_Cauchy(Np,kx000,lambdax,kx(i,:))
            endif
            if(t_corr .eq.1) then          
                w1      = g(3,:)
            elseif(t_corr.eq.2) then
                CALL PDF_Cauchy(Np,w000,lambdax,w1)
            endif            
            !            CALL PDF_Cauchy(Np,w000,tauc,w1)

            ky(i,:) = 1.0*int(ky(i,:)*C2)/C2
            kz(i,:) = 0.0*(C2*ky(i,:)*q00- int(C2*ky(i,:)*q00))/C3 ! ATENTIE, am pus un 0.0 aici fiindca am schimbat reprezentarea
            kz(i,:) = kz(i,:) + dm(i,:)/C3

            kperp(i,:) = (rhoi/R0)*SQRT(kx(i,:)*gees(1,1)*kx(i,:) + 2.0*kx(i,:)*gees(1,2)*ky(i,:) + &
                & 2.0*kx(i,:)*gees(1,3)*kz(i,:) + ky(i,:)*gees(2,2)*ky(i,:)  + &
                & 2.0*ky(i,:)*gees(2,3)*kz(i,:) + kz(i,:)*gees(3,3)*kz(i,:))
            w2 = - Ln*(Te/Ti)*((rhoi/R0/normB)**2)*(kx(i,:)*Vstar1 + ky(i,:)*Vstar2 + kz(i,:)*Vstar3)
            w2 = w2/(1.+Aeff*(Te/Ti)*kperp(i,:)**2.)
            w(i,:) = w1 + 1.0*w2
        END DO

        w = real(USE_freq,dp)*w

    ELSEIF(USE_turb /= ON) THEN

        ph  = 0.0
        phx = 0.0
        phy = 0.0
        kx  = 0.0
        ky  = 0.0
        kz  = 0.0
        w   = 0.0
        wc  = 0.0
        kperp = 0.0


    ENDIF


    ! first trajectory is without turbulence
    ph(:,1)  = 0.0
    phx(:,1) = 0.0
    phy(:,1) = 0.0
    kx(:,1)  = 0.0
    ky(:,1)  = 0.0
    kz(:,1)  = 0.0
    w(:,1)   = 0.0
    wc(:,1)  = 0.0
    kperp(:,1) = 0.0

    END SUBROUTINE wavenum_old
    
    
    SUBROUTINE wavenum(normB, Vstar1, Vstar2, Vstar3, gees)
  USE constants
  IMPLICIT NONE
  REAL(KIND=dp), INTENT(IN) :: normB, Vstar1, Vstar2, Vstar3
  REAL(KIND=dp), INTENT(IN) :: gees(3,3)

  REAL(KIND=dp) :: invR0  , rhoi_over_R0, invC3, Te_over_Ti, R2, ATeTi
  REAL(KIND=dp) :: sx, sz, sw1, swc
  INTEGER :: i
  REAL(KIND=dp), ALLOCATABLE :: tmp4(:,:)  ! 4xNp scratch for PDF_G
    real(kind=dp) :: w000, kx000
    
IF (ALLOCATED(ph)) THEN
    IF (SIZE(ph,1) /= Nc .OR. SIZE(ph,2) /= Np) THEN
        DEALLOCATE(ph, kx, ky, kperp, kz, w, phrmp, phx, phy, wc, dm)
        ALLOCATE(ph(Nc,Np))
        ALLOCATE(kx, ky, kperp, kz, w, phrmp, phx, phy, wc, dm, MOLD=ph)
    END IF
ELSE
    ALLOCATE(ph(Nc,Np))
    ALLOCATE(kx, ky, kperp, kz, w, phrmp, phx, phy, wc, dm, MOLD=ph)
END IF

  IF (USE_turb == ON) THEN
     ! phases + parallel shift (dm)
     CALL RANDOM_NUMBER(ph);     ph     = pi*(2.0_dp*ph     - 1.0_dp)
     CALL RANDOM_NUMBER(phrmp);  phrmp  = pi*(2.0_dp*phrmp  - 1.0_dp)
     CALL RANDOM_NUMBER(phx);    phx    = pi*(2.0_dp*phx    - 1.0_dp)
     CALL RANDOM_NUMBER(phy);    phy    = pi*(2.0_dp*phy    - 1.0_dp)
     CALL RANDOM_NUMBER(dm);     dm     = 2.0_dp*(dm - 0.5_dp)
     dm = INT(q00/lambdaz * dm)          ! integer shifts

     ! constants
     invR0        = 1.0_dp/R0
     rhoi_over_R0 = rhoi*invR0
     invC3        = 1.0_dp/C3
     Te_over_Ti   = Te/Ti
     R2           = (rhoi_over_R0/normB)**2
     ATeTi        = Aeff*Te_over_Ti

     sx  = 1.0_dp/lambdax
     sz  = 1.0_dp/lambdaz
     sw1 = 1.0_dp/tauc
     swc = 1.0_dp/taucc
        w000 = 0.0
        kx000 = 0.0    ! central radial wavenumber in the case of cauchy like kx spectrum

     ! ----- ITG block (1..Nci)
     ALLOCATE(tmp4(4,Np))
     DO i = 1, Nci
        CALL PDF_G(4, Np, tmp4, wavelength=[lambdax, lambdaz, tauc, taucc])  ! uses faster Box-Muller

        kz(i,:) = tmp4(2,:)
        wc(i,:) = tmp4(4,:)
        
            if(x_corr .eq.1) then          
               kx(i,:) = tmp4(1,:)
            elseif(x_corr.eq.2) then
                CALL PDF_Cauchy(Np,kx000,lambdax,kx(i,:))
            endif
            if(t_corr .eq.1) then          
                w(i,:)     = tmp4(3,:)
            elseif(t_corr.eq.2) then
                CALL PDF_Cauchy(Np,w000,lambdax,w(i,:))
            endif
            
        CALL PDF_ky(Np, lambday, k0i, ky(i,:))

        ky(i,:) = REAL(INT(ky(i,:)*C2),dp)/C2
        kz(i,:) = REAL(INT(kz(i,:)),dp)/C3     ! aici e ceva cmplicat: kz sunt gausiene cu lambdaz in unitati de qR; de asta kz echivaleaza cu numere poloidale m;

        ! k_perp and frequency correction
        kperp(i,:) = (rhoi_over_R0)*SQRT( kx(i,:)*gees(1,1)*kx(i,:) + 2.0_dp*kx(i,:)*gees(1,2)*ky(i,:) + &
                                         2.0_dp*kx(i,:)*gees(1,3)*kz(i,:) + ky(i,:)*gees(2,2)*ky(i,:)  + &
                                         2.0_dp*ky(i,:)*gees(2,3)*kz(i,:) + kz(i,:)*gees(3,3)*kz(i,:) )

        w(i,:) = w(i,:) + ( + Ln*R2*(kx(i,:)*Vstar1 + ky(i,:)*Vstar2 + kz(i,:)*Vstar3) ) / (1.0_dp + ATeTi*kperp(i,:)**2)
     END DO

     ! ----- TEM block (Nci+1..Nc)
     DO i = Nci+1, Nc
        CALL PDF_G(4, Np, tmp4, wavelength=[lambdax, lambdaz, tauc, taucc])
        kz(i,:) = tmp4(2,:)
        wc(i,:) = tmp4(4,:)
            if(x_corr .eq.1) then          
               kx(i,:) = tmp4(1,:)
            elseif(x_corr.eq.2) then
                CALL PDF_Cauchy(Np,kx000,lambdax,kx(i,:))
            endif
            if(t_corr .eq.1) then          
                w(i,:)      = tmp4(3,:)
            elseif(t_corr.eq.2) then
                CALL PDF_Cauchy(Np,w000,lambdax,w(i,:))
            endif
        CALL PDF_ky(Np, lambday, k0e, ky(i,:))

        ky(i,:) = REAL(INT(ky(i,:)*C2),dp)/C2
        kz(i,:) = REAL(INT(kz(i,:)),dp)/C3     ! aici e ceva cmplicat: kz sunt gausiene cu lambdaz in unitati de qR; de asta kz echivaleaza cu numere poloidale m;

        kperp(i,:) = (rhoi_over_R0)*SQRT( kx(i,:)*gees(1,1)*kx(i,:) + 2.0_dp*kx(i,:)*gees(1,2)*ky(i,:) + &
                                         2.0_dp*kx(i,:)*gees(1,3)*kz(i,:) + ky(i,:)*gees(2,2)*ky(i,:)  + &
                                         2.0_dp*ky(i,:)*gees(2,3)*kz(i,:) + kz(i,:)*gees(3,3)*kz(i,:) )

        w(i,:) = w(i,:) + ( - Ln*Te_over_Ti*R2*(kx(i,:)*Vstar1 + ky(i,:)*Vstar2 + kz(i,:)*Vstar3) ) / (1.0_dp + ATeTi*kperp(i,:)**2)
     END DO
     DEALLOCATE(tmp4)

     w = real(USE_freq,dp) * w
     q00wrap = real(int(C2*ky*q00), dp)

  ELSE
     ph  = 0.0_dp; phx = 0.0_dp; phy = 0.0_dp
     kx  = 0.0_dp; ky  = 0.0_dp; kz  = 0.0_dp
     w   = 0.0_dp; wc  = 0.0_dp; kperp = 0.0_dp;  q00wrap = 0.0_dp;
  END IF

  ! first trajectory without turbulence
  ph(:,1) = 0.0_dp; phx(:,1) = 0.0_dp; phy(:,1) = 0.0_dp
  kx(:,1) = 0.0_dp; ky(:,1)  = 0.0_dp; kz(:,1)  = 0.0_dp
  w(:,1)  = 0.0_dp; wc(:,1)  = 0.0_dp; kperp(:,1) = 0.0_dp; q00wrap(:,1) = 0.0_dp;
END SUBROUTINE wavenum


    SUBROUTINE Larmor_old(ind,X,Y,Z,mut)
    IMPLICIT NONE
    REAL(KIND=dp), DIMENSION(Np)    :: X,Y,Z,mut,B
!    REAL(KIND=dp)   :: ind
    INTEGER                :: i, ind


    IF(ALLOCATED(L))     THEN
        DEALLOCATE(L)
    ENDIF
    ALLOCATE(L(Nc,Np))

    ! ------------------------------------------------------------------------------------------------------------------------
    IF(ind == ON) THEN
        !     CALL magnetic_norm(X,Y,Z,B)                              ! the norm of the magnetic field evaluated at initial positions
        !$OMP PARALLEL DO &
        !$OMP SHARED (mut,B) &
        !$OMP PRIVATE (i)
        DO i=1, Np
            !         L(:,i) = BESSEL_J0(kperp(:,i)*SQRT(2.0*Aw*abs(mut(i))/Zw**2/B(i)))
            L(:,i) = BESSEL_J0(kperp(:,i)*SQRT(2.0*Aw*abs(mut(i))/Zw**2))
        END DO
        !$OMP END parallel DO
    ELSE
        L = 1.d0
    END IF

    END SUBROUTINE Larmor_old
    
    
    SUBROUTINE Larmor(ind, X, Y, Z, mut)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ind
  REAL(KIND=dp), DIMENSION(Np), INTENT(IN) :: X, Y, Z, mut
  REAL(KIND=dp), DIMENSION(Np) :: B  ! if you later compute norm(B)
  INTEGER :: i
  REAL(KIND=dp) :: coef

  IF (ALLOCATED(L)) DEALLOCATE(L)
  ALLOCATE(L(Nc,Np))

  IF (ind == ON) THEN
     coef = SQRT(2.0_dp*Aw/Zw**2)
     !$omp parallel do default(none) shared(L, kperp, mut, coef, Nc, Np) private(i)
     DO i = 1, Np
        L(:,i) = BESSEL_J0(kperp(:,i) * coef * SQRT(ABS(mut(i))))
     END DO
     !$omp end parallel do
  ELSE
     L = 1.0_dp
  END IF
END SUBROUTINE Larmor


    ! = = = = == ======================   =========================
    !  ========================================================
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SUBROUTINE wavenum
    !!! Generates wavenumbers, phases and frequencies;; all constant for each 'run'
    !!! Should be called in 'main' at the start of every 'run', after parameters(run)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE wavenum_single(normB, Vstar1, Vstar2, Vstar3, gees)
    IMPLICIT NONE
    REAL(KIND=dp), DIMENSION(Nc)    :: w1s,w2s
    REAL(KIND=dp), DIMENSION(5,Nc)  :: g
    INTEGER                         :: i
    REAL(KIND=dp)                   :: normB, Vstar1, Vstar2, Vstar3
    REAL(KIND=dp),  DIMENSION(3,3)  :: gees            ! matrix of kperp^2
    real(kind=dp) :: w000, kx000

    IF(ALLOCATED(phs))     THEN
        IF(SIZE(phs).ne.Nc)  THEN
            DEALLOCATE(phs,kxs,kys,kperps,kzs,ws, phrmps, phxs,phys,wcs,dms)

            ALLOCATE(phs(Nc))
            ALLOCATE(kxs,kys,kperps,kzs,ws, phrmps, phxs,phys,wcs,dms, mold=phs)
        END IF
    ELSE
        ALLOCATE(phs(Nc))
        ALLOCATE(kxs,kys,kperps,kzs,ws, phrmps, phxs,phys,wcs,dms, mold=phs)
    END IF
        w000 = 0.0
        kx000 = 0.0    ! central radial wavenumber in the case of cauchy like kx spectrum

    ! ------------------------------------------------------------------------------------------------------------------------
    IF(USE_turb == ON) THEN

        CALL random_number(phs)
        CALL random_number(phrmps)
        CALL random_number(phxs)
        CALL random_number(phys)
        CALL random_number(dms)            ! this is for further parrallel fluctuations

        dms = 2.0*(dms - 0.5)               ! dm = [R0q0/lambdaz*k_paralel]
        dms = int(q00/lambdaz*dms)
        phs  = pi*(2.0*phs-1.0)
        phrmps = pi*(2.0*phrmps-1.0)
        phxs = pi*(2.0*phxs-1.0)
        phys = pi*(2.0*phys-1.0)

        ! ------------------------------------------------------------------------------------------------------------------------
        ! Wavenumbers Kx, Ky, Kz, frequency w (depending on 'corr'; used in Fields); phx, phy, w (used in colliz)
        ! ------------------------------------------------------------------------------------------------------------------------

        CALL PDF_G(4,Nc,g(1:4,:),wavelength=[lambdax,lambdaz,tauc,taucc])           ! corr = 1 ~ (PDF_G, PDF_G, PDF_w_Bes)

        kzs = g(2,:)

        !      CALL PDF_Cauchy(Np,w000,tauc,w1s)
        wcs = g(4,:)
        
            if(x_corr .eq.1) then          
               kxs = g(1,:)
            elseif(x_corr.eq.2) then
                CALL PDF_Cauchy(Np,kx000,lambdax,kxs)
            endif
            if(t_corr .eq.1) then          
                w1s      = g(3,:)
            elseif(t_corr.eq.2) then
                CALL PDF_Cauchy(Np,w000,lambdax,w1s)
            endif
        CALL PDF_ky(Nc,lambday,k0i,g(5,:))

        DO i = 1, Nci     !ITG
            kys(i) = g(5,i)
            kys(i) = 1.0*int(kys(i)*C2)/C2
            kz(i,:) = REAL(INT(kz(i,:),dp))/C3     ! aici e ceva cmplicat: kz sunt gausiene cu lambdaz in unitati de qR; de asta kz echivaleaza cu numere poloidale m;
            kperps(i) = (rhoi/R0)*SQRT(kxs(i)*gees(1,1)*kxs(i) + 2.0*kxs(i)*gees(1,2)*kys(i) + &
                & 2.0*kxs(i)*gees(1,3)*kzs(i) + kys(i)*gees(2,2)*kys(i)  + &
                & 2.0*kys(i)*gees(2,3)*kzs(i) + kzs(i)*gees(3,3)*kzs(i))
            w2s(i) = + Ln*(Ti/Ti)*((rhoi/R0/normB)**2)*(kxs(i)*Vstar1 + kys(i)*Vstar2 + kzs(i)*Vstar3)
            w2s(i) = w2s(i)/(1.+Aeff*(Te/Ti)*kperps(i)**2.)
            ws(i) = w1s(i) + 1.0*w2s(i)
        END DO

        CALL PDF_ky(Nc,lambday,k0e,g(5,:))

        DO i = Nci+1, Nc    ! TEM

            kys(i) = g(5,i)
            kys(i) = 1.0*int(kys(i)*C2)/C2
            kz(i,:) = REAL(INT(kz(i,:),dp))/C3     ! aici e ceva cmplicat: kz sunt gausiene cu lambdaz in unitati de qR; de asta kz echivaleaza cu numere poloidale m;
            kperps(i) = (rhoi/R0)*SQRT(kxs(i)*gees(1,1)*kxs(i) + 2.0*kxs(i)*gees(1,2)*kys(i) + &
                & 2.0*kxs(i)*gees(1,3)*kzs(i) + kys(i)*gees(2,2)*kys(i)  + &
                & 2.0*kys(i)*gees(2,3)*kzs(i) + kzs(i)*gees(3,3)*kzs(i))
            w2s(i) = - Ln*(Te/Ti)*((rhoi/R0/normB)**2)*(kxs(i)*Vstar1 + kys(i)*Vstar2 + kzs(i)*Vstar3)
            w2s(i) = w2s(i)/(1.+Aeff*(Te/Ti)*kperps(i)**2.)
            ws(i) = w1s(i) + 1.0*w2s(i)
        END DO

        ws = real(USE_freq,dp)*ws
        q00wraps = real(int(C2*kys*q00), dp)

    ELSEIF(USE_turb /= ON) THEN

        phs  = 0.0_dp
        phxs = 0.0_dp
        phys = 0.0_dp
        kxs  = 0.0_dp
        kys  = 0.0_dp
        kzs  = 0.0_dp
        ws   = 0.0_dp
        wcs  = 0.0_dp
        kperps = 0.0_dp
        q00wraps = 0.0_dp
    ENDIF

        phs(1)  = 0.0_dp; phxs(1) = 0.0_dp; phys(1) = 0.0_dp; kxs(1)  = 0.0_dp
        kys(1)  = 0.0_dp; kzs(1)  = 0.0_dp; ws(1)   = 0.0_dp; wcs(1)  = 0.0_dp
        kperps(1) = 0.0_dp; q00wraps(1) = 0.0_dp


    END SUBROUTINE wavenum_single


    SUBROUTINE Larmor_single(ind,X,Y,Z,mut)
    IMPLICIT NONE
    REAL(KIND=dp), DIMENSION(Np)    :: X,Y,Z,mut,B
    INTEGER                :: i, ind

    IF(ALLOCATED(Ls))     THEN
        DEALLOCATE(Ls)
    ENDIF
    ALLOCATE(Ls(Nc,Np))

    ! ------------------------------------------------------------------------------------------------------------------------
    IF(ind == ON) THEN
        DO i=1, Np
            Ls(:,i) = BESSEL_J0(kperps*SQRT(2.0*Aw*abs(mut(i))/Zw**2))
        END DO
    ELSE
        Ls = 1.d0
    END IF

    END SUBROUTINE Larmor_single


    SUBROUTINE Wienner()
    USE constants
    IMPLICIT NONE

    ! Local variables
    REAL(KIND = dp),  DIMENSION(5,Np)                     :: grf!, grf1, grf2
    !      REAL(KIND = dp),  DIMENSION(Nt+1,5,Np)                :: Wien
    REAL(KIND = dp)                                       :: dt
    INTEGER                                               :: i


    IF(ALLOCATED(Wien))     THEN
        IF(SIZE(Wien).ne.(Nt+1))  THEN
            DEALLOCATE(Wien)

            ALLOCATE(Wien(Nt+1,5,Np))
        END IF
    ELSE
        ALLOCATE(Wien(Nt+1,5,Np))
    END IF

    dt = (tmax - t0)/REAL(Nt)                                               ! time step

    Wien = 0.0
    if(USE_coll == ON) then
        do i=2,Nt+1
            CALL PDF_G(5,Np,grf(1:5,:))           ! corr = 1 ~ (PDF_G, PDF_G, PDF_w_Bes)
            Wien(i,:,:) = Wien(i-1,:,:) + grf*sqrt(dt)
        enddo
    endif
    END SUBROUTINE Wienner



    END MODULE random_numbers

