
PROGRAM DNS
   IMPLICIT NONE

!========================================================================================================
!============ VARIABLE DECLARATION ======================================================================
!========================================================================================================
   CHARACTER(LEN=10000) :: CWD                             ! current working directory (for exporting data)
   INTEGER :: lenCWD                          ! length of CWD
   CHARACTER(LEN=:), ALLOCATABLE :: filein, fileout                  ! address where data is exported
   CHARACTER(LEN=35) :: time, shot                      ! auxiliary, for folder definition
   INTEGER :: i, j, nlines, iostat, ngrid1, ngrid2! index of the simulation
   CHARACTER(LEN=:), ALLOCATABLE :: folder                          ! name of the folder for exporting data
   INTEGER :: run ! index of the simulation

   REAL, ALLOCATABLE, DIMENSION(:, :) :: Vx     ! positions, parallel velocity, magnetic moment, probability - for RK4
   REAL, ALLOCATABLE, DIMENSION(:) :: z     ! positions, parallel velocity, magnetic moment, probability - for RK4
   CHARACTER(LEN=80) :: zk                                 ! current working directory (for exporting data)
!========================================================================================================
!============ ADDRESS DECLARATION FOR IMPORT DATA ====================================================
!========================================================================================================
   CALL GetCWD(CWD)
   filein = TRIM(ADJUSTL(CWD))
   fileout = filein
   lenCWD = LEN(filein)
   shot = '54719'
   shot = '77651'
   shot = '68408'

   !============ run = time of data ====================================================
   DO run = 1, 1
      WRITE (time, '(I0)') run
!============ file to import and export ====================================================
      filein = filein(1:lenCWD - 12)//'\'//trim(shot)//'\geqdsk_'//trim(shot)//'_'//trim(time)
      IF (run .LT. 10) THEN
         fileout = fileout(1:lenCWD - 12)//'\'//trim(shot)//'\geqdsk_'//trim(shot)//'_00'//trim(time)//'.dat'
      ELSEIF ((10 .LE. run) .AND. (run .LT. 100)) THEN
         fileout = fileout(1:lenCWD - 12)//'\'//trim(shot)//'\geqdsk_'//trim(shot)//'_0'//trim(time)//'.dat'
      ELSEIF (100 .LE. run) THEN
         fileout = fileout(1:lenCWD - 12)//'\'//trim(shot)//'\geqdsk_'//trim(shot)//'_'//trim(time)//'.dat'
      END IF
      iostat = 0
      OPEN (file=trim(filein), status='old', action='read', iostat=iostat, unit=100)
      OPEN (file=trim(fileout), action='write', iostat=iostat, unit=200)
!============ file dimension ====================================================
      WRITE (*, *) filein
      WRITE (*, *) fileout
!============ read first line to extract grid dimension ====================================================
      READ (100, '(A)', iostat=iostat) zk
      READ (zk(54:56), '(I)', iostat=iostat) ngrid1
      READ (zk(58:61), '(I)', iostat=iostat) ngrid2
!============ see how many lines are there; the first is lost from the previous import ====================================================
      nlines = 0
      DO WHILE (iostat == 0)
         READ (100, '(A)', iostat=iostat) zk
         nlines = nlines + 1
      END DO
      nlines = nlines - 35
!============ start from beggining to import ====================================================
      REWIND (unit=100)

      WRITE (*, *) 'here', nlines
!============ remove the first line ====================================================
      READ (100, '(A)') zk

      DO i = 1, nlines - 2
         READ (100, '(A)') zk
         DO j = 0, 4
            WRITE (200, *) zk(16*j + 1:16*(j + 1))
         END DO
      END DO
      CLOSE (100)
      CLOSE (200)
      WRITE (*, *) run
   END DO
END PROGRAM
