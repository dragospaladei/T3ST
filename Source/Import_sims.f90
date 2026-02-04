!========================================================================================================
! Module: sims
!
! Purpose:
!   Import parameters for many simulations from the file "Sim_XX.dat" into array1.
!   The Sims files are placed in .../Sims/ and generated externally (Mathematica).
!   This module should be called in main only once, before the run loop.
!
! Record of revisions:
!    Date       | Programmer          | Description
!  -----------  | ------------------- | -----------------------------------------
!  --/--/2021   | D. I. Palade        | Initial version
!  --/--/2023   | L. M. Pomarjanschi  | Corrections and development
!  01/01/2024   | D. I. Palade        | Remake
!  01/07/2025   | D. I. Palade        | Extended applicability
!========================================================================================================

MODULE sims
   IMPLICIT NONE
   SAVE

   !---------------------------------------------------------------------------------
   ! Public data
   !---------------------------------------------------------------------------------
   INTEGER, PARAMETER :: dp = selected_real_kind(7)      ! working precision across the entire code

   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: array1 ! storage of simulation parameters
   INTEGER                                    :: rows1  ! number of runs/scenarios in the Sim file
   INTEGER                                    :: cols1  ! number of parameters per run

   CHARACTER(LEN=2) :: noofsimstr
   INTEGER          :: noofsim
   INTEGER, PARAMETER   :: name_len = 64
   CHARACTER(LEN=name_len), ALLOCATABLE :: param_names(:)

  PUBLIC :: sim_values, array1, rows1, cols1, param_index

CONTAINS

   !=====================================================================================================
   ! sim_values
   !   - asks user for Sim number (1..99)
   !   - reads Sims/Sim_XX.dat into array1
   !   - mirrors the file into data/Sim_XX/Sim_XX.dat
   !   - ensures data/Sim_XX and its Run_XXXX folders exist
   !=====================================================================================================
   SUBROUTINE sim_values
      CHARACTER(LEN=10000)          :: CWD          ! current working directory (.../source)
      CHARACTER(LEN=:), ALLOCATABLE :: source       ! path accumulator / trimmed address
      CHARACTER(LEN=35)             :: aux          ! run index as string for folder naming
      INTEGER                       :: lenCWD
      LOGICAL                       :: exists
      INTEGER                       :: i
	   CHARACTER(LEN=10000) :: header_line
      !---------------------------------------------------------------------------------
      ! Build path to Sims folder (sibling of /source)
      !---------------------------------------------------------------------------------
      CALL GetCWD(CWD)
      source = TRIM(ADJUSTL(CWD))
      lenCWD = LEN(source)
      source = source(1:lenCWD - 7)//'/Sims'

      !---------------------------------------------------------------------------------
      ! Ask user for simulation number
      !---------------------------------------------------------------------------------
      PRINT *, 'Please enter the number of the simulation:'
      READ  *, noofsim

      DO WHILE (noofsim /= INT(noofsim) .OR. noofsim < 1 .OR. noofsim > 99)
         PRINT *, 'The number of the simulation does not exist, please try again!'
         READ  *, noofsim
      END DO

      WRITE(noofsimstr, '(I2.2)') noofsim
      source = source//'/Sim_'//noofsimstr//'.dat'   ! .../Sims/Sim_XX.dat

      !---------------------------------------------------------------------------------
      ! Open parameter file in Sims/
      !---------------------------------------------------------------------------------
      OPEN(UNIT=999, FILE=source, STATUS='old', ACTION='read')

      !---------------------------------------------------------------------------------
      ! Prepare output folder: .../data/Sim_XX
      !---------------------------------------------------------------------------------
      source = TRIM(ADJUSTL(CWD))
      lenCWD = LEN(source)
      source = source(1:lenCWD - 7)//'/data/Sim_'//noofsimstr

      INQUIRE(DIRECTORY=source, EXIST=exists)
      IF (.NOT. exists) THEN
         CALL execute_command_line('mkdir '//source, wait=.true.)
      ELSE
         PRINT *, 'The Simulation folder Sim_'//noofsimstr//' already exists'
      END IF

      WRITE(*,*) '------------------------------------------------------------------------------'

      ! Mirror the parameters into data/Sim_XX/Sim_XX.dat
      OPEN(UNIT=998, FILE=source//'/Sim_'//noofsimstr//'.dat', FORM='formatted', ACCESS='stream', STATUS='replace')

      !---------------------------------------------------------------------------------
      ! Read rows, pre-create Run folders
      !---------------------------------------------------------------------------------
      READ (999, *) rows1
      WRITE(998, *) rows1

      DO i = 1, rows1
         WRITE(aux, *) INT(i)
         aux = ADJUSTL(ADJUSTR(aux))

         IF (i < 10) THEN
            INQUIRE(DIRECTORY=source//'/Run_0000'//aux, EXIST=exists)
            IF (.NOT. exists) CALL execute_command_line('mkdir '//source//'/Run_0000'//aux, wait=.true.)

         ELSEIF (i < 100) THEN
            INQUIRE(DIRECTORY=source//'/Run_000'//aux, EXIST=exists)
            IF (.NOT. exists) CALL execute_command_line('mkdir '//source//'/Run_000'//aux, wait=.true.)

         ELSEIF (i < 1000) THEN
            INQUIRE(DIRECTORY=source//'/Run_00'//aux, EXIST=exists)
            IF (.NOT. exists) CALL execute_command_line('mkdir '//source//'/Run_00'//aux, wait=.true.)
         END IF
      END DO

      !---------------------------------------------------------------------------------
      ! Read cols and the full parameter matrix
      !---------------------------------------------------------------------------------
      READ (999, *) cols1
      WRITE(998, *) cols1

      ! NEW: read and mirror the parameter-name header row
      READ (999, '(A)') header_line
      WRITE(998, '(A)') TRIM(header_line)
      CALL parse_header_names(header_line, cols1, param_names)

      ALLOCATE(array1(rows1, cols1))

      DO i = 1, rows1
         READ (999, *) array1(i, :)
         WRITE(998, *) array1(i, :)
      END DO

      CLOSE(999)
      CLOSE(998)

   END SUBROUTINE sim_values


   SUBROUTINE parse_header_names(line, n_expected, names)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: line
      INTEGER, INTENT(IN)          :: n_expected
      CHARACTER(LEN=*), ALLOCATABLE, INTENT(OUT) :: names(:)

      INTEGER :: i, n, p, L, start, finish
      CHARACTER(LEN=LEN(line)) :: s

      s = TRIM(line)
      L = LEN_TRIM(s)

      ALLOCATE(names(n_expected))
      names = ""   ! clear

      n = 0
      p = 1
      DO WHILE (p <= L)
         ! skip spaces
         DO WHILE (p <= L .AND. s(p:p) == ' ')
            p = p + 1
         END DO
         IF (p > L) EXIT

         start = p
         DO WHILE (p <= L .AND. s(p:p) /= ' ')
            p = p + 1
         END DO
         finish = p - 1

         n = n + 1
         IF (n > n_expected) THEN
            ERROR STOP "Sim header has more names than cols1"
         END IF

         names(n) = ADJUSTL(s(start:finish))
      END DO

      IF (n /= n_expected) THEN
         ERROR STOP "Sim header name count does not match cols1"
      END IF
   END SUBROUTINE parse_header_names


   INTEGER FUNCTION param_index(name) RESULT(idx)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER :: i

      IF (.NOT. ALLOCATED(param_names)) THEN
         ERROR STOP "param_names not allocated: did you call sim_values()?"
      END IF

      idx = 0
      DO i = 1, SIZE(param_names)
         IF (TRIM(param_names(i)) == TRIM(name)) THEN
            idx = i
            RETURN
         END IF
      END DO

      ! Fail loud: missing parameter
      ERROR STOP "Parameter not found in Sim header: "//TRIM(name)
   END FUNCTION param_index


END MODULE sims

