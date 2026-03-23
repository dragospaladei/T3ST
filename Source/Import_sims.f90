!========================================================================================================
! Module: sims
!
! Purpose:
!   Import parameters for many simulations from the file "Sim_XX.dat" or "DB_XX.dat"
!   into array1. The files are placed in .../Sims/ and generated externally
!   (T3ST GUI). This module should be called in main only once, before the run loop.
!
! Record of revisions:
!    Date       | Programmer          | Description
!  -----------  | ------------------- | -----------------------------------------
!  --/--/2021   | D. I. Palade        | Initial version
!  --/--/2023   | L. M. Pomarjanschi  | Corrections and development
!  01/01/2024   | D. I. Palade        | Remake
!  01/07/2025   | D. I. Palade        | Extended applicability
!  03/2026      | D. I. Palade        | Added DB file type support; comment-line handling
!========================================================================================================

MODULE sims
   IMPLICIT NONE
   SAVE

   !---------------------------------------------------------------------------------
   ! Public data
   !---------------------------------------------------------------------------------
   INTEGER, PARAMETER :: dp = selected_real_kind(7)      ! working precision across the entire code

   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: array1 ! storage of simulation parameters
   INTEGER                                    :: rows1  ! number of runs/scenarios in the file
   INTEGER                                    :: cols1  ! number of parameters per run

   CHARACTER(LEN=2)  :: noofsimstr
   INTEGER           :: noofsim
   INTEGER, PARAMETER   :: name_len = 64
   CHARACTER(LEN=name_len), ALLOCATABLE :: param_names(:)

   ! 'Sim' or 'DB' — set during sim_values, used by the rest of the program if needed
   CHARACTER(LEN=3)   :: file_type

   ! First-line comment from the .dat file, e.g.:
   !   "Manual | scenario: CycloneBaseCase | sweeps: none | exported: 2026-03-23"
   ! Empty string if the file has no comment line (old-format files).
   ! Accessible anywhere via USE sims.
   CHARACTER(LEN=512) :: sim_comment

   PUBLIC :: sim_values, array1, rows1, cols1, param_index, file_type, sim_comment

CONTAINS

   !=====================================================================================================
   ! sim_values
   !   - asks user whether to load a Sim or DB file
   !   - asks user for file number (1..99)
   !   - reads Sims/Sim_XX.dat  (or DB_XX.dat) into array1
   !   - mirrors the file into data/Sim_XX/  (or data/DB_XX/)
   !   - ensures data/Sim_XX (or data/DB_XX) and its Run_XXXX folders exist
   !=====================================================================================================
   SUBROUTINE sim_values(type_choice)
      CHARACTER(LEN=10000)          :: CWD          ! current working directory (.../source)
      CHARACTER(LEN=:), ALLOCATABLE :: source       ! path accumulator / trimmed address
      CHARACTER(LEN=35)             :: aux          ! run index as string for folder naming
      INTEGER                       :: lenCWD
      LOGICAL                       :: exists
      INTEGER                       :: i
      CHARACTER(LEN=512)            :: first_line   ! first line read from file (may be a comment)
      CHARACTER(LEN=10000)          :: header_line
      INTEGER                       :: type_choice

      !---------------------------------------------------------------------------------
      ! Ask user for file type
      !---------------------------------------------------------------------------------
      PRINT *, 'Please select the file type to load:'
      PRINT *, '  1 - Sim (Sims/Sim_XX.dat)'
      PRINT *, '  2 - DB  (Sims/DB_XX.dat)'
      READ  *, type_choice

      DO WHILE (type_choice /= 1 .AND. type_choice /= 2)
         PRINT *, 'Invalid choice. Please enter 1 (Sim) or 2 (DB):'
         READ  *, type_choice
      END DO

      IF (type_choice == 1) THEN
         file_type = 'Sim'
      ELSE
         file_type = 'DB'
      END IF

      !---------------------------------------------------------------------------------
      ! Build path to Sims/ folder (sibling of /source)
      !---------------------------------------------------------------------------------
      CALL GetCWD(CWD)
      source = TRIM(ADJUSTL(CWD))
      lenCWD = LEN(source)
      source = source(1:lenCWD - 7)//'/Sims'

      !---------------------------------------------------------------------------------
      ! Ask user for file number
      !---------------------------------------------------------------------------------
      PRINT *, 'Please enter the number of the '//TRIM(file_type)//' file:'
      READ  *, noofsim

      DO WHILE (noofsim /= INT(noofsim) .OR. noofsim < 1 .OR. noofsim > 99)
         PRINT *, 'The number does not exist, please try again!'
         READ  *, noofsim
      END DO

      WRITE(noofsimstr, '(I2.2)') noofsim
      source = source//'/'//TRIM(file_type)//'_'//noofsimstr//'.dat'

      !---------------------------------------------------------------------------------
      ! Open parameter file in Sims/
      !---------------------------------------------------------------------------------
      OPEN(UNIT=999, FILE=source, STATUS='old', ACTION='read')

      !---------------------------------------------------------------------------------
      ! Read (and retain) optional comment line
      !
      ! The GUI may prepend a '#'-prefixed comment as the very first line, e.g.:
      !   # Database | db: Full_database | selection: Vary_Ln_only | ...
      ! Old files without a comment line start directly with the integer row count.
      ! We handle both by reading the first line as a character string and checking.
      !---------------------------------------------------------------------------------
      sim_comment = ''
      READ(999, '(A)') first_line
      first_line = ADJUSTL(first_line)

      IF (first_line(1:1) == '#') THEN
         ! Strip the leading '#' and optional space, store for later use
         sim_comment = ADJUSTL(first_line(2:))
         ! Now read rows1 from the next line
         READ(999, *) rows1
      ELSE
         ! No comment — first_line already contains the rows count as a string
         READ(first_line, *) rows1
      END IF

      !---------------------------------------------------------------------------------
      ! Prepare output folder: .../data/Sim_XX  or  .../data/DB_XX
      !---------------------------------------------------------------------------------
      source = TRIM(ADJUSTL(CWD))
      lenCWD = LEN(source)
      source = source(1:lenCWD - 7)//'/data/'//TRIM(file_type)//'_'//noofsimstr

      INQUIRE(DIRECTORY=source, EXIST=exists)
      IF (.NOT. exists) THEN
         CALL execute_command_line('mkdir '//source, wait=.true.)
      ELSE
         PRINT *, 'The folder '//TRIM(file_type)//'_'//noofsimstr//' already exists'
      END IF

      WRITE(*,*) '------------------------------------------------------------------------------'
      IF (LEN_TRIM(sim_comment) > 0) THEN
         WRITE(*,*) 'File origin: '//TRIM(sim_comment)
         WRITE(*,*) '------------------------------------------------------------------------------'
      END IF

      ! Mirror the parameters into data/Sim_XX/Sim_XX.dat  or  data/DB_XX/DB_XX.dat
      OPEN(UNIT=998, FILE=source//'/'//TRIM(file_type)//'_'//noofsimstr//'.dat', &
           FORM='formatted', ACCESS='stream', STATUS='replace')

      ! Write comment line first (if present), then rows1 — mirrors the source format exactly
      IF (LEN_TRIM(sim_comment) > 0) THEN
         WRITE(998, '(A)') '# '//TRIM(sim_comment)
      END IF
      WRITE(998, *) rows1

      !---------------------------------------------------------------------------------
      ! Read rows, pre-create Run folders
      !---------------------------------------------------------------------------------
      DO i = 1, rows1
         WRITE(aux, *) INT(i)
         aux = ADJUSTL(ADJUSTR(aux))

         IF (i < 10) THEN
            INQUIRE(DIRECTORY=source//'/Run_0000'//TRIM(aux), EXIST=exists)
            IF (.NOT. exists) CALL execute_command_line('mkdir '//source//'/Run_0000'//TRIM(aux), wait=.true.)

         ELSEIF (i < 100) THEN
            INQUIRE(DIRECTORY=source//'/Run_000'//TRIM(aux), EXIST=exists)
            IF (.NOT. exists) CALL execute_command_line('mkdir '//source//'/Run_000'//TRIM(aux), wait=.true.)

         ELSEIF (i < 1000) THEN
            INQUIRE(DIRECTORY=source//'/Run_00'//TRIM(aux), EXIST=exists)
            IF (.NOT. exists) CALL execute_command_line('mkdir '//source//'/Run_00'//TRIM(aux), wait=.true.)
         END IF
      END DO

      !---------------------------------------------------------------------------------
      ! Read cols and the full parameter matrix
      !---------------------------------------------------------------------------------
      READ (999, *) cols1
      WRITE(998, *) cols1

      ! Read and mirror the parameter-name header row
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
      CHARACTER(LEN=*), INTENT(IN)               :: line
      INTEGER, INTENT(IN)                        :: n_expected
      CHARACTER(LEN=*), ALLOCATABLE, INTENT(OUT) :: names(:)

      INTEGER :: n, p, L, start, finish
      CHARACTER(LEN=LEN(line)) :: s

      s = TRIM(line)
      L = LEN_TRIM(s)

      ALLOCATE(names(n_expected))
      names = ""

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
            ERROR STOP "Header has more names than cols1"
         END IF

         names(n) = ADJUSTL(s(start:finish))
      END DO

      IF (n /= n_expected) THEN
         ERROR STOP "Header name count does not match cols1"
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

      ERROR STOP "Parameter not found in header: "//TRIM(name)
   END FUNCTION param_index


END MODULE sims
