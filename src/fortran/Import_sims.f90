!========================================================================================================
! Module: sims
!
! Purpose:
!   Import parameters for many simulations from the file "Sim_XXX.dat" or "DB_XXX.dat"
!   into array1. The files are placed in .../Sims/ and generated externally
!   (T3ST GUI). This module should be called in main only once, before the run loop.
!
! Expected file layout:
!   optional comment line beginning with '#'
!   number of runs
!   number of parameter columns
!   one header row containing parameter names
!   one numeric row per run
!
! Example:
!   # mode=manual | scenario=CycloneBaseCase | ...
!   2
!   73
!   USE_larmor USE_coll ... Lns Lts
!   0 0 ... 0 0
!   0 0 ... 1 0
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
   ! Public data shared with the main program and constants module.
   !
   ! array1(run, parameter_index) stores the numeric values read from the selected
   ! Sim_XXX.dat or DB_XXX.dat file. The constants module later accesses values by
   ! parameter name through parameter_index(), so column order can change safely as long
   ! as the header row is correct.
   !---------------------------------------------------------------------------------
   INTEGER, PARAMETER :: rp = selected_real_kind(7)

   REAL(KIND=rp), ALLOCATABLE, DIMENSION(:, :) :: array1
   INTEGER                                     :: rows1
   INTEGER                                     :: cols1

   CHARACTER(LEN=3)                    :: noofsimstr
   INTEGER                             :: noofsim
   INTEGER, PARAMETER                  :: name_len = 64
   CHARACTER(LEN=name_len), ALLOCATABLE :: param_names(:)
   CHARACTER(LEN=:), ALLOCATABLE      :: project_root

   ! 'Sim' or 'DB'. This selects Sims/Sim_XXX.dat or Sims/DB_XXX.dat and also
   ! determines the output folder name under data/.
   CHARACTER(LEN=3)   :: file_type

   ! First-line comment from the .dat file, e.g.:
   !   "Manual | scenario: CycloneBaseCase | sweeps: none | exported: 2026-03-23"
   ! Empty string if the file has no comment line (old-format files).
   ! Accessible anywhere via USE sims.
   CHARACTER(LEN=512) :: sim_comment

   PUBLIC :: load_simulation_matrix, array1, rows1, cols1, parameter_index, parameter_exists, file_type, sim_comment, project_root

CONTAINS

   !=====================================================================================================
   ! load_simulation_matrix
   !
   ! Main entry point for this module.
   !
   ! It supports two modes:
   !   1. Interactive mode: main calls load_simulation_matrix(type_choice), and this routine asks
   !      the user which file type and file number to load.
   !   2. Command-line mode: main passes type_choice_arg and sim_number_arg, usually
   !      from T3ST 1 1 or T3ST 2 3.
   !
   ! On return:
   !   - array1 contains all run rows from the selected file.
   !   - rows1 is the number of runs in that file.
   !   - cols1 is the number of parameters per run.
   !   - param_names stores the header names, used for name-based lookup.
   !   - data/Sim_XXX or data/DB_XXX and Run_XXXXX folders exist.
   !=====================================================================================================
   SUBROUTINE load_simulation_matrix(type_choice, type_choice_arg, sim_number_arg)
      CHARACTER(LEN=10000)          :: CWD
      CHARACTER(LEN=:), ALLOCATABLE :: source
      CHARACTER(LEN=:), ALLOCATABLE :: output_dir
      CHARACTER(LEN=:), ALLOCATABLE :: run_dir
      CHARACTER(LEN=35)             :: aux
      LOGICAL                       :: exists
      INTEGER                       :: i
      CHARACTER(LEN=512)            :: first_line
      CHARACTER(LEN=10000)          :: header_line
      INTEGER, INTENT(OUT)          :: type_choice
      INTEGER, INTENT(IN), OPTIONAL :: type_choice_arg
      INTEGER, INTENT(IN), OPTIONAL :: sim_number_arg

      !---------------------------------------------------------------------------------
      ! Get file type from command-line argument or interactive prompt.
      !---------------------------------------------------------------------------------
      IF (PRESENT(type_choice_arg)) THEN
         type_choice = type_choice_arg
      ELSE
         PRINT *, 'Please select the file type to load:'
         PRINT *, '  1 - Sim (Sims/Sim_XXX.dat)'
         PRINT *, '  2 - DB  (Sims/DB_XXX.dat)'
         READ  *, type_choice
      END IF

      DO WHILE (type_choice /= 1 .AND. type_choice /= 2)
         IF (PRESENT(type_choice_arg)) THEN
            ERROR STOP 'Invalid file type argument. Use 1 for Sim or 2 for DB.'
         END IF
         PRINT *, 'Invalid choice. Please enter 1 (Sim) or 2 (DB):'
         READ  *, type_choice
      END DO

      IF (type_choice == 1) THEN
         file_type = 'Sim'
      ELSE
         file_type = 'DB'
      END IF

      !---------------------------------------------------------------------------------
      ! Locate project root, independent of whether the executable is run
      ! from scripts/, src/fortran/, build/, or the repository root.
      !---------------------------------------------------------------------------------
      CALL GetCWD(CWD)
      project_root = find_project_root(TRIM(ADJUSTL(CWD)))

      !---------------------------------------------------------------------------------
      ! Get file number from command-line argument or interactive prompt.
      !---------------------------------------------------------------------------------
      IF (PRESENT(sim_number_arg)) THEN
         noofsim = sim_number_arg
      ELSE
         PRINT *, 'Please enter the number of the '//TRIM(file_type)//' file:'
         READ  *, noofsim
      END IF

      DO WHILE (noofsim /= INT(noofsim) .OR. noofsim < 1 .OR. noofsim > 999)
         IF (PRESENT(sim_number_arg)) THEN
            ERROR STOP 'Invalid file number argument. Use an integer from 1 to 999.'
         END IF
         PRINT *, 'The number does not exist, please try again!'
         READ  *, noofsim
      END DO

      WRITE(noofsimstr, '(I3.3)') noofsim
      source = project_root//'/Sims/'//TRIM(file_type)//'_'//noofsimstr//'.dat'

      !---------------------------------------------------------------------------------
      ! Open the selected parameter file in Sims/.
      !
      ! Example:
      !   type_choice = 1, noofsim = 1  ->  Sims/Sim_001.dat
      !   type_choice = 2, noofsim = 4  ->  Sims/DB_004.dat
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
      ! Prepare output folder: .../data/Sim_XXX  or  .../data/DB_XXX.
      !
      ! The simulation outputs are written into this folder. The original parameter
      ! file is also mirrored here, which makes each data folder self-describing even
      ! if the original Sims/ file is later changed.
      !---------------------------------------------------------------------------------
      output_dir = project_root//'/data/'//TRIM(file_type)//'_'//noofsimstr

      INQUIRE(DIRECTORY=output_dir, EXIST=exists)
      IF (.NOT. exists) THEN
         CALL ensure_directory(output_dir)
      ELSE
         PRINT *, 'The folder '//TRIM(file_type)//'_'//noofsimstr//' already exists'
      END IF

      WRITE(*,*) '------------------------------------------------------------------------------'
      IF (LEN_TRIM(sim_comment) > 0) THEN
         WRITE(*,*) 'File origin: '//TRIM(sim_comment)
         WRITE(*,*) '------------------------------------------------------------------------------'
      END IF

      ! Mirror the parameters into data/Sim_XXX/Sim_XXX.dat or data/DB_XXX/DB_XXX.dat.
      OPEN(UNIT=998, FILE=output_dir//'/'//TRIM(file_type)//'_'//noofsimstr//'.dat', &
           FORM='formatted', ACCESS='stream', STATUS='replace')

      ! Write comment line first (if present), then rows1. The remaining file content
      ! is mirrored below as it is read.
      IF (LEN_TRIM(sim_comment) > 0) THEN
         WRITE(998, '(A)') '# '//TRIM(sim_comment)
      END IF
      WRITE(998, *) rows1

      !---------------------------------------------------------------------------------
      ! Pre-create one output folder for each run row.
      !
      ! Run row 1 -> Run_00001, row 2 -> Run_00002, etc.
      ! T3ST_main later opens binary outputs inside these folders.
      !---------------------------------------------------------------------------------
      DO i = 1, rows1
         WRITE(aux, '(I5.5)') i
         run_dir = output_dir//'/Run_'//TRIM(aux)
         CALL ensure_directory(run_dir)
      END DO

      !---------------------------------------------------------------------------------
      ! Read cols1, header names, and the full parameter matrix.
      !
      ! cols1 is checked against the number of names found in the header. This avoids
      ! silent misalignment between numeric values and parameter names.
      !---------------------------------------------------------------------------------
      READ (999, *) cols1
      WRITE(998, *) cols1

      ! Read and mirror the parameter-name header row
      READ (999, '(A)') header_line
      WRITE(998, '(A)') TRIM(header_line)
      CALL parse_parameter_header(header_line, cols1, param_names)

      ALLOCATE(array1(rows1, cols1))

      ! Read each simulation row. The main program loops over rows1 and calls
      ! load_run_parameters(run), which copies array1(run, :) into named constants.
      DO i = 1, rows1
         READ (999, *) array1(i, :)
         WRITE(998, *) array1(i, :)
      END DO

      CLOSE(999)
      CLOSE(998)

   END SUBROUTINE load_simulation_matrix


   FUNCTION find_project_root(start_dir) RESULT(root)
      !---------------------------------------------------------------------------------
      ! Walk upward from start_dir until the repository root is found.
      !
      ! This lets the executable be launched from the project root, scripts/,
      ! build/, or src/fortran/ without hard-coding absolute paths.
      !---------------------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: start_dir
      CHARACTER(LEN=:), ALLOCATABLE :: root

      CHARACTER(LEN=10000) :: current, parent
      LOGICAL :: has_config, has_sims, has_src
      INTEGER :: slash

      current = TRIM(start_dir)

      DO
         ! The project root is identified by the three directories that define this
         ! codebase layout.
         INQUIRE(DIRECTORY=TRIM(current)//'/config', EXIST=has_config)
         INQUIRE(DIRECTORY=TRIM(current)//'/Sims', EXIST=has_sims)
         INQUIRE(DIRECTORY=TRIM(current)//'/src', EXIST=has_src)

         IF (has_config .AND. has_sims .AND. has_src) THEN
            root = TRIM(current)
            RETURN
         END IF

         ! Move one directory upward. If there is no parent left, fail below.
         slash = INDEX(TRIM(current), '/', BACK=.TRUE.)
         IF (slash <= 1) EXIT

         parent = current(1:slash - 1)
         IF (TRIM(parent) == TRIM(current)) EXIT
         current = TRIM(parent)
      END DO

      ERROR STOP 'Could not locate project root containing config/, Sims/, and src/'
   END FUNCTION find_project_root


   SUBROUTINE ensure_directory(path)
      !---------------------------------------------------------------------------------
      ! Create a directory only if it does not already exist.
      !
      ! mkdir -p is used because it also creates missing parent directories.
      !---------------------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: path
      LOGICAL :: exists

      INQUIRE(DIRECTORY=TRIM(path), EXIST=exists)
      IF (.NOT. exists) THEN
         CALL execute_command_line('mkdir -p "'//TRIM(path)//'"', wait=.true.)
      END IF
   END SUBROUTINE ensure_directory


   SUBROUTINE parse_parameter_header(line, n_expected, names)
      !---------------------------------------------------------------------------------
      ! Split the parameter-name header row into names(:).
      !
      ! The header is space-separated, for example:
      !   USE_larmor USE_coll ... annulus_width pitch_type
      !
      ! The resulting names array is the bridge between human-readable parameter
      ! names and numeric columns in array1.
      !---------------------------------------------------------------------------------
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
         ! Skip any spaces before the next name.
         DO WHILE (p <= L .AND. s(p:p) == ' ')
            p = p + 1
         END DO
         IF (p > L) EXIT

         start = p

         ! Advance until the next space or the end of the line. The substring
         ! start:finish is one parameter name.
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
   END SUBROUTINE parse_parameter_header


   INTEGER FUNCTION parameter_index(name) RESULT(idx)
      !---------------------------------------------------------------------------------
      ! Return the column index associated with a parameter name.
      !
      ! This is used throughout Run_parameters.f90 as:
      !   Ti = pp(parameter_index("Ti"))
      !
      ! If a required parameter is absent, stop immediately. A missing required
      ! parameter would otherwise shift physics silently.
      !---------------------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER :: i

      IF (.NOT. ALLOCATED(param_names)) THEN
         ERROR STOP "param_names not allocated: did you call load_simulation_matrix()?"
      END IF

      idx = 0
      DO i = 1, SIZE(param_names)
         IF (TRIM(param_names(i)) == TRIM(name)) THEN
            idx = i
            RETURN
         END IF
      END DO

      ERROR STOP "Parameter not found in header: "//TRIM(name)
   END FUNCTION parameter_index


   LOGICAL FUNCTION parameter_exists(name) RESULT(found)
      !---------------------------------------------------------------------------------
      ! Check whether an optional parameter exists in the input header.
      !
      ! This is useful when adding a new parameter while still allowing older
      ! Sim_XXX.dat files to run with a default value.
      !---------------------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER :: i

      found = .FALSE.

      IF (.NOT. ALLOCATED(param_names)) THEN
         ERROR STOP "param_names not allocated: did you call load_simulation_matrix()?"
      END IF

      DO i = 1, SIZE(param_names)
         IF (TRIM(param_names(i)) == TRIM(name)) THEN
            found = .TRUE.
            RETURN
         END IF
      END DO
   END FUNCTION parameter_exists


END MODULE sims
