!========================================================================================================
! File: Thread_rng.f90
! Thread-private random-number generator used where OpenMP workers need independent streams.
! Module name kept as rng_mod for compatibility.
!========================================================================================================

MODULE rng_mod
   USE, INTRINSIC :: iso_fortran_env, ONLY: INT64, REAL64
   USE omp_lib
   IMPLICIT NONE(TYPE, EXTERNAL)

   INTEGER(INT64), PARAMETER :: a03 = 2862933555777941757_INT64
   INTEGER(INT64), PARAMETER :: c03 = 3037000493_INT64

   INTEGER(INT64) :: seed03
!$omp threadprivate(seed03)

CONTAINS

   SUBROUTINE rng_init(global_seed)
      INTEGER(INT64), INTENT(in) :: global_seed
!$omp parallel
      seed03 = global_seed + int(omp_get_thread_num(), INT64)*1234567_INT64
!$omp end parallel
   END SUBROUTINE rng_init

   REAL(REAL64) FUNCTION rng_uniform()
      ! NOTE: relies on 64-bit wraparound overflow (non-standard Fortran).
      seed03 = a03*seed03 + c03
      rng_uniform = 0.5_REAL64 + REAL(seed03, REAL64)*2.0_REAL64**(-64)
   END FUNCTION rng_uniform

END MODULE rng_mod
