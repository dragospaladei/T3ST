module rng_mod
  use, intrinsic :: iso_fortran_env, only: int64, real64
  use omp_lib
  implicit none (type, external)

  integer(int64), parameter :: a03 = 2862933555777941757_int64
  integer(int64), parameter :: c03 = 3037000493_int64

  integer(int64) :: seed03
!$omp threadprivate(seed03)

contains

  subroutine rng_init(global_seed)
    integer(int64), intent(in) :: global_seed
!$omp parallel
    seed03 = global_seed + int(omp_get_thread_num(), int64) * 1234567_int64
!$omp end parallel
  end subroutine rng_init

  real(real64) function rng_uniform()
    ! NOTE: relies on 64-bit wraparound overflow (non-standard Fortran).
    seed03 = a03*seed03 + c03
    rng_uniform = 0.5_real64 + real(seed03, real64) * 2.0_real64**(-64)
  end function rng_uniform

end module rng_mod

