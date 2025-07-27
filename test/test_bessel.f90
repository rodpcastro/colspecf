module test_bessel
! Test of Bessel functions.

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use csf_kinds, only: wp
  use csf_constants, only: nan, ninf, pinf
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use csf_bessel, only: besselj0, besselj1, bessely0, bessely1
  use specfun_evaluation, only: eval_write

  implicit none
  private
  public :: collect_bessel_tests
  
contains

  subroutine collect_bessel_tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("test_besselj0", test_besselj0), &
      new_unittest("test_besselj0_extremes", test_besselj0_extremes), &
      new_unittest("test_besselj1", test_besselj1), &
      new_unittest("test_besselj1_extremes", test_besselj1_extremes), &
      new_unittest("test_bessely0", test_bessely0), &
      new_unittest("test_bessely0_extremes", test_bessely0_extremes), &
      new_unittest("test_bessely1", test_bessely1), &
      new_unittest("test_bessely1_extremes", test_bessely1_extremes) &
    ]
  end subroutine collect_bessel_tests

  subroutine test_besselj0(error)
    type(error_type), allocatable, intent(out) :: error
    logical, allocatable :: specfun_ic(:)
    character(len=100) :: file
    
    file = 'bessel_j0.csv'
    call eval_write(besselj0, file, specfun_ic)

    call check(error, all(specfun_ic))
    if (allocated(error)) return
  end subroutine test_besselj0

  subroutine test_besselj0_extremes(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: besselj0w

    besselj0w = besselj0(ninf())
    call check(error, besselj0w, 0.0_wp)
    if (allocated(error)) return

    besselj0w = besselj0(pinf())
    call check(error, besselj0w, 0.0_wp)
    if (allocated(error)) return
  end subroutine test_besselj0_extremes

  subroutine test_besselj1(error)
    type(error_type), allocatable, intent(out) :: error
    logical, allocatable :: specfun_ic(:)
    character(len=100) :: file
    
    file = 'bessel_j1.csv'
    call eval_write(besselj1, file, specfun_ic)

    call check(error, all(specfun_ic))
    if (allocated(error)) return
  end subroutine test_besselj1

  subroutine test_besselj1_extremes(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: besselj1w

    besselj1w = besselj1(ninf())
    call check(error, besselj1w, 0.0_wp)
    if (allocated(error)) return

    besselj1w = besselj1(pinf())
    call check(error, besselj1w, 0.0_wp)
    if (allocated(error)) return
  end subroutine test_besselj1_extremes

  subroutine test_bessely0(error)
    type(error_type), allocatable, intent(out) :: error
    logical, allocatable :: specfun_ic(:)
    character(len=100) :: file
    
    file = 'bessel_y0.csv'
    call eval_write(bessely0, file, specfun_ic)

    call check(error, all(specfun_ic))
    if (allocated(error)) return
  end subroutine test_bessely0

  subroutine test_bessely0_extremes(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: bessely0w

    bessely0w = bessely0(-1.0_wp)
    call check(error, ieee_is_nan(bessely0w))
    if (allocated(error)) return

    bessely0w = bessely0(0.0_wp)
    call check(error, bessely0w, ninf())
    if (allocated(error)) return

    bessely0w = bessely0(pinf())
    call check(error, bessely0w, 0.0_wp)
    if (allocated(error)) return
  end subroutine test_bessely0_extremes

  subroutine test_bessely1(error)
    type(error_type), allocatable, intent(out) :: error
    logical, allocatable :: specfun_ic(:)
    character(len=100) :: file
    
    file = 'bessel_y1.csv'
    call eval_write(bessely1, file, specfun_ic)

    call check(error, all(specfun_ic))
    if (allocated(error)) return
  end subroutine test_bessely1

  subroutine test_bessely1_extremes(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: bessely1w

    bessely1w = bessely1(-1.0_wp)
    call check(error, ieee_is_nan(bessely1w))
    if (allocated(error)) return

    bessely1w = bessely1(0.0_wp)
    call check(error, bessely1w, ninf())
    if (allocated(error)) return

    bessely1w = bessely1(pinf())
    call check(error, bessely1w, 0.0_wp)
    if (allocated(error)) return
  end subroutine test_bessely1_extremes

end module test_bessel
