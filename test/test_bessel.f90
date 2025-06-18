module test_bessel
! Test of Bessel functions.

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use csf_kinds, only: wp
  use csf_constants, only: nan, ninf, pinf
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use csf_bessel, only: j0x, j1x, y0x, y1x
  use specfun_evaluation, only: eval_write

  implicit none
  private
  public :: collect_bessel_tests
  
contains

  subroutine collect_bessel_tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("test_j0x", test_j0x), &
      new_unittest("test_j0x_extremes", test_j0x_extremes), &
      new_unittest("test_j1x", test_j1x), &
      new_unittest("test_j1x_extremes", test_j1x_extremes), &
      new_unittest("test_y0x", test_y0x), &
      new_unittest("test_y0x_extremes", test_y0x_extremes), &
      new_unittest("test_y1x", test_y1x), &
      new_unittest("test_y1x_extremes", test_y1x_extremes) &
    ]
  end subroutine collect_bessel_tests

  subroutine test_j0x(error)
    type(error_type), allocatable, intent(out) :: error
    logical, allocatable :: specfun_ic(:)
    character(len=100) :: file
    
    file = 'bessel_j0x.csv'
    call eval_write(j0x, file, specfun_ic)

    call check(error, all(specfun_ic))
    if (allocated(error)) return
  end subroutine test_j0x

  subroutine test_j0x_extremes(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: j0xw

    j0xw = j0x(ninf())
    call check(error, j0xw, 0.0_wp)
    if (allocated(error)) return

    j0xw = j0x(pinf())
    call check(error, j0xw, 0.0_wp)
    if (allocated(error)) return
  end subroutine test_j0x_extremes

  subroutine test_j1x(error)
    type(error_type), allocatable, intent(out) :: error
    logical, allocatable :: specfun_ic(:)
    character(len=100) :: file
    
    file = 'bessel_j1x.csv'
    call eval_write(j1x, file, specfun_ic)

    call check(error, all(specfun_ic))
    if (allocated(error)) return
  end subroutine test_j1x

  subroutine test_j1x_extremes(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: j1xw

    j1xw = j1x(ninf())
    call check(error, j1xw, 0.0_wp)
    if (allocated(error)) return

    j1xw = j1x(pinf())
    call check(error, j1xw, 0.0_wp)
    if (allocated(error)) return
  end subroutine test_j1x_extremes

  subroutine test_y0x(error)
    type(error_type), allocatable, intent(out) :: error
    logical, allocatable :: specfun_ic(:)
    character(len=100) :: file
    
    file = 'bessel_y0x.csv'
    call eval_write(y0x, file, specfun_ic)

    call check(error, all(specfun_ic))
    if (allocated(error)) return
  end subroutine test_y0x

  subroutine test_y0x_extremes(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: y0xw

    y0xw = y0x(-1.0_wp)
    call check(error, ieee_is_nan(y0xw))
    if (allocated(error)) return

    y0xw = y0x(0.0_wp)
    call check(error, y0xw, ninf())
    if (allocated(error)) return

    y0xw = y0x(pinf())
    call check(error, y0xw, 0.0_wp)
    if (allocated(error)) return
  end subroutine test_y0x_extremes

  subroutine test_y1x(error)
    type(error_type), allocatable, intent(out) :: error
    logical, allocatable :: specfun_ic(:)
    character(len=100) :: file
    
    file = 'bessel_y1x.csv'
    call eval_write(y1x, file, specfun_ic)

    call check(error, all(specfun_ic))
    if (allocated(error)) return
  end subroutine test_y1x

  subroutine test_y1x_extremes(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: y1xw

    y1xw = y1x(-1.0_wp)
    call check(error, ieee_is_nan(y1xw))
    if (allocated(error)) return

    y1xw = y1x(0.0_wp)
    call check(error, y1xw, ninf())
    if (allocated(error)) return

    y1xw = y1x(pinf())
    call check(error, y1xw, 0.0_wp)
    if (allocated(error)) return
  end subroutine test_y1x_extremes

end module test_bessel
