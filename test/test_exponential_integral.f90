module test_exponential_integral
! Test of exponential integrals Ei and E1.

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use csf_kinds, only: wp
  use csf_constants, only: pi, ninf, pinf
  use csf_exponential_integral, only: ei, e1
  use specfun_evaluation, only: eval_write

  implicit none
  private
  public :: collect_exponential_integral_tests
  
contains

  subroutine collect_exponential_integral_tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("test_ei", test_ei), &
      new_unittest("test_ei_extremes", test_ei_extremes), &
      new_unittest("test_e1x", test_e1x), &
      new_unittest("test_e1x_extremes", test_e1x_extremes), &
      new_unittest("test_e1z", test_e1z), &
      new_unittest("test_e1z_extremes", test_e1z_extremes) &
    ]
  end subroutine collect_exponential_integral_tests

  subroutine test_ei(error)
    type(error_type), allocatable, intent(out) :: error
    logical, allocatable :: specfun_ic(:)
    character(len=100) :: file
    
    file = 'exponential_integral_ei.csv'
    call eval_write(ei, file, specfun_ic)

    call check(error, all(specfun_ic))
    if (allocated(error)) return
  end subroutine test_ei

  subroutine test_ei_extremes(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: eiw

    eiw = ei(0.0_wp)
    call check(error, eiw, ninf())
    if (allocated(error)) return

    eiw = ei(ninf())
    call check(error, eiw, 0.0_wp)
    if (allocated(error)) return

    eiw = ei(pinf())
    call check(error, eiw, pinf())
    if (allocated(error)) return
  end subroutine test_ei_extremes

  subroutine test_e1x(error)
    type(error_type), allocatable, intent(out) :: error
    logical, allocatable :: specfun_ic(:)
    character(len=100) :: file
    
    file = 'exponential_integral_e1x.csv'
    call eval_write(e1x, file, specfun_ic)

    call check(error, all(specfun_ic))
    if (allocated(error)) return

  contains

    real(wp) function e1x(x) 
      real(wp), intent(in) :: x
      e1x = e1(x)
    end function e1x
  end subroutine test_e1x

  subroutine test_e1x_extremes(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: e1xw

    e1xw = e1(0.0_wp)
    call check(error, e1xw, pinf())
    if (allocated(error)) return

    e1xw = e1(pinf())
    call check(error, e1xw, 0.0_wp)
    if (allocated(error)) return
  end subroutine test_e1x_extremes

  subroutine test_e1z(error)
    type(error_type), allocatable, intent(out) :: error
    logical, allocatable :: specfun_ic(:)
    character(len=100) :: file
    
    file = 'exponential_integral_e1z.csv'
    call eval_write(e1z, file, specfun_ic)

    call check(error, all(specfun_ic))
    if (allocated(error)) return

  contains

    complex(wp) function e1z(z) 
      complex(wp), intent(in) :: z
      e1z = e1(z)
    end function e1z
  end subroutine test_e1z

  subroutine test_e1z_extremes(error)
    type(error_type), allocatable, intent(out) :: error
    complex(wp) :: z, e1zw, e1zr

    z = (0.0_wp, 0.0_wp)
    e1zw = e1(z)
    e1zr = cmplx(pinf(), -pi)
    call check(error, e1zw, e1zr)
    if (allocated(error)) return

    z = cmplx(ninf(), 0.0_wp)
    e1zw = e1(z)
    e1zr = cmplx(ninf(), -pi)
    call check(error, e1zw, e1zr)
    if (allocated(error)) return

    z = cmplx(pinf(), 0.0_wp)
    e1zw = e1(z)
    e1zr = (0.0_wp, 0.0_wp)
    call check(error, e1zw, e1zr)
    if (allocated(error)) return

    z = cmplx(0.0_wp, ninf())
    e1zw = e1(z)
    e1zr = (0.0_wp, 0.0_wp)
    call check(error, e1zw, e1zr)
    if (allocated(error)) return

    z = cmplx(0.0_wp, pinf())
    e1zw = e1(z)
    e1zr = (0.0_wp, 0.0_wp)
    call check(error, e1zw, e1zr)
    if (allocated(error)) return

    z = cmplx(pinf(), pinf())
    e1zw = e1(z)
    e1zr = (0.0_wp, 0.0_wp)
    call check(error, e1zw, e1zr)
    if (allocated(error)) return

    z = cmplx(pinf(), ninf())
    e1zw = e1(z)
    e1zr = (0.0_wp, 0.0_wp)
    call check(error, e1zw, e1zr)
    if (allocated(error)) return
  end subroutine test_e1z_extremes

end module test_exponential_integral
