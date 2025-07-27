module test_numerror
! Test of numerror functions.

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use csf_kinds, only: wp
  use csf_numerror, only: ismall, isclose 

  implicit none
  private
  public :: collect_numerror_tests
  
contains

  subroutine collect_numerror_tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("test_ismall", test_ismall), &
      new_unittest("test_isclose", test_isclose) &
    ]
  end subroutine collect_numerror_tests


  subroutine test_ismall(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: a, b 

    a = 1.0e-16_wp
    call check(error, ismall(a))
    if (allocated(error)) return

    a = 1.0e-16_wp
    b = 1.0e-1_wp
    call check(error, .not. ismall(a, b))
    if (allocated(error)) return
  end subroutine test_ismall


  subroutine test_isclose(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: a, b 
    complex(wp) :: c, d

    ! For real numbers a and b

    a = 1.0e-8_wp
    b = 2.0e-8_wp
    call check(error, isclose(a, b))
    if (allocated(error)) return

    a = 1.0e-8_wp
    b = 2.1e-8_wp
    call check(error, .not. isclose(a, b))
    if (allocated(error)) return

    call check(error, isclose(a, b, abs_tol=1.0e-7_wp))
    if (allocated(error)) return

    a = 123456789.0_wp
    b = 123456788.0_wp
    call check(error, isclose(a, b))
    if (allocated(error)) return

    a = 123456789.0_wp
    b = 123456787.0_wp
    call check(error, .not. isclose(a, b))
    if (allocated(error)) return

    call check(error, isclose(a, b, rel_tol=1.0e-7_wp))
    if (allocated(error)) return

    ! For complex numbers c and d

    c = (1.0e-8_wp, 2.0e-8_wp)
    d = (2.0e-8_wp, 3.0e-8_wp)
    call check(error, isclose(c, d))
    if (allocated(error)) return

    c = (1.0e-8_wp, 2.0e-8_wp)
    d = (2.1e-8_wp, 3.1e-8_wp)
    call check(error, .not. isclose(c, d))
    if (allocated(error)) return

    call check(error, isclose(c, d, abs_tol=1.0e-7_wp))
    if (allocated(error)) return

    c = (123456789.0_wp, 987654321.0_wp)
    d = (123456788.0_wp, 987654322.0_wp)
    call check(error, isclose(c, d))
    if (allocated(error)) return

    c = (123456789.0_wp, 987654321.0_wp)
    d = (123456787.0_wp, 987654323.0_wp)
    call check(error, .not. isclose(c, d))
    if (allocated(error)) return

    call check(error, isclose(c, d, rel_tol=1.0e-7_wp))
    if (allocated(error)) return
  end subroutine test_isclose

end module test_numerror
