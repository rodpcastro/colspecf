module test_struve
! Test of Struve functions.

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use csf_kinds, only: wp
  use csf_constants, only: ninf, pinf
  use csf_struve, only: struveh0, struveh1
  use specfun_evaluation, only: eval_write

  implicit none
  private
  public :: collect_struve_tests
  
contains

  subroutine collect_struve_tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("test_struveh0", test_struveh0), &
      new_unittest("test_struveh0_extremes", test_struveh0_extremes), &
      new_unittest("test_struveh1", test_struveh1), &
      new_unittest("test_struveh1_extremes", test_struveh1_extremes) &
    ]
  end subroutine collect_struve_tests


  subroutine test_struveh0(error)
    type(error_type), allocatable, intent(out) :: error
    logical, allocatable :: specfun_ic(:)
    character(len=100) :: file
    
    file = 'struve_h0.csv'
    call eval_write(struveh0, file, specfun_ic)

    call check(error, all(specfun_ic))
    if (allocated(error)) return
  end subroutine test_struveh0


  subroutine test_struveh0_extremes(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: struveh0w

    struveh0w = struveh0(ninf())
    call check(error, struveh0w, 0.0_wp)
    if (allocated(error)) return

    struveh0w = struveh0(pinf())
    call check(error, struveh0w, 0.0_wp)
    if (allocated(error)) return
  end subroutine test_struveh0_extremes


  subroutine test_struveh1(error)
    type(error_type), allocatable, intent(out) :: error
    logical, allocatable :: specfun_ic(:)
    character(len=100) :: file
    
    file = 'struve_h1.csv'
    call eval_write(struveh1, file, specfun_ic)

    call check(error, all(specfun_ic))
    if (allocated(error)) return
  end subroutine test_struveh1


  subroutine test_struveh1_extremes(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: struveh1w

    struveh1w = struveh1(ninf())
    call check(error, struveh1w, 0.0_wp)
    if (allocated(error)) return

    struveh1w = struveh1(pinf())
    call check(error, struveh1w, 0.0_wp)
    if (allocated(error)) return
  end subroutine test_struveh1_extremes

end module test_struve
