module test_hypergeometric
! Test of hypergeometric functions.

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use csf_kinds, only: wp
  use csf_numerror, only: isclose
  use csf_hypergeometric, only: hyp2f1
  use readwrite, only: read_test_points, write_test_points

  implicit none
  private
  public :: collect_hypergeometric_tests
  
contains

  subroutine collect_hypergeometric_tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("test_hyp2f1", test_hyp2f1) &
    ]
  end subroutine collect_hypergeometric_tests


  subroutine test_hyp2f1(error)
    type(error_type), allocatable, intent(out) :: error
    logical, allocatable :: specfun_ic(:)
    complex(wp), allocatable :: ref_x(:), ref_y(:), specfun_y(:)
    complex(wp) :: a, b, c, z
    integer :: npts, fileunit, i
    character(len=100) :: file, filename
    
    a = (0.5_wp, 0.0_wp)
    c = (1.5_wp, 0.0_wp)

    file = 'hypergeometric_hyp2f1.csv'

    filename = 'test/test_points/' // file
    call read_test_points(filename, ref_x, ref_y, npts)
    allocate(specfun_y(npts))
    allocate(specfun_ic(npts))

    do i = 1, npts
      b = cmplx(ref_x(i)%re, 0.0_wp)
      z = cmplx(ref_x(i)%im, 0.0_wp)
      specfun_y(i) = hyp2f1(a, b, c, z)
      specfun_ic(i) = isclose(ref_y(i), specfun_y(i), 1.0e-3_wp, 1.0e-3_wp)
    end do
 
    filename = 'test/test_specfun/' // file
    call write_test_points(filename, ref_x, ref_y, specfun_y, specfun_ic)

    call check(error, all(specfun_ic))
    if (allocated(error)) return
  end subroutine test_hyp2f1

end module test_hypergeometric
