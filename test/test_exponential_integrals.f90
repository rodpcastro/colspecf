module test_exponential_integrals
! Simple test of special functions.
!
! Author
! ------
! Rodrigo Castro (GitHub: rodpcastro)
!
! History
! -------
! 22-05-2025 - Rodrigo Castro - Original code

  use, intrinsic :: iso_fortran_env, only: real64
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use exponential_integrals, only: ei, e1
  use numerror, only: isclose
  use readwrite, only: read_test_points, write_test_points

  implicit none
  private
  public :: collect_ei_tests
  
contains

  subroutine collect_ei_tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("test_ei", test_ei), &
      new_unittest("test_e1x", test_e1x), &
      new_unittest("test_e1z", test_e1z) &
    ]
  end subroutine collect_ei_tests

  subroutine test_ei(error)
    type(error_type), allocatable, intent(out) :: error
    real(real64), allocatable :: ref_x(:), ref_y(:), specfun_y(:)
    logical, allocatable :: specfun_ic(:)
    character(len=100) :: file, filename
    integer :: npts, fileunit, i
    
    file = 'exponential_integrals_ei.csv'

    filename = 'test/test_points/' // file
    call read_test_points(filename, ref_x, ref_y, npts)
    allocate(specfun_y(npts))
    allocate(specfun_ic(npts))

    do i = 1, npts
      specfun_y(i) = ei(ref_x(i))
      specfun_ic(i) = isclose(ref_y(i), specfun_y(i))
    end do
 
    filename = 'test/test_specfun/' // file
    call write_test_points(filename, ref_x, ref_y, specfun_y, specfun_ic)

    call check(error, all(specfun_ic))
    if (allocated(error)) return
  end subroutine test_ei

  subroutine test_e1x(error)
    type(error_type), allocatable, intent(out) :: error
    real(real64), allocatable :: ref_x(:), ref_y(:), specfun_y(:)
    logical, allocatable :: specfun_ic(:)
    character(len=100) :: file, filename
    integer :: npts, fileunit, i
    
    file = 'exponential_integrals_e1x.csv'

    filename = 'test/test_points/' // file
    call read_test_points(filename, ref_x, ref_y, npts)
    allocate(specfun_y(npts))
    allocate(specfun_ic(npts))

    do i = 1, npts
      specfun_y(i) = e1(ref_x(i))
      specfun_ic(i) = isclose(ref_y(i), specfun_y(i))
    end do
 
    filename = 'test/test_specfun/' // file
    call write_test_points(filename, ref_x, ref_y, specfun_y, specfun_ic)

    call check(error, all(specfun_ic))
    if (allocated(error)) return
  end subroutine test_e1x

  subroutine test_e1z(error)
    type(error_type), allocatable, intent(out) :: error
    complex(real64), allocatable :: ref_x(:), ref_y(:), specfun_y(:)
    logical, allocatable :: specfun_ic(:)
    character(len=100) :: file, filename
    integer :: npts, fileunit, i
    
    file = 'exponential_integrals_e1z.csv'

    filename = 'test/test_points/' // file
    call read_test_points(filename, ref_x, ref_y, npts)
    allocate(specfun_y(npts))
    allocate(specfun_ic(npts))

    do i = 1, npts
      specfun_y(i) = e1(ref_x(i))
      specfun_ic(i) = isclose(ref_y(i), specfun_y(i))
    end do
 
    filename = 'test/test_specfun/' // file
    call write_test_points(filename, ref_x, ref_y, specfun_y, specfun_ic)

    call check(error, all(specfun_ic))
    if (allocated(error)) return
  end subroutine test_e1z

end module test_exponential_integrals
