module specfun_evaluation
! Evaluation of special function for test points.
!
! Author
! ------
! Rodrigo Castro (GitHub: rodpcastro)
!
! History
! -------
! 25-05-2025 - Rodrigo Castro - Original code

  use, intrinsic :: iso_fortran_env, only: real64
  use readwrite, only: read_test_points, write_test_points
  use numerror, only: isclose

  implicit none
  private
  public :: eval_write

  abstract interface
    ! Dummy function of single real variable.

    function specfunx(x) result(y)
      import :: real64
      real(real64), intent(in) :: x
      real(real64) :: y
    end function specfunx
  end interface

  abstract interface
    ! Dummy function of single complex variable.

    function specfunz(z) result(w)
      import :: real64
      complex(real64), intent(in) :: z
      complex(real64) :: w
    end function specfunz
  end interface

  interface eval_write
    ! Evaluates special function for test points and write results.
    ! Works with both real and complex data.

    module procedure evalx_write, evalz_write
  end interface eval_write

contains

  subroutine evalx_write(fx, file, specfun_ic)
    procedure(specfunx) :: fx
    character(*), intent(in) :: file
    logical, intent(out), allocatable :: specfun_ic(:)
    real(real64), allocatable :: ref_x(:), ref_y(:), specfun_y(:)
    integer :: npts, fileunit, i
    character(len=100) :: filename

    filename = 'test/test_points/' // file
    call read_test_points(filename, ref_x, ref_y, npts)
    allocate(specfun_y(npts))
    allocate(specfun_ic(npts))

    do i = 1, npts
      specfun_y(i) = fx(ref_x(i))
      specfun_ic(i) = isclose(ref_y(i), specfun_y(i))
    end do
 
    filename = 'test/test_specfun/' // file
    call write_test_points(filename, ref_x, ref_y, specfun_y, specfun_ic)
  end subroutine evalx_write

  subroutine evalz_write(fz, file, specfun_ic)
    procedure(specfunz) :: fz
    character(*), intent(in) :: file
    logical, intent(out), allocatable :: specfun_ic(:)
    complex(real64), allocatable :: ref_x(:), ref_y(:), specfun_y(:)
    integer :: npts, fileunit, i
    character(len=100) :: filename

    filename = 'test/test_points/' // file
    call read_test_points(filename, ref_x, ref_y, npts)
    allocate(specfun_y(npts))
    allocate(specfun_ic(npts))

    do i = 1, npts
      specfun_y(i) = fz(ref_x(i))
      specfun_ic(i) = isclose(ref_y(i), specfun_y(i))
    end do
 
    filename = 'test/test_specfun/' // file
    call write_test_points(filename, ref_x, ref_y, specfun_y, specfun_ic)
  end subroutine evalz_write

end module specfun_evaluation
