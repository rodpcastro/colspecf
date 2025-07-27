module specfun_evaluation
! Evaluation of special functions for test points.

  use csf_kinds, only: wp
  use csf_numerror, only: isclose
  use readwrite, only: read_test_points, write_test_points

  implicit none
  private
  public :: eval_write

  abstract interface
    ! Dummy function of single real variable.
    function specfunx(x) result(y)
      import :: wp
      real(wp), intent(in) :: x
      real(wp) :: y
    end function specfunx
  end interface

  abstract interface
    ! Dummy function of single complex variable.
    function specfunz(z) result(w)
      import :: wp
      complex(wp), intent(in) :: z
      complex(wp) :: w
    end function specfunz
  end interface

  interface eval_write
    ! Evaluates real or complex function for test points and write results.
    module procedure evalx_write, evalz_write
  end interface eval_write

contains

  subroutine evalx_write(fx, file, specfun_ic)
    ! Evaluates real function for test points and write results.

    procedure(specfunx) :: fx
    character(*), intent(in) :: file
    logical, intent(out), allocatable :: specfun_ic(:)
    real(wp), allocatable :: ref_x(:), ref_y(:), specfun_y(:)
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
    ! Evaluates complex function for test points and write results.

    procedure(specfunz) :: fz
    character(*), intent(in) :: file
    logical, intent(out), allocatable :: specfun_ic(:)
    complex(wp), allocatable :: ref_x(:), ref_y(:), specfun_y(:)
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
