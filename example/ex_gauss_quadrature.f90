module example_gauss_quadrature
! Simple tests for Gauss-Legendre quadrature contained in CALGO 683.

  use csf_kinds, only: wp
  use csf_constants, only: pi
  use calgo_683, only: g8, gaus8

  implicit none
  private
  public :: example_g8, example_gaus8
 
contains

  real(wp) function funx(x)
    real(wp), intent(in) :: x
    funx = sin(x)
  end function funx

  subroutine example_g8()
    real(wp) :: a, b, x, h

    a = 0.0_wp
    b = pi
    x = (a + b) / 2.0_wp
    h = (b - a) / 2.0_wp

    print '(a)', '----------------------------------'
    print '(a)', 'Gauss-Legendre quadrature g8(f(x))'
    print '(a)', '----------------------------------'
    print '(a, es22.15)', 'Int(f(x)) = ', g8(funx, x, h)
  end subroutine example_g8

  subroutine example_gaus8()
    real(wp) :: a, b, err, ans
    integer :: ierr

    a = 0.0_wp
    b = pi
    err = 1.0e-8_wp
    call gaus8(funx, a, b, err, ans, ierr)

    print '(a)', '----------------------------------------------'
    print '(a)', 'Adaptive Gauss-Legendre quadrature gaus8(f(x))'
    print '(a)', '----------------------------------------------'
    print '(a, es22.15)', 'Int(f(x)) = ', ans
  end subroutine example_gaus8

end module example_gauss_quadrature

