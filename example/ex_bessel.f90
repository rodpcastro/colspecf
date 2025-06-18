module example_bessel
! Simple tests for Bessel functions.

  use csf_kinds, only: wp
  use csf_constants, only: nan, ninf, pinf
  use csf_numerror, only: isclose, eps_wp
  use csf_bessel, only: j0x, j1x, y0x, y1x

  implicit none
  private
  public :: example_j0x, example_j1x, example_y0x, example_y1x
 
contains

  subroutine example_j0x()
    real(wp) :: j0x_m

    j0x_m = 0.7651976865579665514_wp
    if (.not. isclose(j0x(1.0_wp), j0x_m)) error stop 'error: example j0x'

    print '(a)', '---------------------'
    print '(a)', 'Bessel function J0(x)'
    print '(a)', '---------------------'
    print '(a, es22.15)', 'j0x(0.0) = ', j0x(0.0_wp)
    print '(a, es22.15)', 'j0x(1.0) = ', j0x(1.0_wp)
    print '(a, es22.15)', 'j0x(-1.0) = ', j0x(-1.0_wp)
    print '(a, es22.15)', 'j0x(1.0e32) = ', j0x(1.0e31_wp)
    print '(a, es22.15)', 'j0x(1.0e33) = ', j0x(1.0e32_wp)
    print '(a, es22.15)', 'j0x(1.0e-7) = ', j0x(1.0e-7_wp)
    print '(a, es22.15)', 'j0x(1.0e-8) = ', j0x(1.0e-8_wp)
    print '(a, sp, g0)', 'j0x(+Inf) = ', j0x(pinf())
    print '(a, sp, g0)', 'j0x(-Inf) = ', j0x(ninf())
  end subroutine example_j0x

  subroutine example_j1x()
    real(wp) :: j1x_m

    j1x_m = 0.440050585744933516_wp
    if (.not. isclose(j1x(1.0_wp), j1x_m)) error stop 'error: example j1x'

    print '(a)', '---------------------'
    print '(a)', 'Bessel function J1(x)'
    print '(a)', '---------------------'
    print '(a, es22.15)', 'j1x(0.0) = ', j1x(0.0_wp)
    print '(a, es22.15)', 'j1x(1.0) = ', j1x(1.0_wp)
    print '(a, es22.15)', 'j1x(-1.0) = ', j1x(-1.0_wp)
    print '(a, es22.15)', 'j1x(1.0e32) = ', j1x(1.0e31_wp)
    print '(a, es22.15)', 'j1x(1.0e33) = ', j1x(1.0e32_wp)
    print '(a, es22.15)', 'j1x(1.0e-who) = ', j1x(1.0e-323_wp)
    print '(a, es22.15)', 'j1x(1.0e-who) = ', j1x(1.0e-324_wp)
  end subroutine example_j1x

  subroutine example_y0x()
    real(wp) :: y0x_m

    y0x_m = 0.08825696421567695744_wp
    if (.not. isclose(y0x(1.0_wp), y0x_m)) error stop 'error: example y0x'

    print '(a)', '---------------------'
    print '(a)', 'Bessel function Y0(x)'
    print '(a)', '---------------------'
    print '(a, es22.15)', 'y0x(0.0) = ', y0x(0.0_wp)
    print '(a, es22.15)', 'y0x(1.0) = ', y0x(1.0_wp)
    print '(a, es22.15)', 'y0x(-1.0) = ', y0x(-1.0_wp)
    print '(a, es22.15)', 'y0x(1.0e32) = ', y0x(1.0e31_wp)
    print '(a, es22.15)', 'y0x(1.0e33) = ', y0x(1.0e32_wp)
    print '(a, es22.15)', 'y0x(1.0e-323_wp) = ', y0x(1.0e-323_wp)
    print '(a, es22.15)', 'y0x(1.0e-324_wp) = ', y0x(1.0e-324_wp)
  end subroutine example_y0x

  subroutine example_y1x()
    real(wp) :: y1x_m

    y1x_m = -0.7812128213002887123_wp
    if (.not. isclose(y1x(1.0_wp), y1x_m)) error stop 'error: example y1x'

    print '(a)', '---------------------'
    print '(a)', 'Bessel function Y1(x)'
    print '(a)', '---------------------'
    print '(a, es22.15)', 'y1x(0.0) = ', y1x(0.0_wp)
    print '(a, es22.15)', 'y1x(1.0) = ', y1x(1.0_wp)
    print '(a, es22.15)', 'y1x(-1.0) = ', y1x(-1.0_wp)
    print '(a, es22.15)', 'y1x(1.0e32) = ', y1x(1.0e31_wp)
    print '(a, es22.15)', 'y1x(1.0e33) = ', y1x(1.0e32_wp)
    print '(a, es22.15)', 'y1x(1.0e-308_wp) = ', y1x(1.0e-308_wp)
    print '(a, es22.15)', 'y1x(1.0e-309_wp) = ', y1x(1.0e-309_wp)
  end subroutine example_y1x

end module example_bessel
