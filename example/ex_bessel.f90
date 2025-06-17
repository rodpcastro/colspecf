module example_bessel
! Simple tests for Bessel functions.

  use csf_kinds, only: wp
  use csf_constants, only: ninf, pinf
  use csf_bessel, only: j0x, j1x, y0x, y1x

  implicit none
  private
  public :: example_j0x, example_j1x, example_y0x, example_y1x
 
contains

  subroutine example_j0x()
    real(wp) :: j0x_t1

    j0x_t1 = abs(0.7651976865579665514_wp - j0x(1.0_wp))

    print '(a)', '---------------------'
    print '(a)', 'Bessel function J0(x)'
    print '(a)', '---------------------'
    print '(a, es22.15)', 'j0x_t1 = ', j0x_t1
  end subroutine example_j0x

  subroutine example_j1x()
    real(wp) :: j1x_t1

    j1x_t1 = abs(0.440050585744933516_wp - j1x(1.0_wp))

    print '(a)', '---------------------'
    print '(a)', 'Bessel function J1(x)'
    print '(a)', '---------------------'
    print '(a, es22.15)', 'j1x_t1 = ', j1x_t1
  end subroutine example_j1x

  subroutine example_y0x()
    real(wp) :: y0x_t1

    y0x_t1 = abs(0.08825696421567695744_wp - y0x(1.0_wp))

    print '(a)', '---------------------'
    print '(a)', 'Bessel function Y0(x)'
    print '(a)', '---------------------'
    print '(a, es22.15)', 'y0x_t1 = ', y0x_t1
  end subroutine example_y0x

  subroutine example_y1x()
    real(wp) :: y1x_t1

    y1x_t1 = abs(-0.7812128213002887123_wp - y1x(1.0_wp))

    print '(a)', '---------------------'
    print '(a)', 'Bessel function Y1(x)'
    print '(a)', '---------------------'
    print '(a, es22.15)', 'y1x_t1 = ', y1x_t1
  end subroutine example_y1x

end module example_bessel
