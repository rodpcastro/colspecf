module example_bessel
! Simple tests for Bessel functions.

  use csf_kinds, only: wp
  use csf_constants, only: nan, ninf, pinf
  use csf_numerror, only: isclose, eps_wp
  use csf_bessel, only: besselj0, besselj1, bessely0, bessely1

  implicit none
  private
  public :: example_besselj0, example_besselj1, example_bessely0, example_bessely1
 
contains

  subroutine example_besselj0()
    real(wp) :: besselj0_ref

    besselj0_ref = 0.7651976865579665514_wp
    if (.not. isclose(besselj0(1.0_wp), besselj0_ref)) then
      error stop 'error: example besselj0'
    end if

    print '(a)', '---------------------'
    print '(a)', 'Bessel function J0(x)'
    print '(a)', '---------------------'
    print '(a, es22.15)', 'besselj0(0.0) = ', besselj0(0.0_wp)
    print '(a, es22.15)', 'besselj0(1.0) = ', besselj0(1.0_wp)
    print '(a, es22.15)', 'besselj0(-1.0) = ', besselj0(-1.0_wp)
    print '(a, es22.15)', 'besselj0(1.0e32) = ', besselj0(1.0e31_wp)
    print '(a, es22.15)', 'besselj0(1.0e33) = ', besselj0(1.0e32_wp)
    print '(a, es22.15)', 'besselj0(1.0e-7) = ', besselj0(1.0e-7_wp)
    print '(a, es22.15)', 'besselj0(1.0e-8) = ', besselj0(1.0e-8_wp)
    print '(a, sp, g0)', 'besselj0(+Inf) = ', besselj0(pinf())
    print '(a, sp, g0)', 'besselj0(-Inf) = ', besselj0(ninf())
  end subroutine example_besselj0

  subroutine example_besselj1()
    real(wp) :: besselj1_ref

    besselj1_ref = 0.440050585744933516_wp
    if (.not. isclose(besselj1(1.0_wp), besselj1_ref)) then
      error stop 'error: example besselj1'
    end if

    print '(a)', '---------------------'
    print '(a)', 'Bessel function J1(x)'
    print '(a)', '---------------------'
    print '(a, es22.15)', 'besselj1(0.0) = ', besselj1(0.0_wp)
    print '(a, es22.15)', 'besselj1(1.0) = ', besselj1(1.0_wp)
    print '(a, es22.15)', 'besselj1(-1.0) = ', besselj1(-1.0_wp)
    print '(a, es22.15)', 'besselj1(1.0e32) = ', besselj1(1.0e31_wp)
    print '(a, es22.15)', 'besselj1(1.0e33) = ', besselj1(1.0e32_wp)
    print '(a, es22.15)', 'besselj1(1.0e-who) = ', besselj1(1.0e-323_wp)
    print '(a, es22.15)', 'besselj1(1.0e-who) = ', besselj1(1.0e-324_wp)
  end subroutine example_besselj1

  subroutine example_bessely0()
    real(wp) :: bessely0_ref

    bessely0_ref = 0.08825696421567695744_wp
    if (.not. isclose(bessely0(1.0_wp), bessely0_ref)) then
      error stop 'error: example bessely0'
    end if

    print '(a)', '---------------------'
    print '(a)', 'Bessel function Y0(x)'
    print '(a)', '---------------------'
    print '(a, es22.15)', 'bessely0(0.0) = ', bessely0(0.0_wp)
    print '(a, es22.15)', 'bessely0(1.0) = ', bessely0(1.0_wp)
    print '(a, es22.15)', 'bessely0(-1.0) = ', bessely0(-1.0_wp)
    print '(a, es22.15)', 'bessely0(1.0e32) = ', bessely0(1.0e31_wp)
    print '(a, es22.15)', 'bessely0(1.0e33) = ', bessely0(1.0e32_wp)
    print '(a, es22.15)', 'bessely0(1.0e-323_wp) = ', bessely0(1.0e-323_wp)
    print '(a, es22.15)', 'bessely0(1.0e-324_wp) = ', bessely0(1.0e-324_wp)
  end subroutine example_bessely0

  subroutine example_bessely1()
    real(wp) :: bessely1_ref

    bessely1_ref = -0.7812128213002887123_wp

    if (.not. isclose(bessely1(1.0_wp), bessely1_ref)) then
      error stop 'error: example bessely1'
    end if

    print '(a)', '---------------------'
    print '(a)', 'Bessel function Y1(x)'
    print '(a)', '---------------------'
    print '(a, es22.15)', 'bessely1(0.0) = ', bessely1(0.0_wp)
    print '(a, es22.15)', 'bessely1(1.0) = ', bessely1(1.0_wp)
    print '(a, es22.15)', 'bessely1(-1.0) = ', bessely1(-1.0_wp)
    print '(a, es22.15)', 'bessely1(1.0e32) = ', bessely1(1.0e31_wp)
    print '(a, es22.15)', 'bessely1(1.0e33) = ', bessely1(1.0e32_wp)
    print '(a, es22.15)', 'bessely1(1.0e-308_wp) = ', bessely1(1.0e-308_wp)
    print '(a, es22.15)', 'bessely1(1.0e-309_wp) = ', bessely1(1.0e-309_wp)
  end subroutine example_bessely1

end module example_bessel
