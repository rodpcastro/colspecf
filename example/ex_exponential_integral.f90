module example_exponential_integral
! Simple tests for exponential integrals Ei and E1.

  use csf_kinds, only: wp
  use csf_constants, only: ninf, pinf
  use csf_exponential_integral, only: ei, e1, enz

  implicit none
  private
  public :: example_ei, example_e1x, example_e1z, example_enz
 
contains

  subroutine example_ei()
    real(wp) :: ei_t1, ei_t2, ei_root
    ei_root = 0.3725074107813666_wp

    ei_t1 = abs(-0.30266853926582593_wp - ei(0.3_wp)) / abs(-0.30266853926582593_wp)
    ei_t2 = abs(1.058563689713169e+20_wp - ei(50.0_wp)) / abs(1.058563689713169e+20_wp)

    print '(a)', '--------------------------'
    print '(a)', 'Exponential Integral Ei(x)'
    print '(a)', '--------------------------'
    print '(a, es22.15)', 'ei_t1 = ', ei_t1
    print '(a, es22.15)', 'ei_t2 = ', ei_t2
    print '(a, es22.15)', '-ei(1.0) = ', -ei(1.0_wp)
    print '(a, es22.15)', 'ei(-1.5) = ', ei(-1.5_wp)
    print '(a, es22.15)', 'ei(root) = ', ei(ei_root)
    print '(a, sp, g0)', 'ei(0.0) = ', ei(0.0_wp)
    print '(a, sp, g0)', 'ei(+Inf) = ', ei(pinf())
    print '(a, sp, g0)', 'ei(+Inf) - ei(709.0) = ', ei(pinf()) - ei(709.0_wp)
    print '(a, sp, g0)', 'ei(+Inf) - ei(710.0) = ', ei(pinf()) - ei(710.0_wp)
    print '(a, sp, g0)', 'ei(709.0) = ', ei(709.0_wp)
    print '(a, sp, g0)', 'ei(710.0) = ', ei(710.0_wp)
    print '(a, sp, g0)', 'ei(-Inf) = ', ei(ninf())
    print '(a, es22.15)', 'ei(-738.0) = ', ei(-738.0_wp)
    print '(a, es22.15)', 'ei(-739.0) = ', ei(-739.0_wp)
  end subroutine example_ei

  subroutine example_e1x()
    real(wp) :: e1x_t1, e1x_t2

    e1x_t1 = abs(0.55977359477616084_wp - e1(0.5_wp)) / abs(0.55977359477616084_wp)
    e1x_t2 = abs(0.10001958240663265_wp - e1(1.5_wp)) / abs(0.10001958240663265_wp)

    print '(a)', '--------------------------'
    print '(a)', 'Exponential Integral E1(x)'
    print '(a)', '--------------------------'
    print '(a, es22.15)', 'e1x_t1 = ', e1x_t1
    print '(a, es22.15)', 'e1x_t2 = ', e1x_t2
    print '(a, es22.15)', '-e1x(1.5) = ', -e1(1.5_wp)
    print '(a, sp, g0)' , 'e1x(0.0) = ', e1(0.0_wp)
    print '(a, sp, g0)', 'e1x(738.0) = ', e1(738.0_wp)
    print '(a, sp, g0)', 'e1x(739.0) = ', e1(739.0_wp)
    print '(a, sp, g0)', 'e1x(+Inf) = ', e1(pinf())
  end subroutine example_e1x

  subroutine example_e1z()
    real(wp) :: e1z_t1, e1z_t2, e1z_t3
    complex(wp) :: zre_ninf, zre_pinf, zim_ninf, zim_pinf, z_ppinf, z_pninf

    zre_ninf = cmplx(ninf(), 0.0_wp, kind=wp)
    zre_pinf = cmplx(pinf(), 0.0_wp, kind=wp)
    zim_ninf = cmplx(0.0_wp, ninf(), kind=wp)
    zim_pinf = cmplx(0.0_wp, pinf(), kind=wp)
    z_ppinf = cmplx(pinf(), pinf(), kind=wp)
    z_pninf = cmplx(pinf(), ninf(), kind=wp)

    e1z_t1 = abs((-0.014529959529202443_wp, -0.015866824826503003_wp) - e1((2.5_wp, 1.8_wp))) / &
             abs((-0.014529959529202443_wp, -0.015866824826503003_wp))
    e1z_t2 = abs((0.75128638206377318_wp, -57.45565336638456_wp) - e1((-6.4_wp, -8.9_wp))) / &
             abs((0.75128638206377318_wp, -57.45565336638456_wp))
    e1z_t3 = abs((-4.2432089533906728e-6_wp, -7.7188970166104144e-6_wp) - e1((9.1_wp, 7.7_wp))) / &
             abs((-4.2432089533906728e-6_wp, -7.7188970166104144e-6_wp))

    print '(a)', '--------------------------'
    print '(a)', 'Exponential Integral E1(z)'
    print '(a)', '--------------------------'
    print '(a, es22.15)', 'e1z_t1 = ', e1z_t1
    print '(a, es22.15)', 'e1z_t2 = ', e1z_t2
    print '(a, es22.15)', 'e1z_t3 = ', e1z_t3
    print '(a, 2(es22.15, 1x))', 'e1z((-1.0, 0.0)) = ', e1((-1.0_wp, 0.0_wp))
    print '(a, sp, g0, es22.15)', 'e1z((0.0, 0.0)) = ', e1((0.0_wp, 0.0_wp))
    print '(a, sp, g0, es22.15)', 'e1z((-Inf, 0.0) = ', e1(zre_ninf)
    print '(a, sp, g0, es22.15)', 'e1z((+Inf, 0.0)) = ', e1(zre_pinf)
    print '(a, sp, g0, es22.15)', 'e1z((0.0, -Inf)) = ', e1(zim_ninf)
    print '(a, sp, g0, es22.15)', 'e1z((0.0, +Inf)) = ', e1(zim_pinf)
    print '(a, sp, g0, es22.15)', 'e1z((+Inf, +Inf)) = ', e1(z_ppinf)
    print '(a, sp, g0, es22.15)', 'e1z((+Inf, -Inf)) = ', e1(z_pninf)
    print *, '------------------------------'
    print '(a, 2(es22.15, 1x))', 'e1z((-700.0, 0.0)) = ', e1((-700.0_wp, 0.0_wp))
    print '(a, 2(es22.15, 1x))', 'e1z((-701.0, 0.0)) = ', e1((-701.0_wp, 0.0_wp))
    print '(a, 2(es22.15, 1x))', 'e1z((738.0, 0.0)) = ', e1((738.0_wp, 0.0_wp))
    print '(a, 2(es22.15, 1x))', 'e1z((739.0, 0.0)) = ', e1((739.0_wp, 0.0_wp))
    print '(a, 2(es22.15, 1x))', 'e1z((0.0, -1.0e15)) = ', e1((0.0_wp, -1.0e15_wp))
    print '(a, 2(es22.15, 1x))', 'e1z((0.0, -1.0e16)) = ', e1((0.0_wp, -1.0e16_wp))
    print '(a, 2(es22.15, 1x))', 'e1z((0.0, +1.0e15)) = ', e1((0.0_wp, +1.0e15_wp))
    print '(a, 2(es22.15, 1x))', 'e1z((0.0, +1.0e16)) = ', e1((0.0_wp, +1.0e16_wp))
    print '(a, 2(es22.15, 1x))', 'e1z((+1.0e15, +1.0e15)) = ', e1((1.0e15_wp, 1.0e15_wp))
    print '(a, 2(es22.15, 1x))', 'e1z((+1.0e16, +1.0e16)) = ', e1((1.0e16_wp, 1.0e16_wp))
    print '(a, 2(es22.15, 1x))', 'e1z((+1.0e15, -1.0e15)) = ', e1((1.0e15_wp, -1.0e15_wp))
    print '(a, 2(es22.15, 1x))', 'e1z((+1.0e16, -1.0e16)) = ', e1((1.0e16_wp, -1.0e16_wp))
  end subroutine example_e1z

  subroutine example_enz()
    real(wp) :: e2zt, e3zt, e4zt
    complex(wp) :: z, e2z, e3z, e4z

    z = (7.8_wp, -4.1_wp)
    e2z = enz(2, z)
    e3z = enz(3, z)
    e4z = enz(4, z)
    e2zt = abs((-7.957699432525390302e-6, -0.00003818901915888152331) - e2z)
    e3zt = abs((-8.440221762387431442e-6, -0.00003501448656784751492) - e3z)
    e4zt = abs((-8.710784325030124611e-6, -0.00003225622296606488706) - e4z)

    print '(a)', '--------------------------'
    print '(a)', 'Exponential Integral En(x)'
    print '(a)', '--------------------------'
    print '(a, es22.15)', 'e2zt = ', e2zt
    print '(a, es22.15)', 'e3zt = ', e3zt
    print '(a, es22.15)', 'e4zt = ', e4zt
  end subroutine example_enz

end module example_exponential_integral
