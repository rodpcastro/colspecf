module example_hypergeometric
! Simple tests for hypergeometric functions.

  use csf_kinds, only: wp
  use csf_constants, only: nan, ninf, pinf
  use csf_numerror, only: isclose, eps_wp
  use csf_hypergeometric, only: hyp2f1

  implicit none
  private
  public :: example_hyp2f1
 
contains

  subroutine example_hyp2f1()
    complex(wp) :: a, b, c, z, hc, hm
    real(wp) :: abs_err, rel_err, err

    print '(a)', '---------------------------------------------'
    print '(a)', 'Gauss hypergeometric function hyp2f1(a,b,c,z)'
    print '(a)', '---------------------------------------------'

    a = (0.5_wp, 0.0_wp)
    b = (-0.5_wp, 0.0_wp)
    c = (1.5_wp, 0.0_wp)
    z = (1.5_wp, 0.0_wp)
    hc = hyp2f1(a, b, c, z)
    hm = (0.6412749150809320536_wp, 0.08473048557705339066_wp)
    abs_err = abs(hm - hc)
    rel_err = abs_err / abs(hm)
    err = min(abs_err, rel_err)

    print '(a, es22.15)', 'err_hyp2f1(0.5, -0.5, 1.5, 1.5) = ', err

    a = (0.5_wp, 0.0_wp)
    b = (-5.0_wp, 0.0_wp)
    c = (1.5_wp, 0.0_wp)
    z = (1.5_wp, 0.0_wp)
    hc = hyp2f1(a, b, c, z)
    hm = (0.3007305194805194801_wp, 0.0_wp)
    abs_err = abs(hm - hc)
    rel_err = abs_err / abs(hm)
    err = min(abs_err, rel_err)

    print '(a, es22.15)', 'err_hyp2f1(0.5, -5.0, 1.5, 1.5) = ', err
  end subroutine example_hyp2f1

end module example_hypergeometric
