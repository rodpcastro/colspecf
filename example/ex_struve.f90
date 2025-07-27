module example_struve
! Simple tests for Struve functions.

  use csf_kinds, only: wp
  use csf_constants, only: ninf, pinf
  use csf_numerror, only: isclose
  use csf_struve, only: struveh0, struveh1

  implicit none
  private
  public :: example_struveh0, example_struveh1
 
contains

  subroutine example_struveh0()
    real(wp) :: h0_m

    h0_m = 0.5686566270482879548_wp
    if (.not. isclose(struveh0(1.0_wp), h0_m)) error stop 'error: example struveh0'

    print '(a)', '---------------------'
    print '(a)', 'Struve function H0(x)'
    print '(a)', '---------------------'
    print '(a, es22.15)', 'struveh0(0.0) = ', struveh0(0.0_wp)
    print '(a, es22.15)', 'struveh0(1.0) = ', struveh0(1.0_wp)
    print '(a, es22.15)', 'struveh0(-1.0) = ', struveh0(-1.0_wp)
    print '(a, sp, g0)', 'struveh0(+Inf) = ', struveh0(pinf())
    print '(a, sp, g0)', 'struveh0(-Inf) = ', struveh0(ninf())
  end subroutine example_struveh0

  subroutine example_struveh1()
    real(wp) :: h1_m

    h1_m = 0.1984573362019444003_wp
    if (.not. isclose(struveh1(1.0_wp), h1_m)) error stop 'error: example struveh1'

    print '(a)', '---------------------'
    print '(a)', 'Struve function H1(x)'
    print '(a)', '---------------------'
    print '(a, es22.15)', 'struveh1(0.0) = ', struveh1(0.0_wp)
    print '(a, es22.15)', 'struveh1(1.0) = ', struveh1(1.0_wp)
    print '(a, es22.15)', 'struveh1(-1.0) = ', struveh1(-1.0_wp)
    print '(a, sp, g0)', 'struveh1(+Inf) = ', struveh1(pinf())
    print '(a, sp, g0)', 'struveh1(-Inf) = ', struveh1(ninf())
  end subroutine example_struveh1

end module example_struve
