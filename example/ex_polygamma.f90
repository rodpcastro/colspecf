module example_polygamma
! Simple test for Digamma function contained in CALGO 683.

  use csf_kinds, only: wp
  use csf_constants, only: gm 
  use calgo_683, only: psixn

  implicit none
  private
  public :: example_psixn
 
contains

  subroutine example_psixn()
    print '(a)', '----------------'
    print '(a)', 'Digamma function'
    print '(a)', '----------------'
    print '(a, es22.15)', '-gamma =   ', -gm
    print '(a, es22.15)', 'psixn(1) = ', psixn(1)
  end subroutine example_psixn

end module example_polygamma
