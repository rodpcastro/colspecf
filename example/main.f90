program examples
! Simple tests for internal procedures and special functions.

  use example_gauss_quadrature
  use example_polygamma
  use example_exponential_integral
  use example_bessel
  use example_hypergeometric

  implicit none
  
  ! call example_g8()
  ! call example_gaus8()
  ! call example_psixn()
  call example_ei()
  call example_e1x()
  call example_e1z()
  ! call example_enz()
  call example_j0x()
  call example_j1x()
  call example_y0x()
  call example_y1x()
  call example_hyp2f1()

end program examples
