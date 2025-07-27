program examples
! Simple tests for internal procedures and special functions.

  use example_exponential_integral
  use example_bessel
  use example_hypergeometric
  use example_struve

  implicit none
  
  call example_ei()
  call example_e1x()
  call example_e1z()
  call example_enz()
  call example_besselj0()
  call example_besselj1()
  call example_bessely0()
  call example_bessely1()
  call example_hyp2f1()
  call example_struveh0()
  call example_struveh1()

end program examples
