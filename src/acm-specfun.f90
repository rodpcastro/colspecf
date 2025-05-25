module acm_specfun
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, acm-specfun!"
  end subroutine say_hello
end module acm_specfun
