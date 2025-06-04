!  ┏┓┏┓┏┓  Licensed under the MIT License
!  ┃ ┗┓┣   Copyright (c) 2025 Rodrigo Castro 
!  ┗┛┗┛┻   https://github.com/rodpcastro/colspecf

module csf_template
!* Module template. The space after * is important for Ford `docmark_alt`.
!
! References
! ----------
! [1] Author. Year. Title.
!*    URL

  use csf_kinds, only: i1, wp

  implicit none
  private
  public :: phi, fibonacci

  real(wp), parameter :: phi = 1.61803398874989484820458683436563811e+0_wp  !! Golden ratio φ

contains

  pure function fibonacci(n) result(fib)
    !! n first terms of the Fibonacci sequence.

    integer(i1), intent(in) :: n  !! number of terms
    integer(i1) :: fib(n)
    integer(i1) :: i

    do i = 1, n
      if (i == 1) then
        fib(i) = 0
      elseif (i == 2) then
        fib(i) = 1
      else
        fib(i) = fib(i-1) + fib(i-2)
      end if
    end do
  end function fibonacci

end module csf_template
