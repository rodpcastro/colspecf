!  ┏┓┏┓┏┓  Licensed under the MIT License
!  ┃ ┗┓┣   Copyright (c) 2025 Rodrigo Castro 
!  ┗┛┗┛┻   https://github.com/rodpcastro/colspecf

module csf_bessel
!* # Bessel functions
! Bessel functions.
!
! Procedures:
!
! - `besselj0`: Bessel function of the first kind of order zero \(J_0(x)\)
! - `besselj1`: Bessel function of the first kind of order one \(J_1(x)\)
! - `bessely0`: Bessel function of the second kind of order zero \(Y_0(x)\)
! - `bessely1`: Bessel function of the second kind of order one \(Y_0(x)\)
!
! ## References
! 1. W. J. Cody. 1993. Algorithm 715: SPECFUN–a portable FORTRAN package of special
!    function routines and test drivers. ACM Trans. Math. Softw. 19, 1 (March 1993),
!*   22–30. <https://doi.org/10.1145/151271.151273>

  use csf_kinds, only: wp
  use csf_constants, only: nan, ninf, pinf
  use calgo_715, only: caljy0, caljy1

  implicit none
  private
  public :: besselj0, besselj1, bessely0, bessely1

contains

  real(wp) function besselj0(x)
    !! Bessel function of the first kind of order zero \(J_0(x)\).
    !
    !! \(x \in \mathbb{R}\)

    real(wp), intent(in) :: x

    call caljy0(x, besselj0, 0)
  end function besselj0

  real(wp) function besselj1(x)
    !! Bessel function of the first kind of order one \(J_1(x)\).
    !
    !! \(x \in \mathbb{R}\)

    real(wp), intent(in) :: x

    call caljy1(x, besselj1, 0)
  end function besselj1

  real(wp) function bessely0(x)
    !! Bessel function of the second kind of order zero \(Y_0(x)\).
    !
    !! \(\lbrace x \in \mathbb{R} \mid x \gt 0 \rbrace\)

    real(wp), intent(in) :: x

    if (x < 0) then
      bessely0 = nan()
    else
      call caljy0(x, bessely0, 1)
      if (bessely0 <= -huge(0.0_wp)) then
        bessely0 = ninf()
      end if
    end if
  end function bessely0

  real(wp) function bessely1(x)
    !! Bessel function of the second kind of order one \(Y_1(x)\).
    !
    !! \(\lbrace x \in \mathbb{R} \mid x \gt 0 \rbrace\)

    real(wp), intent(in) :: x

    if (x < 0) then
      bessely1 = nan()
    else
      call caljy1(x, bessely1, 1)
      if (bessely1 <= -huge(0.0_wp)) then
        bessely1 = ninf()
      end if
    end if
  end function bessely1

end module csf_bessel
