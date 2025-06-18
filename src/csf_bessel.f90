!  ┏┓┏┓┏┓  Licensed under the MIT License
!  ┃ ┗┓┣   Copyright (c) 2025 Rodrigo Castro 
!  ┗┛┗┛┻   https://github.com/rodpcastro/colspecf

module csf_bessel
!* # Bessel functions
! Bessel functions.
!
! Procedures:
!
! - `j0x`: Bessel function of the first kind of order zero \(J_0(x)\)
! - `j1x`: Bessel function of the first kind of order one \(J_1(x)\)
! - `y0x`: Bessel function of the second kind of order zero \(Y_0(x)\)
! - `y1x`: Bessel function of the second kind of order one \(Y_0(x)\)
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
  public :: j0x, j1x, y0x, y1x

contains

  real(wp) function j0x(x)
    !! Bessel function of the first kind of order zero \(J_0(x)\).

    real(wp), intent(in) :: x

    if (abs(x) <= 1.0e-8_wp) then
      j0x = 1.0_wp
    else if (abs(x) >= 1.0e32_wp) then
      j0x = 0.0_wp
    else
      call caljy0(x, j0x, 0)
    end if
  end function j0x

  real(wp) function j1x(x)
    !! Bessel function of the first kind of order one \(J_1(x)\).

    real(wp), intent(in) :: x

    if (abs(x) < 1.0e-323_wp) then
      j1x = 0.0_wp
    else if (abs(x) >= 1.0e32_wp) then
      j1x = 0.0_wp
    else
      call caljy1(x, j1x, 0)
    end if
  end function j1x

  real(wp) function y0x(x)
    !! Bessel function of the second kind of order zero \(Y_0(x)\).

    real(wp), intent(in) :: x

    if (x < 0) then
      y0x = nan()
    else if (x < 1.0e-323_wp) then
      y0x = ninf()
    else if (x >= 1.0e32_wp) then
      y0x = 0.0_wp
    else
      call caljy0(x, y0x, 1)
    end if
  end function y0x

  real(wp) function y1x(x)
    !! Bessel function of the second kind of order one \(Y_1(x)\).

    real(wp), intent(in) :: x

    if (x < 0) then
      y1x = nan()
    else if (x < 1.0e-308_wp) then
      y1x = ninf()
    else if (x >= 1.0e32_wp) then
      y1x = 0.0_wp
    else
      call caljy1(x, y1x, 1)
    end if
  end function y1x

end module csf_bessel
