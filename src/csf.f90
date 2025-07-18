!  ┏┓┏┓┏┓  Licensed under the MIT License
!  ┃ ┗┓┣   Copyright (c) 2025 Rodrigo Castro 
!  ┗┛┗┛┻   https://github.com/rodpcastro/colspecf
!          ASCII Art (Font Tmplr) by https://patorjk.com

module csf
!* # ColSpecF
! Collected Special Functions.
!
! In the following list of procedures, \(x \in \mathbb{R}\) and \(z \in \mathbb{C}\).
!
! Procedures:
!
! - `ei`: Exponential integral \(\mathrm{Ei}(x)\)
! - `e1`: Exponential integral \(\mathrm{E}_1(x)\) or \(\mathrm{E}_1(z)\)
! - `j0x`: Bessel function of the first kind of order zero \(J_0(x)\)
! - `j1x`: Bessel function of the first kind of order one \(J_1(x)\)
! - `y0x`: Bessel function of the second kind of order zero \(Y_0(x)\)
! - `y1x`: Bessel function of the second kind of order one \(Y_0(x)\)
! - `hyp2f1`: Gauss hypergeometric function \({}_2F_1(a, b; c; z)\)
!*

  use csf_exponential_integral
  use csf_bessel
  use csf_hypergeometric

  implicit none
  public

end module csf
