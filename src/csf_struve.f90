!  ┏┓┏┓┏┓  Licensed under the MIT License
!  ┃ ┗┓┣   Copyright (c) 2025 Rodrigo Castro 
!  ┗┛┗┛┻   https://github.com/rodpcastro/colspecf

module csf_struve
!* # Struve functions
! Struve functions.
!
! Procedures:
!
! - `struveh0`: Struve function \(\mathbf{H}_0(x)\)
! - `struveh1`: Struve function \(\mathbf{H}_1(x)\)
!
! ## References
! 1. Allan J. MacLeod. 1996. Algorithm 757: MISCFUN, a software package to compute
!    uncommon special functions. ACM Trans. Math. Softw. 22, 3 (Sept. 1996), 288–301.
!*   <https://doi.org/10.1145/232826.232846>

  use csf_kinds, only: wp
  use calgo_757, only: strvh0, strvh1

  implicit none
  private
  public :: struveh0, struveh1

contains

  real(wp) function struveh0(x)
    !! Struve function \(\mathbf{H}_0(x)\).
    !
    !! \(x \in \mathbb{R}\)

    real(wp), intent(in) :: x

    struveh0 = strvh0(x)
  end function struveh0

  real(wp) function struveh1(x)
    !! Struve function \(\mathbf{H}_1(x)\).
    !
    !! \(x \in \mathbb{R}\)

    real(wp), intent(in) :: x

    struveh1 = strvh1(x)
  end function struveh1

end module csf_struve
