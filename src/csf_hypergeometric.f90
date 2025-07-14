!  ┏┓┏┓┏┓  Licensed under the MIT License
!  ┃ ┗┓┣   Copyright (c) 2025 Rodrigo Castro 
!  ┗┛┗┛┻   https://github.com/rodpcastro/colspecf

module csf_hypergeometric
!* # Hypergeometric
! Hypergeometric functions.
!
! Procedures:
!
! - `hyp2f1`: Gauss hypergeometric function \({}_2F_1(a, b; c; z)\)
!
! ## References
! 1. N. Michel and M. V. Stoitsov. 2008. Fast computation of the Gauss hypergeometric 
!    function with all its parameters complex with application to the Pöschl-Teller-
!    Ginocchio potential wave functions. Computer Physics Communications 178, 7 (April
!*   2008), 535–551. <https://doi.org/10.1016/J.CPC.2007.11.007>

  use csf_kinds, only: wp
  use cpc_michel, only: hyp_2f1

  implicit none
  private
  public :: hyp2f1

contains

  complex(wp) function hyp2f1(a, b, c, z)
    !! Gauss hypergeometric function \({}_2F_1(a, b; c; z)\)

    complex(wp), intent(in) :: a, b, c, z

    hyp2f1 = hyp_2f1(a, b, c, z)

  end function hyp2f1

end module csf_hypergeometric
