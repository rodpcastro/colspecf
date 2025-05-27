module acm_exponential_integral
! Exponential integrals Ei, E1 and En.
!
! Author
! ------
! Rodrigo Castro (GitHub: rodpcastro)
!
! History
! -------
! 27-05-2025 - Rodrigo Castro - Original code
!
! References
! ----------
! [1] Kathleen A. Paciorek. 1970. Algorithm 385: Exponential integral Ei(x).
!     Commun. ACM 13, 7 (July 1970), 446–447. https://doi.org/10.1145/362686.362696
! [2] Donald E. Amos. 1990. Algorithms 683: a portable FORTRAN subroutine for
!     exponential integrals of a complex argument. ACM Trans. Math. Softw. 16,
!     2 (June 1990), 178–182. https://doi.org/10.1145/78928.78934

  use, intrinsic :: iso_fortran_env, only: int32, real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, &
                                           ieee_positive_inf, &
                                           ieee_negative_inf
  use numerror, only: eps64
  use calgo_385, only: dei
  use calgo_556, only: expint
  use calgo_683, only: cexint

  implicit none
  private
  public :: ei, e1, e1x, e1z, enz

  interface e1
    module procedure e1x, e1z
  end interface e1

contains

  real(real64) function ei(x)
    ! Exponential integral Ei(x).
    !
    ! Parameters
    ! ----------
    ! x : real(real64)
    !   Real number.
    !    
    ! Returns
    ! -------
    ! ei : real(real64) 
    !   Exponential integral Ei(x).

    real(real64), intent(in) :: x

    if (x == 0.0d0) then
      ei = ieee_value(1.0_real64, ieee_negative_inf)
    else
      ei = dei(x)
    end if
  end function ei

  real(real64) function e1x(x)
    ! Exponential integral E1(x).
    !
    ! Parameters
    ! ----------
    ! x : real(real64)
    !   Real number ≥ 0.
    !    
    ! Returns
    ! -------
    ! e1x : real(real64) 
    !   Exponential integral E1(x).

    real(real64), intent(in) :: x
    complex(real64) :: z, e1z

    if (x == 0.0d0) then
      e1x = ieee_value(1.0_real64, ieee_positive_inf)
    else
      z = cmplx(x, 0.0d0)
      e1z = enz(1, z)
      e1x = e1z%re
    end if
  end function e1x

  complex(real64) function e1z(z)
    ! Exponential integral E1(z).
    !
    ! Parameters
    ! ----------
    ! z : complex(real64)
    !   Complex number.
    !    
    ! Returns
    ! -------
    ! e1z : complex(real64) 
    !   Exponential integral E1(z).

    complex(real64), intent(in) :: z
    
    if (abs(z) == 0.0d0) then
      e1z = cmplx(ieee_value(1.0_real64, ieee_positive_inf), 0.0_real64)
    else
      e1z = enz(1, z)
    end if
  end function e1z

  complex(real64) function enz(n, z)
    ! Exponential integral En(z).
    !
    ! Parameters
    ! ----------
    ! n : integer(int32)
    !   Exponential integral order.
    ! z : complex(real64)
    !   Complex number.
    !    
    ! Returns
    ! -------
    ! enz : complex(real64) 
    !   Exponential integral En(z).

    integer(int32), intent(in) :: n
    complex(real64), intent(in) :: z
    real(real64), parameter :: tol = eps64
    integer(int32), parameter :: m = 1, kode = 1
    complex(real64) :: cy(m)
    integer(int32) :: ierr
    
    call cexint(z, n, kode, tol, m, cy, ierr)  ! TODO: Write error messages for each ierr value.
    enz = cy(1)
  end function enz

end module acm_exponential_integral
