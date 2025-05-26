module acm_exponential_integral
! Exponential integrals Ei, E1 and En.
!
! Author
! ------
! Rodrigo Castro (GitHub: rodpcastro)
!
! History
! -------
! 26-05-2025 - Rodrigo Castro - Original code
!
! References
! ----------
! [1] Shanjie Zhang, Jianming Jin. 1996. Computation of Special Functions.
!     Wiley, New York, NY.

  use, intrinsic :: iso_fortran_env, only: int32, real64
  use numerror, only: eps64
  use calgo_683, only: cexint

  implicit none
  private
  public :: e1z, enz

contains

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
    
    e1z = enz(1, z)
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
