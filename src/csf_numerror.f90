!  ┏┓┏┓┏┓  Licensed under the MIT License
!  ┃ ┗┓┣   Copyright (c) 2025 Rodrigo Castro 
!  ┗┛┗┛┻   https://github.com/rodpcastro/colspecf

module csf_numerror
!* # Numerical error
! Numerical error handling.
!
! Constants:
!
! - `eps_wp`: Working precision machine epsilon
! - `tol_wp`: Tolerance for approximations
!
! Procedures:
!
! - `ismall`: Evaluates the smallness of a variable compared to a reference value.
! - `isclose`: Evaluates the closeness between two real or complex variables.
!
! ## References
! 1. Christopher Barker. 2015. PEP 485 – A Function for testing
!*   approximate equality. <https://peps.python.org/pep-0485/>

  use csf_kinds, only: wp

  implicit none
  private
  public :: eps_wp, tol_wp, isclose, ismall

  real(wp), parameter :: eps_wp = epsilon(1.0_wp)  !! Working precision machine epsilon
  real(wp), parameter :: tol_wp = 1.0e-8_wp  !! Tolerance for approximations

  ! For reference:
  ! eps_sp = 1.19209290E-07
  ! eps_dp = 2.2204460492503131E-016

  interface isclose
    !! Evaluates the closeness between two real or complex variables.
    module procedure isclose_real, isclose_complex
  end interface isclose

contains

  pure logical function ismall(x, ref)
    !* Evaluates the smallness of a variable compared to a reference value.
    !
    ! `.true.` if `x` is small compared to `ref` according to `eps_wp`, and
    !* `.false.` otherwise.

    real(wp), intent(in) :: x
    real(wp), intent(in), optional :: ref  !! Reference value, default=`1.0`
    real(wp) :: ref_

    ref_ = 1.0_wp

    if (present(ref)) then
      ref_ = ref
    end if

    ismall = abs(x) < eps_wp * abs(ref_)
  end function ismall

  pure logical function isclose_real(a, b, rel_tol, abs_tol) 
    !* Evaluates the closeness between two real variables.
    ! 
    ! `.true.` if `a` and `b` are close to each other according to a tolerance, and
    ! `.false.` otherwise. Tolerance is defined by `rel_tol` and `abs_tol`, where
    !* `abs_tol` is used for comparing values close no zero.

    real(wp), intent(in) :: a, b
    real(wp), intent(in), optional :: rel_tol  !! Relative tolerance, default=`tol_wp`
    real(wp), intent(in), optional :: abs_tol  !! Absolute tolerace, default=`tol_wp`
    real(wp) :: rel_tol_, abs_tol_

    rel_tol_ = tol_wp
    abs_tol_ = tol_wp

    if (present(rel_tol)) then
      rel_tol_ = rel_tol
    end if

    if (present(abs_tol)) then
      abs_tol_ = abs_tol
    end if

    isclose_real = abs(a - b) <= max(rel_tol_ * max(abs(a), abs(b)), abs_tol_)
  end function isclose_real

  pure logical function isclose_complex(a, b, rel_tol, abs_tol) 
    !* Evaluates the closeness between two complex variables.
    ! 
    ! `.true.` if `a` and `b` are close to each other according to a tolerance, and
    ! `.false.` otherwise. Tolerance is defined by `rel_tol` and `abs_tol`, where
    !* `abs_tol` is used for comparing values close no zero.

    complex(wp), intent(in) :: a, b
    real(wp), intent(in), optional :: rel_tol  !! Relative tolerance, default=`tol_wp`
    real(wp), intent(in), optional :: abs_tol  !! Absolute tolerace, default=`tol_wp`
    real(wp) :: are, aim, bre, bim
    real(wp) :: rel_tol_, abs_tol_

    rel_tol_ = tol_wp
    abs_tol_ = tol_wp

    if (present(rel_tol)) then
      rel_tol_ = rel_tol
    end if

    if (present(abs_tol)) then
      abs_tol_ = abs_tol
    end if

    are = a%re
    aim = a%im
    bre = b%re
    bim = b%im

    isclose_complex = isclose_real(are, bre, rel_tol_, abs_tol_) .and. &
                      isclose_real(aim, bim, rel_tol_, abs_tol_)
  end function isclose_complex

end module csf_numerror
