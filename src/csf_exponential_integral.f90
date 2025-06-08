!  ┏┓┏┓┏┓  Licensed under the MIT License
!  ┃ ┗┓┣   Copyright (c) 2025 Rodrigo Castro 
!  ┗┛┗┛┻   https://github.com/rodpcastro/colspecf

module csf_exponential_integral
!* # Exponential integral
! Exponential integrals.
!
! Procedures:
!
! - `ei`: Exponential integral \(\mathrm{Ei}(x)\)
! - `e1`: Exponential integral \(\mathrm{E}_1(x)\) or \(\mathrm{E}_1(z)\)
!
! Untested procedures:  
!
! - `enz`: Exponential integral \(\mathrm{E}_n(z)\), untested for \(n \neq 1\)
!
! ## References
! 1. Kathleen A. Paciorek. 1970. Algorithm 385: Exponential integral Ei(x). Commun.
!    ACM 13, 7 (July 1970), 446–447. <https://doi.org/10.1145/362686.362696>
! 2. Donald E. Amos. 1990. Algorithms 683: a portable FORTRAN subroutine for
!    exponential integrals of a complex argument. ACM Trans. Math. Softw. 16,
!*   2 (June 1990), 178–182. <https://doi.org/10.1145/78928.78934>

  use, intrinsic :: iso_fortran_env, only: stderr => error_unit
  use csf_kinds, only: i4, wp
  use csf_constants, only: pi, ninf, pinf
  use csf_numerror, only: eps_wp
  use calgo_385, only: dei
  use calgo_683, only: cexint

  implicit none
  private
  public :: ei, e1, enz

  interface e1
    !! Exponential integral \(\mathrm{E}_1(x)\) or \(\mathrm{E}_1(z)\).
    module procedure e1x, e1z
  end interface e1

contains

  real(wp) function ei(x)
    !! Exponential integral \(\mathrm{Ei}(x)\).
    !
    !! \(\lbrace x \in \mathbb{R} \mid x \neq 0 \rbrace\)

    real(wp), intent(in) :: x

    if (x == 0.0_wp) then
      ei = ninf()
    else if (x < -738.0_wp) then
      ei = 0.0_wp
    else if (x <= 709.0_wp) then
      ei = dei(x)
    else
      ei = pinf()
    end if
  end function ei

  real(wp) function e1x(x)
    !! Exponential integral \(\mathrm{E}_1(x)\).
    !
    !! \(\lbrace x \in \mathbb{R} \mid x \neq 0 \rbrace\)

    real(wp), intent(in) :: x
    complex(wp) :: z, e1z

    if (x == 0.0_wp) then
      e1x = pinf()
    else if (x <= 738.0_wp) then
      z = cmplx(x, 0.0_wp, kind=wp)
      e1z = enz(1, z)
      e1x = e1z%re
    else
      e1x = 0.0_wp
    end if
  end function e1x

  complex(wp) function e1z(z)
    !! Exponential integral \(\mathrm{E}_1(z)\).
    !
    !! \(z \in \mathbb{C} \setminus \left( \lbrace z \in \mathbb{C} \mid \Re(z) \lt 0,
    !! \thinspace 0 \lt |\Im(z)| \lt 10^{-6} \rbrace \cup \lbrace 0 \rbrace \right)\)

    complex(wp), intent(in) :: z
    real(wp) :: zabs

    zabs = abs(z)

    if (zabs == 0.0_wp) then
      e1z = cmplx(pinf(), -pi, kind=wp)
      return
    else if (z%re < -700.0_wp .and. z%im == 0.0_wp) then
      e1z = cmplx(ninf(), -pi, kind=wp)
      return
    else if (z%re > 738.0_wp .or. (z%re >= 0 .and. abs(z%im) > 1.0e15_wp)) then
      ! For Re(z) ≥ 0 and Im(z) > 1e15, CALGO 683 throws IERR = 5.
      e1z = (0.0_wp, 0.0_wp)
      return
    else
      e1z = enz(1, z, .false.)
    end if

    if (z%re < 0.0_wp .and. z%im == 0.0_wp) then
      e1z = cmplx(e1z%re, -pi, kind=wp)
    end if
  end function e1z

  complex(wp) function enz(n, z, show_warning)
    !* Exponential integral \(\mathrm{E}_n(z)\).
    ! 
    ! \(n \geq 1,\thinspace z \in \mathbb{C},\thinspace -\pi \lt \arg(z) \leq \pi \)
    !
    ! If `show_warning = .true.`, a warning message is displayed when `cexint` returns
    !* `IERR = 2` or `IERR = 4`.

    integer(i4), intent(in) :: n
    complex(wp), intent(in) :: z
    logical, intent(in), optional :: show_warning  !! Default=`.true.`
    real(wp), parameter :: tol = eps_wp
    integer(i4), parameter :: m = 1
    complex(wp) :: cy(m)
    integer(i4) :: ierr

    logical :: show_warning_
    
    ! write(stderr, *) 'show_warning = ', show_warning

    if (present(show_warning)) then
      show_warning_ = show_warning
    else
      show_warning_ = .true.
    end if

    ! write(stderr, *) 'show_warning = ', show_warning_

    ! Computing En(z) with no scaling (KODE = 1).
    call cexint(z, n, 1, tol, m, cy, ierr)
    enz = cy(1)

    select case (ierr)
      case (1)
        error stop 'CALGO 683 IERR = 1: An input error. No computation.'
      case (2)
        ! Computing E1(z) with scaling (KODE = 2).
        call cexint(z, n, 2, tol, m, cy, ierr)
        enz = exp(-z) * cy(1)

        if (show_warning_) then
          write(stderr, '(l, a)') show_warning_, 'CALGO 683 IERR = 2: Underflow. ' // &
                               'En(z) = (0.0, 0.0). Real(z) > 0.0 ' // &
                               'too large on KODE = 1. Recomputed with KODE = 2.'
        end if
      case (3)
        error stop 'CALGO 683 IERR = 3: Overflow. No computation. ' // &
                   'Real(z) < 0.0 too small on KODE = 1.'
      case (4)
        if (show_warning_ .eqv. .true.) then
          write(stderr, '(a)') 'CALGO 683 IERR = 4: |z| or n large. '  // &
                               'Computation done but losses of significance by ' // &
                               'argument reduction may exceed half precision.'
        end if
      case (5)
        error stop 'CALGO 683 IERR = 5: |z| or n large. No computation. ' // &
                   'All loss of significance by argument reduction has occurred.'
      case (6)
        error stop 'CALGO 683 IERR = 6: Convergence error. No computation. ' // &
                   'Algorithm termination condition not met.'
      case (7)
        error stop 'CALGO 683 IERR = 7: Discrimination error. No computation. ' // &
                   'This condition should never occur.'
    end select
  end function enz

  ! complex(wp) function enz(n, z)
  !   !! Exponential integral \(\mathrm{E}_n(z)\).
  !   ! 
  !   !! \(n \geq 1,\thinspace z \in \mathbb{C},\thinspace -\pi \lt \arg(z) \leq \pi \)
  !
  !   integer(i4), intent(in) :: n
  !   complex(wp), intent(in) :: z
  !   real(wp), parameter :: tol = eps_wp
  !   integer(i4), parameter :: m = 1
  !   complex(wp) :: cy(m)
  !   integer(i4) :: ierr
  !
  !   ! Computing En(z) with no scaling (KODE = 1).
  !   call cexint(z, n, 1, tol, m, cy, ierr)
  !   enz = cy(1)
  !
  !   select case (ierr)
  !     case (1)
  !       error stop 'CALGO 683 IERR = 1: An input error. No computation.'
  !     case (2)
  !       if (n == 1) then
  !         ! Original behavior has been overwritten for E1(z).
  !         ! Computing E1(z) with scaling (KODE = 2).
  !         call cexint(z, n, 2, tol, m, cy, ierr)
  !         enz = exp(-z) * cy(1)
  !       else
  !         write(stderr, '(a)') 'CALGO 683 IERR = 2: Underflow. ' // &
  !                              'En(z) = (0.0, 0.0). Real(z) > 0.0 ' // &
  !                              'too large on KODE = 1.'
  !       end if
  !     case (3)
  !       error stop 'CALGO 683 IERR = 3: Overflow. No computation. ' // &
  !                  'Real(z) < 0.0 too small on KODE = 1.'
  !     case (4)
  !       write(stderr, '(a)') 'CALGO 683 IERR = 4: |z| or n large. '  // &
  !                            'Computation done but losses of significance by ' // &
  !                            'argument reduction may exceed half precision.'
  !     case (5)
  !       if (n == 1) then
  !         ! Original behavior has been overwritten for E1(z).
  !         write(stderr, '(a)') 'CALGO 683 IERR = 5: |z| large. All loss of ' // &
  !                              'significance by argument reduction has ' // &
  !                              'occurred. Returning E1(z) = (0.0, 0.0).'
  !         enz = (0.0_wp, 0.0_wp)
  !       else
  !         error stop 'CALGO 683 IERR = 5: |z| or n large. No computation. ' // &
  !                    'All loss of significance by argument reduction has occurred.'
  !       end if
  !     case (6)
  !       error stop 'CALGO 683 IERR = 6: Convergence error. No computation. ' // &
  !                  'Algorithm termination condition not met.'
  !     case (7)
  !       error stop 'CALGO 683 IERR = 7: Discrimination error. No computation. ' // &
  !                  'This condition should never occur.'
  !   end select
  ! end function enz

end module csf_exponential_integral
