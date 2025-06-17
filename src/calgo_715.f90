! Licensed under the ACM Software License Agreement
! Copyright © 1970–2012 Association for Computing Machinery (ACM).
! https://www.acm.org/publications/policies/software-copyright-notice

module calgo_715
!* # CALGO 715
! Algorithm 715. 
!
! Procedures:
!
! - `caljy0`: Bessel functions \(J_0(x)\) and \(Y_0(x)\)
! - `caljy1`: Bessel functions \(J_1(x)\) and \(Y_1(x)\)
!
! ## Author
! W. J. Cody
!
! ## History
! - 1992-03-15 - W. J. Cody
!     - Original code.
! - 2000-03-30 - W. J. Cody
!     - F77 code distributed by ACM: <https://calgo.acm.org/715.zip>
! - 2003-01-14 - Alan Miller
!     - F90 code adaptation by Alan Miller:
!       <https://jblevins.org/mirror/amiller/toms715.zip>
! - 2025-06-17 - Rodrigo Castro (GitHub: rodpcastro)
!     - Retained only subroutines `caljy0` and `caljy1`; additional subroutines will
!       be included as required.
!     - Replaced `dp` (double precision) by `wp` (working precision)
!     - Replaced array constructor `(/.../)` by the less verbose `[...]`
!
! ## References
! 1. W. J. Cody. 1993. Algorithm 715: SPECFUN–a portable FORTRAN package of special
!    function routines and test drivers. ACM Trans. Math. Softw. 19, 1 (March 1993),
!*   22–30. <https://doi.org/10.1145/151271.151273>

  use csf_kinds, only: wp

  implicit none
  private
  public :: caljy0, caljy1

contains

  SUBROUTINE caljy0(arg, result, jint)
    !* CALGO 715 Bessel functions \(J_0(x)\) and \(Y_0(x)\).
    !
    ! To obtain:
    ! 
    ! - `y` = \(J_0(x)\), call `caljy0(x, y, 0)`
    ! - `y` = \(Y_0(x)\), call `caljy0(x, y, 1)`
    !*

    !-------------------------------------------------------------------
    ! This packet computes zero-order Bessel functions of the first and
    ! second kind (J0 and Y0), for real arguments X, where 0 < X <= XMAX
    ! for Y0, and |X| <= XMAX for J0.  It contains two function-type
    ! subprograms, BESJ0  and  BESY0, and one subroutine-type
    ! subprogram, CALJY0.  The calling statements for the primary
    ! entries are:
    !
    ! Y = BESJ0(X)
    !
    ! and
    !
    ! Y = BESY0(X),
    !
    ! where the entry points correspond to the functions J0(X) and Y0(X),
    ! respectively. The routine CALJY0 is intended for internal packet
    ! use only, all computations within the packet being concentrated in
    ! this one routine. The function subprograms invoke CALJY0 with
    ! the statement
    !
    ! CALL CALJY0(ARG, RESULT, JINT),
    !
    ! where the parameter usage is as follows:
    !
    ! Function     Parameters for CALJY0
    ! call         ARG                    RESULT    JINT
    !
    ! BESJ0(ARG)   |ARG| .LE. XMAX        J0(ARG)   0
    ! BESY0(ARG)   0 .LT. ARG .LE. XMAX   Y0(ARG)   1
    !
    ! The main computation uses unpublished minimax rational
    ! approximations for X .LE. 8.0, and an approximation from the
    ! book  Computer Approximations  by Hart, et. al., Wiley and Sons,
    ! New York, 1968, for arguments larger than 8.0  Part of this
    ! transportable packet is patterned after the machine-dependent
    ! FUNPACK program BESJ0(X), but cannot match that version for
    ! efficiency or accuracy. This version uses rational functions
    ! that are theoretically accurate to at least 18 significant decimal
    ! digits for X <= 8, and at least 18 decimal places for X > 8. The
    ! accuracy achieved depends on the arithmetic system, the compiler,
    ! the intrinsic functions, and proper selection of the machine-
    ! dependent constants.
    !
    !-------------------------------------------------------------------
    !
    ! The following machine-dependent constants must be declared in
    ! DATA statements. IEEE values are provided as a default.
    !
    ! XINF   = largest positive machine number
    ! XMAX   = largest acceptable argument. The functions AINT, SIN
    !          and COS must perform properly for ABS(X) .LE. XMAX.
    !          We recommend that XMAX be a small integer multiple of
    !          sqrt(1/eps), where eps is the smallest positive number
    !          such that  1+eps > 1.
    ! XSMALL = positive argument such that 1.0-(X/2)**2 = 1.0
    !          to machine precision for all ABS(X) .LE. XSMALL.
    !          We recommend that  XSMALL < sqrt(eps)/beta, where beta
    !          is the floating-point radix (usually 2 or 16).
    !
    ! Approximate values for some important machines are
    !
    !                       eps       XMAX      XSMALL    XINF
    !
    ! CDC 7600      (S.P.)  7.11E-15  1.34E+08  2.98E-08  1.26E+322
    ! CRAY-1        (S.P.)  7.11E-15  1.34E+08  2.98E-08  5.45E+2465
    ! IBM PC (8087) (S.P.)  5.96E-08  8.19E+03  1.22E-04  3.40E+38
    ! IBM PC (8087) (D.P.)  1.11D-16  2.68D+08  3.72D-09  1.79D+308
    ! IBM 195       (D.P.)  2.22D-16  6.87D+09  9.09D-13  7.23D+75
    ! UNIVAC 1108   (D.P.)  1.73D-18  4.30D+09  2.33D-10  8.98D+307
    ! VAX 11/780    (D.P.)  1.39D-17  1.07D+09  9.31D-10  1.70D+38
    !
    !-------------------------------------------------------------------
    !
    ! Error Returns
    !
    ! The program returns the value zero for X .GT. XMAX, and returns
    ! -XINF when BESLY0 is called with a negative or zero argument.
    !-------------------------------------------------------------------

    real(wp), INTENT(IN)  :: arg
    real(wp), INTENT(OUT) :: result
    INTEGER, INTENT(IN)   :: jint

    INTEGER  :: i
    real(wp) :: ax, down, prod, resj, r0, r1, up, w, wsq, xden, xnum, xy, z, zsq
    !-------------------------------------------------------------------
    ! Mathematical constants
    ! CONS = ln(.5) + Euler's gamma
    !-------------------------------------------------------------------
    real(wp), PARAMETER :: zero = 0.0_wp, one = 1.0_wp, three = 3.0_wp, &
                           four = 4.0_wp, eight = 8.0_wp, five5 = 5.5_wp, &
                           sixty4 = 64.0_wp, oneov8 = 0.125_wp, &
                           p17 = 0.1716_wp, two56 = 256.0_wp, &
                           cons = -1.1593151565841244881e-1_wp, &
                           pi2 = 6.3661977236758134308e-1_wp, &
                           twopi = 6.2831853071795864769_wp, &
                           twopi1 = 6.28125_wp, twopi2 = 1.9353071795864769253e-3_wp
    !-------------------------------------------------------------------
    ! Machine-dependent constants
    !-------------------------------------------------------------------
    real(wp), PARAMETER :: XMAX = 2.68e+8_wp, XSMALL = 3.72e-9_wp, XINF = 1.79e+308_wp
    !-------------------------------------------------------------------
    ! Zeroes of Bessel functions
    !-------------------------------------------------------------------
    real(wp), PARAMETER :: XJ0 = 2.4048255576957727686_wp, &
                           XJ1 = 5.5200781102863106496_wp, &
                           XY0 = 8.9357696627916752158e-1_wp, &
                           XY1 = 3.9576784193148578684_wp, &
                           XY2 = 7.0860510603017726976_wp, &
                           XJ01 =  616.0_wp, XJ02 = -1.4244423042272313784e-3_wp, &
                           XJ11 = 1413.0_wp, XJ12 =  5.4686028631064959660e-4_wp, &
                           XY01 =  228.0_wp, XY02 =  2.9519662791675215849e-3_wp, &
                           XY11 = 1013.0_wp, XY12 =  6.4716931485786837568e-4_wp, &
                           XY21 = 1814.0_wp, XY22 =  1.1356030177269762362e-4_wp
    !-------------------------------------------------------------------
    ! Coefficients for rational approximation to ln(x/a)
    !--------------------------------------------------------------------
    real(wp), PARAMETER :: PLG(4) = [ &
      -2.4562334077563243311e+1_wp, 2.3642701335621505212e+2_wp, &
      -5.4989956895857911039e+2_wp, 3.5687548468071500413e+2_wp &
    ]
    real(wp), PARAMETER :: QLG(4) = [ &
      -3.5553900764052419184e+1_wp, 1.9400230218539473193e+2_wp, &
      -3.3442903192607538956e+2_wp, 1.7843774234035750207e+2_wp &
    ]
    !-------------------------------------------------------------------
    ! Coefficients for rational approximation of
    ! J0(X) / (X**2 - XJ0**2), XSMALL < |X| <= 4.0
    !--------------------------------------------------------------------
    real(wp), PARAMETER :: PJ0(7) = [ &
      6.6302997904833794242e+6_wp, -6.2140700423540120665e+8_wp, &
      2.7282507878605942706e+10_wp, -4.1298668500990866786e+11_wp, &
      -1.2117036164593528341e-1_wp, 1.0344222815443188943e+2_wp, &
      -3.6629814655107086448e+4_wp &
    ]
    real(wp), PARAMETER :: QJ0(5) = [ &
      4.5612696224219938200e+5_wp, 1.3985097372263433271e+8_wp, &
      2.6328198300859648632e+10_wp, 2.3883787996332290397e+12_wp, &
      9.3614022392337710626e+2_wp &
    ]
    !-------------------------------------------------------------------
    ! Coefficients for rational approximation of
    ! J0(X) / (X**2 - XJ1**2), 4.0 < |X| <= 8.0
    !-------------------------------------------------------------------
    real(wp), PARAMETER :: PJ1(8) = [ &
      4.4176707025325087628e+3_wp, 1.1725046279757103576e+4_wp, &
      1.0341910641583726701e+4_wp, -7.2879702464464618998e+3_wp, &
      -1.2254078161378989535e+4_wp, -1.8319397969392084011e+3_wp, &
      4.8591703355916499363e+1_wp, 7.4321196680624245801e+2_wp &
    ]
    real(wp), PARAMETER :: QJ1(7) = [ &
      3.3307310774649071172e+2_wp, -2.9458766545509337327e+3_wp, &
      1.8680990008359188352e+4_wp, -8.4055062591169562211e+4_wp, &
      2.4599102262586308984e+5_wp, -3.5783478026152301072e+5_wp, &
      -2.5258076240801555057e+1_wp &
    ]
    !-------------------------------------------------------------------
    ! Coefficients for rational approximation of
    ! (Y0(X) - 2 LN(X/XY0) J0(X)) / (X**2 - XY0**2), XSMALL < |X| <= 3.0
    !--------------------------------------------------------------------
    real(wp), PARAMETER :: PY0(6) = [ &
      1.0102532948020907590e+4_wp, -2.1287548474401797963e+6_wp, &
      2.0422274357376619816e+8_wp, -8.3716255451260504098e+9_wp, &
      1.0723538782003176831e+11_wp, -1.8402381979244993524e+1_wp &
    ]
    real(wp), PARAMETER :: QY0(5) = [ &
      6.6475986689240190091e+2_wp, 2.3889393209447253406e+5_wp, &
      5.5662956624278251596e+7_wp, 8.1617187777290363573e+9_wp, &
      5.8873865738997033405e+11_wp &
    ]
    !-------------------------------------------------------------------
    ! Coefficients for rational approximation of
    ! (Y0(X) - 2 LN(X/XY1) J0(X)) / (X**2 - XY1**2), 3.0 < |X| <= 5.5
    !--------------------------------------------------------------------
    real(wp), PARAMETER :: PY1(7) = [ &
      -1.4566865832663635920e+4_wp, 4.6905288611678631510e+6_wp, &
      -6.9590439394619619534e+8_wp, 4.3600098638603061642e+10_wp, &
      -5.5107435206722644429e+11_wp, -2.2213976967566192242e+13_wp, &
      1.7427031242901594547e+1_wp &
    ]
    real(wp), PARAMETER :: QY1(6) = [ &
      8.3030857612070288823e+2_wp, 4.0669982352539552018e+5_wp, &
      1.3960202770986831075e+8_wp, 3.4015103849971240096e+10_wp, &
      5.4266824419412347550e+12_wp, 4.3386146580707264428e+14_wp &
    ]
    !-------------------------------------------------------------------
    ! Coefficients for rational approximation of
    ! (Y0(X) - 2 LN(X/XY2) J0(X)) / (X**2 - XY2**2), 5.5 < |X| <= 8.0
    !--------------------------------------------------------------------
    real(wp), PARAMETER :: PY2(8) = [ &
      2.1363534169313901632e+4_wp, -1.0085539923498211426e+7_wp, &
      2.1958827170518100757e+9_wp, -1.9363051266772083678e+11_wp, &
      -1.2829912364088687306e+11_wp, 6.7016641869173237784e+14_wp, &
      -8.0728726905150210443e+15_wp, -1.7439661319197499338e+1_wp &
    ]
    real(wp), PARAMETER :: QY2(7) = [ &
      8.7903362168128450017e+2_wp, 5.3924739209768057030e+5_wp, &
      2.4727219475672302327e+8_wp, 8.6926121104209825246e+10_wp, &
      2.2598377924042897629e+13_wp, 3.9272425569640309819e+15_wp, &
      3.4563724628846457519e+17_wp &
    ]
    !-------------------------------------------------------------------
    ! Coefficients for Hart's approximation, |X| > 8.0
    !-------------------------------------------------------------------
    real(wp), PARAMETER :: P0(6) = [ &
      3.4806486443249270347e+3_wp, 2.1170523380864944322e+4_wp, &
      4.1345386639580765797e+4_wp, 2.2779090197304684302e+4_wp, &
      8.8961548424210455236e-1_wp, 1.5376201909008354296e+2_wp &
    ]
    real(wp), PARAMETER :: Q0(5) = [ &
      3.5028735138235608207e+3_wp, 2.1215350561880115730e+4_wp, &
      4.1370412495510416640e+4_wp, 2.2779090197304684318e+4_wp, &
      1.5711159858080893649e+2_wp &
    ]
    real(wp), PARAMETER :: P1(6) = [ &
      -2.2300261666214198472e+1_wp, -1.1183429920482737611e+2_wp, &
      -1.8591953644342993800e+2_wp, -8.9226600200800094098e+1_wp, &
      -8.8033303048680751817e-3_wp, -1.2441026745835638459_wp &
    ]
    real(wp), PARAMETER :: Q1(5) = [ &
      1.4887231232283756582e+3_wp, 7.2642780169211018836e+3_wp, &
      1.1951131543434613647e+4_wp, 5.7105024128512061905e+3_wp, &
      9.0593769594993125859e+1_wp &
    ]
    !-------------------------------------------------------------------
    ! Check for error conditions
    !-------------------------------------------------------------------
    ax = ABS(arg)
    IF (jint == 1 .AND. arg <= zero) THEN
      result = -xinf
      GO TO 80
    ELSE IF (ax > xmax) THEN
      result = zero
      GO TO 80
    END IF
    IF (ax <= eight) THEN
      IF (ax <= xsmall) THEN
        IF (jint == 0) THEN
          result = one
        ELSE
          result = pi2 * (LOG(ax)+cons)
        END IF
        GO TO 80
      END IF
    !-------------------------------------------------------------------
    ! Calculate J0 for appropriate interval, preserving
    ! accuracy near the zero of J0
    !-------------------------------------------------------------------
      zsq = ax * ax
      IF (ax <= four) THEN
        xnum = (pj0(5)*zsq + pj0(6)) * zsq + pj0(7)
        xden = zsq + qj0(5)
        DO  i = 1, 4
          xnum = xnum * zsq + pj0(i)
          xden = xden * zsq + qj0(i)
        END DO
        prod = ((ax-xj01/two56)-xj02) * (ax+xj0)
      ELSE
        wsq = one - zsq / sixty4
        xnum = pj1(7) * wsq + pj1(8)
        xden = wsq + qj1(7)
        DO  i = 1, 6
          xnum = xnum * wsq + pj1(i)
          xden = xden * wsq + qj1(i)
        END DO
        prod = (ax+xj1) * ((ax-xj11/two56)-xj12)
      END IF
      result = prod * xnum / xden
      IF (jint == 0) GO TO 80
    !-------------------------------------------------------------------
    ! Calculate Y0.  First find  RESJ = pi/2 ln(x/xn) J0(x),
    ! where xn is a zero of Y0
    !-------------------------------------------------------------------
      IF (ax <= three) THEN
        up = (ax-xy01/two56) - xy02
        xy = xy0
      ELSE IF (ax <= five5) THEN
        up = (ax-xy11/two56) - xy12
        xy = xy1
      ELSE
        up = (ax-xy21/two56) - xy22
        xy = xy2
      END IF
      down = ax + xy
      IF (ABS(up) < p17*down) THEN
        w = up / down
        wsq = w * w
        xnum = plg(1)
        xden = wsq + qlg(1)
        DO  i = 2, 4
          xnum = xnum * wsq + plg(i)
          xden = xden * wsq + qlg(i)
        END DO
        resj = pi2 * result * w * xnum / xden
      ELSE
        resj = pi2 * result * LOG(ax/xy)
      END IF
    !-------------------------------------------------------------------
    ! Now calculate Y0 for appropriate interval, preserving
    ! accuracy near the zero of Y0
    !-------------------------------------------------------------------
      IF (ax <= three) THEN
        xnum = py0(6) * zsq + py0(1)
        xden = zsq + qy0(1)
        DO  i = 2, 5
          xnum = xnum * zsq + py0(i)
          xden = xden * zsq + qy0(i)
        END DO
      ELSE IF (ax <= five5) THEN
        xnum = py1(7) * zsq + py1(1)
        xden = zsq + qy1(1)
        DO  i = 2, 6
          xnum = xnum * zsq + py1(i)
          xden = xden * zsq + qy1(i)
        END DO
      ELSE
        xnum = py2(8) * zsq + py2(1)
        xden = zsq + qy2(1)
        DO  i = 2, 7
          xnum = xnum * zsq + py2(i)
          xden = xden * zsq + qy2(i)
        END DO
      END IF
      result = resj + up * down * xnum / xden
    ELSE
    !-------------------------------------------------------------------
    ! Calculate J0 or Y0 for |ARG| > 8.0
    !-------------------------------------------------------------------
      z = eight / ax
      w = ax / twopi
      w = AINT(w) + oneov8
      w = (ax-w*twopi1) - w * twopi2
      zsq = z * z
      xnum = p0(5) * zsq + p0(6)
      xden = zsq + q0(5)
      up = p1(5) * zsq + p1(6)
      down = zsq + q1(5)
      DO  i = 1, 4
        xnum = xnum * zsq + p0(i)
        xden = xden * zsq + q0(i)
        up = up * zsq + p1(i)
        down = down * zsq + q1(i)
      END DO
      r0 = xnum / xden
      r1 = up / down
      IF (jint == 0) THEN
        result = SQRT(pi2/ax) * (r0*COS(w) - z*r1*SIN(w))
      ELSE
        result = SQRT(pi2/ax) * (r0*SIN(w) + z*r1*COS(w))
      END IF
    END IF
  80 RETURN
  END SUBROUTINE caljy0

  SUBROUTINE caljy1(arg, result, jint)
    !* CALGO 715 Bessel functions \(J_1(x)\) and \(Y_1(x)\).
    !
    ! To obtain:
    ! 
    ! - `y` = \(J_1(x)\), call `caljy1(x, y, 0)`
    ! - `y` = \(Y_1(x)\), call `caljy1(x, y, 1)`
    !*

    !-------------------------------------------------------------------
    ! This packet computes first-order Bessel functions of the first and
    ! second kind (J1 and Y1), for real arguments X, where 0 < X <= XMAX
    ! for Y1, and |X| <= XMAX for J1. It contains two function-type
    ! subprograms, BESJ1 and BESY1, and one subroutine-type
    ! subprogram, CALJY1. The calling statements for the primary
    ! entries are:
    !
    ! Y = BESJ1(X)
    !
    !   and
    !
    ! Y = BESY1(X),
    !
    ! where the entry points correspond to the functions J1(X) and Y1(X),
    ! respectively.  The routine  CALJY1  is intended for internal packet
    ! use only, all computations within the packet being concentrated in
    ! this one routine.  The function subprograms invoke  CALJY1  with
    ! the statement
    !
    ! CALL CALJY1(ARG, RESULT, JINT),
    !
    ! where the parameter usage is as follows:
    !
    ! Function     Parameters for CALJY1
    ! call         ARG                    RESULT    JINT
    !
    ! BESJ1(ARG)   |ARG| .LE. XMAX        J1(ARG)   0
    ! BESY1(ARG)   0 .LT. ARG .LE. XMAX   Y1(ARG)   1
    !
    ! The main computation uses unpublished minimax rational
    ! approximations for X .LE. 8.0, and an approximation from the
    ! book Computer Approximations  by Hart, et. al., Wiley and Sons,
    ! New York, 1968, for arguments larger than 8.0 Part of this
    ! transportable packet is patterned after the machine-dependent
    ! FUNPACK program BESJ1(X), but cannot match that version for
    ! efficiency or accuracy. This version uses rational functions
    ! that are theoretically accurate to at least 18 significant decimal
    ! digits for X <= 8, and at least 18 decimal places for X > 8. The
    ! accuracy achieved depends on the arithmetic system, the compiler,
    ! the intrinsic functions, and proper selection of the machine-
    ! dependent constants.
    !
    !-------------------------------------------------------------------
    !
    ! The following machine-dependent constants must be declared in
    ! DATA statements.  IEEE values are provided as a default.
    !
    ! XINF   = largest positive machine number
    ! XMAX   = largest acceptable argument. The functions AINT, SIN
    !          and COS must perform properly for ABS(X) .LE. XMAX.
    !          We recommend that XMAX be a small integer multiple of
    !          sqrt(1/eps), where eps is the smallest positive number
    !          such that  1+eps > 1.
    ! XSMALL = positive argument such that 1.0-(1/2)(X/2)**2 = 1.0
    !          to machine precision for all ABS(X) .LE. XSMALL.
    !          We recommend that  XSMALL < sqrt(eps)/beta, where beta
    !          is the floating-point radix (usually 2 or 16).
    !
    ! Approximate values for some important machines are
    !
    !                       eps       XMAX      XSMALL    XINF
    !
    ! CDC 7600      (S.P.)  7.11E-15  1.34E+08  2.98E-08  1.26E+322
    ! CRAY-1        (S.P.)  7.11E-15  1.34E+08  2.98E-08  5.45E+2465
    ! IBM PC (8087) (S.P.)  5.96E-08  8.19E+03  1.22E-04  3.40E+38
    ! IBM PC (8087) (D.P.)  1.11D-16  2.68D+08  3.72D-09  1.79D+308
    ! IBM 195       (D.P.)  2.22D-16  6.87D+09  9.09D-13  7.23D+75
    ! UNIVAC 1108   (D.P.)  1.73D-18  4.30D+09  2.33D-10  8.98D+307
    ! VAX 11/780    (D.P.)  1.39D-17  1.07D+09  9.31D-10  1.70D+38
    !
    !--------------------------------------------------------------------
    !
    ! Error Returns
    !
    ! The program returns the value zero for  X .GT. XMAX, and returns
    ! -XINF when BESLY1 is called with a negative or zero argument.
    !--------------------------------------------------------------------

    real(wp), INTENT(IN)  :: arg
    real(wp), INTENT(OUT) :: result
    INTEGER, INTENT(IN)   :: jint

    INTEGER  :: i
    real(wp) :: ax, down, prod, resj, r0, r1, up, w, wsq, xden, xnum, xy, z, zsq
    !-------------------------------------------------------------------
    ! Mathematical constants
    !-------------------------------------------------------------------
    real(wp), PARAMETER :: eight = 8.0_wp, four = 4.0_wp, half = 0.5_wp, &
                           throv8 = 0.375_wp, p17 = 0.1716_wp, &
                           pi2 = 6.3661977236758134308e-1_wp, zero = 0.0_wp, &
                           twopi = 6.2831853071795864769_wp, twopi1 = 6.28125_wp, &
                           twopi2 = 1.9353071795864769253e-3_wp, two56 = 256.0_wp, &
                           rtpi2 = 7.9788456080286535588e-1_wp
    !-------------------------------------------------------------------
    ! Machine-dependent constants
    !-------------------------------------------------------------------
    real(wp), PARAMETER :: XMAX = 2.68e+8_wp, XSMALL = 3.72e-9_wp, XINF = 1.79e+308_wp
    !-------------------------------------------------------------------
    ! Zeroes of Bessel functions
    !-------------------------------------------------------------------
    real(wp), PARAMETER :: & 
      XJ0 = 3.8317059702075123156_wp, XJ1 = 7.0155866698156187535_wp, &
      XY0 = 2.1971413260310170351_wp, XY1 = 5.4296810407941351328_wp, &
      XJ01 = 981.0_wp, XJ02 = -3.2527979248768438556e-4_wp, &
      XJ11 = 1796.0_wp, XJ12 = -3.8330184381246462950e-5_wp, &
      XY01 = 562.0_wp, XY02 = 1.8288260310170351490e-3_wp, &
      XY11 = 1390.0_wp, XY12 = -6.4592058648672279948e-6_wp
    !-------------------------------------------------------------------
    ! Coefficients for rational approximation to ln(x/a)
    !--------------------------------------------------------------------
    real(wp), PARAMETER :: PLG(4) = [ &
      -2.4562334077563243311e+1_wp, 2.3642701335621505212e+2_wp, &
      -5.4989956895857911039e+2_wp, 3.5687548468071500413e+2_wp &
    ]
    real(wp), PARAMETER :: QLG(4) = [ &
      -3.5553900764052419184e+1_wp, 1.9400230218539473193e+2_wp, &
      -3.3442903192607538956e+2_wp, 1.7843774234035750207e+2_wp &
    ]
    !-------------------------------------------------------------------
    ! Coefficients for rational approximation of
    ! J1(X) / (X * (X**2 - XJ0**2)), XSMALL < |X| <= 4.0
    !--------------------------------------------------------------------
    real(wp), PARAMETER :: PJ0(7) = [ &
      9.8062904098958257677e+5_wp, -1.1548696764841276794e+8_wp, &
      6.6781041261492395835e+9_wp, -1.4258509801366645672e+11_wp, &
      -4.4615792982775076130e+3_wp, 1.0650724020080236441e+1_wp, &
      -1.0767857011487300348e-2_wp &
    ]
    real(wp), PARAMETER :: QJ0(5) = [ &
      5.9117614494174794095e+5_wp, 2.0228375140097033958e+8_wp, &
      4.2091902282580133541e+10_wp, 4.1868604460820175290e+12_wp, &
      1.0742272239517380498e+3_wp &
    ]
    !-------------------------------------------------------------------
    ! Coefficients for rational approximation of
    ! J1(X) / (X * (X**2 - XJ1**2)), 4.0 < |X| <= 8.0
    !-------------------------------------------------------------------
    real(wp), PARAMETER :: PJ1(8) = [ &
      4.6179191852758252280_wp, -7.1329006872560947377e+3_wp, &
      4.5039658105749078904e+6_wp, -1.4437717718363239107e+9_wp, &
      2.3569285397217157313e+11_wp, -1.6324168293282543629e+13_wp, &
      1.1357022719979468624e+14_wp, 1.0051899717115285432e+15_wp &
    ]
    real(wp), PARAMETER :: QJ1(7) = [ &
      1.1267125065029138050e+6_wp, 6.4872502899596389593e+8_wp, &
      2.7622777286244082666e+11_wp, 8.4899346165481429307e+13_wp, &
      1.7128800897135812012e+16_wp, 1.7253905888447681194e+18_wp, &
      1.3886978985861357615e+3_wp &
    ]
    !-------------------------------------------------------------------
    ! Coefficients for rational approximation of
    ! (Y1(X) - 2 LN(X/XY0) J1(X)) / (X**2 - XY0**2), XSMALL < |X| <= 4.0
    !--------------------------------------------------------------------
    real(wp), PARAMETER :: PY0(7) = [ &
      2.2157953222280260820e+5_wp, -5.9157479997408395984e+7_wp, &
      7.2144548214502560419e+9_wp, -3.7595974497819597599e+11_wp, &
      5.4708611716525426053e+12_wp, 4.0535726612579544093e+13_wp, &
      -3.1714424660046133456e+2_wp &
    ]
    real(wp), PARAMETER :: QY0(6) = [ &
      8.2079908168393867438e+2_wp, 3.8136470753052572164e+5_wp, &
      1.2250435122182963220e+8_wp, 2.7800352738690585613e+10_wp, &
      4.1272286200406461981e+12_wp, 3.0737873921079286084e+14_wp &
    ]
    !--------------------------------------------------------------------
    ! Coefficients for rational approximation of
    ! (Y1(X) - 2 LN(X/XY1) J1(X)) / (X**2 - XY1**2), 4.0 < |X| <= 8.0
    !--------------------------------------------------------------------
    real(wp), PARAMETER :: PY1(9) = [ &
      1.9153806858264202986e+6_wp, -1.1957961912070617006e+9_wp, &
      3.7453673962438488783e+11_wp, -5.9530713129741981618e+13_wp, &
      4.0686275289804744814e+15_wp, -2.3638408497043134724e+16_wp, &
      -5.6808094574724204577e+18_wp, 1.1514276357909013326e+19_wp, &
      -1.2337180442012953128e+3_wp &
    ]
    real(wp), PARAMETER :: QY1(8) = [ &
      1.2855164849321609336e+3_wp, 1.0453748201934079734e+6_wp, &
      6.3550318087088919566e+8_wp, 3.0221766852960403645e+11_wp, &
      1.1187010065856971027e+14_wp, 3.0837179548112881950e+16_wp, &
      5.6968198822857178911e+18_wp, 5.3321844313316185697e+20_wp &
    ]
    !-------------------------------------------------------------------
    ! Coefficients for Hart's approximation, |X| > 8.0
    !-------------------------------------------------------------------
    real(wp), PARAMETER :: P0(6) = [ &
      -1.0982405543459346727e+5_wp, -1.5235293511811373833e+6_wp, &
      -6.6033732483649391093e+6_wp, -9.9422465050776411957e+6_wp, &
      -4.4357578167941278571e+6_wp, -1.6116166443246101165e+3_wp &
    ]
    real(wp), PARAMETER :: Q0(6) = [ &
      -1.0726385991103820119e+5_wp,-1.5118095066341608816e+6_wp, &
      -6.5853394797230870728e+6_wp,-9.9341243899345856590e+6_wp, &
      -4.4357578167941278568e+6_wp,-1.4550094401904961825e+3_wp &
    ]
    real(wp), PARAMETER :: P1(6) = [ &
      1.7063754290207680021e+3_wp, 1.8494262873223866797e+4_wp, &
      6.6178836581270835179e+4_wp, 8.5145160675335701966e+4_wp, &
      3.3220913409857223519e+4_wp, 3.5265133846636032186e+1_wp &
    ]
    real(wp), PARAMETER :: Q1(6) = [ &
      3.7890229745772202641e+4_wp, 4.0029443582266975117e+5_wp, &
      1.4194606696037208929e+6_wp, 1.8194580422439972989e+6_wp, &
      7.0871281941028743574e+5_wp, 8.6383677696049909675e+2_wp &
    ]
    !-------------------------------------------------------------------
    ! Check for error conditions
    !-------------------------------------------------------------------
    ax = ABS(arg)
    IF (jint == 1 .AND. (arg <= zero .OR. (arg < half .AND. ax*xinf < pi2))) THEN
      result = -xinf
      GO TO 80
    ELSE IF (ax > xmax) THEN
      result = zero
      GO TO 80
    END IF
    IF (ax > eight) THEN
      GO TO 60
    ELSE IF (ax <= xsmall) THEN
      IF (jint == 0) THEN
        result = arg * half
      ELSE
        result = -pi2 / ax
      END IF
      GO TO 80
    END IF
    !-------------------------------------------------------------------
    ! Calculate J1 for appropriate interval, preserving accuracy near
    ! the zero of J1
    !-------------------------------------------------------------------
    zsq = ax * ax
    IF (ax <= four) THEN
      xnum = (pj0(7)*zsq+pj0(6)) * zsq + pj0(5)
      xden = zsq + qj0(5)
      DO  i = 1, 4
        xnum = xnum * zsq + pj0(i)
        xden = xden * zsq + qj0(i)
      END DO
      prod = arg * ((ax-xj01/two56)-xj02) * (ax+xj0)
    ELSE
      xnum = pj1(1)
      xden = (zsq+qj1(7)) * zsq + qj1(1)
      DO  i = 2, 6
        xnum = xnum * zsq + pj1(i)
        xden = xden * zsq + qj1(i)
      END DO
      xnum = xnum * (ax-eight) * (ax+eight) + pj1(7)
      xnum = xnum * (ax-four) * (ax+four) + pj1(8)
      prod = arg * ((ax-xj11/two56)-xj12) * (ax+xj1)
    END IF
    result = prod * (xnum/xden)
    IF (jint == 0) GO TO 80
    !-------------------------------------------------------------------
    ! Calculate Y1.  First find  RESJ = pi/2 ln(x/xn) J1(x),
    ! where xn is a zero of Y1
    !-------------------------------------------------------------------
    IF (ax <= four) THEN
      up = (ax-xy01/two56) - xy02
      xy = xy0
    ELSE
      up = (ax-xy11/two56) - xy12
      xy = xy1
    END IF
    down = ax + xy
    IF (ABS(up) < p17*down) THEN
      w = up / down
      wsq = w * w
      xnum = plg(1)
      xden = wsq + qlg(1)
      DO  i = 2, 4
        xnum = xnum * wsq + plg(i)
        xden = xden * wsq + qlg(i)
      END DO
      resj = pi2 * result * w * xnum / xden
    ELSE
      resj = pi2 * result * LOG(ax/xy)
    END IF
    !-------------------------------------------------------------------
    ! Now calculate Y1 for appropriate interval, preserving
    ! accuracy near the zero of Y1
    !-------------------------------------------------------------------
    IF (ax <= four) THEN
      xnum = py0(7) * zsq + py0(1)
      xden = zsq + qy0(1)
      DO  i = 2, 6
        xnum = xnum * zsq + py0(i)
        xden = xden * zsq + qy0(i)
      END DO
    ELSE
      xnum = py1(9) * zsq + py1(1)
      xden = zsq + qy1(1)
      DO  i = 2, 8
        xnum = xnum * zsq + py1(i)
        xden = xden * zsq + qy1(i)
      END DO
    END IF
    result = resj + (up*down/ax) * xnum / xden
    GO TO 80
    !-------------------------------------------------------------------
    ! Calculate J1 or Y1 for |ARG| > 8.0
    !-------------------------------------------------------------------
  60 z = eight / ax
    w = AINT(ax/twopi) + throv8
    w = (ax-w*twopi1) - w * twopi2
    zsq = z * z
    xnum = p0(6)
    xden = zsq + q0(6)
    up = p1(6)
    down = zsq + q1(6)
    DO  i = 1, 5
      xnum = xnum * zsq + p0(i)
      xden = xden * zsq + q0(i)
      up = up * zsq + p1(i)
      down = down * zsq + q1(i)
    END DO
    r0 = xnum / xden
    r1 = up / down
    IF (jint == 0) THEN
      result = (rtpi2/SQRT(ax)) * (r0*COS(w)-z*r1*SIN(w))
    ELSE
      result = (rtpi2/SQRT(ax)) * (r0*SIN(w)+z*r1*COS(w))
    END IF
    IF ((jint == 0) .AND. (arg < zero)) result = -result
  80 RETURN
  END SUBROUTINE caljy1

end module calgo_715
