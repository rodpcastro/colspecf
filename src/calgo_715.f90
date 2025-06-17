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
!
! ## References
! 1. W. J. Cody. 1993. Algorithm 715: SPECFUN–a portable FORTRAN package of special
!    function routines and test drivers. ACM Trans. Math. Softw. 19, 1 (March 1993),
!*   22–30. <https://doi.org/10.1145/151271.151273>

  use csf_kinds, only: wp, dp

  implicit none
  private
  public :: caljy0, caljy1

contains

SUBROUTINE caljy0(arg, result, jint)
! This packet computes zero-order Bessel functions of the first and
! second kind (J0 and Y0), for real arguments X, where 0 < X <= XMAX
! for Y0, and |X| <= XMAX for J0.  It contains two function-type
! subprograms,  BESJ0  and  BESY0,  and one subroutine-type
! subprogram,  CALJY0.  The calling statements for the primary
! entries are:

!           Y = BESJ0(X)
!   and
!           Y = BESY0(X),

!   where the entry points correspond to the functions J0(X) and Y0(X),
!   respectively.  The routine  CALJY0  is intended for internal packet
!   use only, all computations within the packet being concentrated in
!   this one routine.  The function subprograms invoke  CALJY0  with
!   the statement
!           CALL CALJY0(ARG,RESULT,JINT),
!   where the parameter usage is as follows:

!      Function                  Parameters for CALJY0
!       call              ARG             RESULT          JINT

!     BESJ0(ARG)     |ARG| .LE. XMAX       J0(ARG)          0
!     BESY0(ARG)   0 .LT. ARG .LE. XMAX    Y0(ARG)          1

!   The main computation uses unpublished minimax rational
!   approximations for X .LE. 8.0, and an approximation from the
!   book  Computer Approximations  by Hart, et. al., Wiley and Sons,
!   New York, 1968, for arguments larger than 8.0   Part of this
!   transportable packet is patterned after the machine-dependent
!   FUNPACK program BESJ0(X), but cannot match that version for
!   efficiency or accuracy.  This version uses rational functions
!   that are theoretically accurate to at least 18 significant decimal
!   digits for X <= 8, and at least 18 decimal places for X > 8.  The
!   accuracy achieved depends on the arithmetic system, the compiler,
!   the intrinsic functions, and proper selection of the machine-
!   dependent constants.

!*******************************************************************

! The following machine-dependent constants must be declared in
!   DATA statements.  IEEE values are provided as a default.

!   XINF   = largest positive machine number
!   XMAX   = largest acceptable argument.  The functions AINT, SIN
!            and COS must perform properly for  ABS(X) .LE. XMAX.
!            We recommend that XMAX be a small integer multiple of
!            sqrt(1/eps), where eps is the smallest positive number
!            such that  1+eps > 1.
!   XSMALL = positive argument such that  1.0-(X/2)**2 = 1.0
!            to machine precision for all  ABS(X) .LE. XSMALL.
!            We recommend that  XSMALL < sqrt(eps)/beta, where beta
!            is the floating-point radix (usually 2 or 16).

!     Approximate values for some important machines are

!                          eps      XMAX     XSMALL      XINF

!  CDC 7600      (S.P.)  7.11E-15  1.34E+08  2.98E-08  1.26E+322
!  CRAY-1        (S.P.)  7.11E-15  1.34E+08  2.98E-08  5.45E+2465
!  IBM PC (8087) (S.P.)  5.96E-08  8.19E+03  1.22E-04  3.40E+38
!  IBM PC (8087) (D.P.)  1.11D-16  2.68D+08  3.72D-09  1.79D+308
!  IBM 195       (D.P.)  2.22D-16  6.87D+09  9.09D-13  7.23D+75
!  UNIVAC 1108   (D.P.)  1.73D-18  4.30D+09  2.33D-10  8.98D+307
!  VAX 11/780    (D.P.)  1.39D-17  1.07D+09  9.31D-10  1.70D+38

!*******************************************************************

! Error Returns

!  The program returns the value zero for  X .GT. XMAX, and returns
!    -XINF when BESLY0 is called with a negative or zero argument.


! Intrinsic functions required are:

!     ABS, AINT, COS, LOG, SIN, SQRT


!  Latest modification: March 13, 1992

!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439, USA

!--------------------------------------------------------------------

REAL (dp), INTENT(IN)   :: arg
REAL (dp), INTENT(OUT)  :: result
INTEGER, INTENT(IN)     :: jint

INTEGER    :: i
REAL (dp)  :: ax, down, prod, resj, r0, r1, up, w, wsq,  &
              xden, xnum, xy, z, zsq
!-------------------------------------------------------------------
!  Mathematical constants
!    CONS = ln(.5) + Euler's gamma
!-------------------------------------------------------------------
REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, three = 3.0_dp,  &
                         four = 4.0_dp, eight = 8.0_dp, five5 = 5.5_dp,  &
                         sixty4 = 64.0_dp, oneov8 = 0.125_dp,    &
                         p17 = 0.1716_dp, two56 = 256.0_dp,  &
                         cons = -1.1593151565841244881D-1,  &
                         pi2 = 6.3661977236758134308D-1,  &
                         twopi = 6.2831853071795864769_dp,  &
                         twopi1 = 6.28125_dp, twopi2 = 1.9353071795864769253D-3
!-------------------------------------------------------------------
!  Machine-dependent constants
!-------------------------------------------------------------------
REAL (dp), PARAMETER  :: XMAX = 2.68D+08, XSMALL = 3.72D-09, XINF = 1.79D+308
!-------------------------------------------------------------------
!  Zeroes of Bessel functions
!-------------------------------------------------------------------
REAL (dp), PARAMETER  :: XJ0 = 2.4048255576957727686_dp,   &
                         XJ1 = 5.5200781102863106496_dp,   &
                         XY0 = 8.9357696627916752158D-1,   &
                         XY1 = 3.9576784193148578684_dp,   &
                         XY2 = 7.0860510603017726976_dp,   &
                         XJ01 =  616.0_dp, XJ02 = -1.4244423042272313784D-03, &
                         XJ11 = 1413.0_dp, XJ12 =  5.4686028631064959660D-04, &
                         XY01 =  228.0_dp, XY02 =  2.9519662791675215849D-03, &
                         XY11 = 1013.0_dp, XY12 =  6.4716931485786837568D-04, &
                         XY21 = 1814.0_dp, XY22 =  1.1356030177269762362D-04
!-------------------------------------------------------------------
!  Coefficients for rational approximation to ln(x/a)
!--------------------------------------------------------------------
REAL (dp), PARAMETER  :: PLG(4) = (/ -2.4562334077563243311D+01,  &
     2.3642701335621505212D+02, -5.4989956895857911039D+02,  &
     3.5687548468071500413D+02 /)
REAL (dp), PARAMETER  :: QLG(4) = (/ -3.5553900764052419184D+01,  &
     1.9400230218539473193D+02, -3.3442903192607538956D+02,  &
     1.7843774234035750207D+02 /)
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!  J0(X) / (X**2 - XJ0**2),  XSMALL  <  |X|  <=  4.0
!--------------------------------------------------------------------
REAL (dp), PARAMETER  :: PJ0(7) =  &
           (/ 6.6302997904833794242D+06, -6.2140700423540120665D+08,  &
              2.7282507878605942706D+10, -4.1298668500990866786D+11,  &
             -1.2117036164593528341D-01,  1.0344222815443188943D+02,  &
             -3.6629814655107086448D+04 /)
REAL (dp), PARAMETER  :: QJ0(5) =  &
           (/ 4.5612696224219938200D+05, 1.3985097372263433271D+08,  &
              2.6328198300859648632D+10, 2.3883787996332290397D+12,  &
              9.3614022392337710626D+02 /)
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!  J0(X) / (X**2 - XJ1**2),  4.0  <  |X|  <=  8.0
!-------------------------------------------------------------------
REAL (dp), PARAMETER  :: PJ1(8) =  &
           (/ 4.4176707025325087628D+03,  1.1725046279757103576D+04,  &
              1.0341910641583726701D+04, -7.2879702464464618998D+03,  &
             -1.2254078161378989535D+04, -1.8319397969392084011D+03,  &
              4.8591703355916499363D+01,  7.4321196680624245801D+02 /)
REAL (dp), PARAMETER  :: QJ1(7) =  &
           (/ 3.3307310774649071172D+02, -2.9458766545509337327D+03,  &
              1.8680990008359188352D+04, -8.4055062591169562211D+04,  &
              2.4599102262586308984D+05, -3.5783478026152301072D+05,  &
             -2.5258076240801555057D+01 /)
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y0(X) - 2 LN(X/XY0) J0(X)) / (X**2 - XY0**2),
!        XSMALL  <  |X|  <=  3.0
!--------------------------------------------------------------------
REAL (dp), PARAMETER  :: PY0(6) =  &
           (/ 1.0102532948020907590D+04, -2.1287548474401797963D+06,  &
              2.0422274357376619816D+08, -8.3716255451260504098D+09,  &
              1.0723538782003176831D+11, -1.8402381979244993524D+01 /)
REAL (dp), PARAMETER  :: QY0(5) =  &
           (/ 6.6475986689240190091D+02, 2.3889393209447253406D+05,  &
              5.5662956624278251596D+07, 8.1617187777290363573D+09,  &
              5.8873865738997033405D+11 /)
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y0(X) - 2 LN(X/XY1) J0(X)) / (X**2 - XY1**2),
!        3.0  <  |X|  <=  5.5
!--------------------------------------------------------------------
REAL (dp), PARAMETER  :: PY1(7) =  &
           (/ -1.4566865832663635920D+04,  4.6905288611678631510D+06,  &
              -6.9590439394619619534D+08,  4.3600098638603061642D+10,  &
              -5.5107435206722644429D+11, -2.2213976967566192242D+13,  &
               1.7427031242901594547D+01 /)
REAL (dp), PARAMETER  :: QY1(6) =  &
           (/ 8.3030857612070288823D+02, 4.0669982352539552018D+05,  &
              1.3960202770986831075D+08, 3.4015103849971240096D+10,  &
              5.4266824419412347550D+12, 4.3386146580707264428D+14 /)
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y0(X) - 2 LN(X/XY2) J0(X)) / (X**2 - XY2**2),
!        5.5  <  |X|  <=  8.0
!--------------------------------------------------------------------
REAL (dp), PARAMETER  :: PY2(8) =  &
           (/ 2.1363534169313901632D+04, -1.0085539923498211426D+07,  &
              2.1958827170518100757D+09, -1.9363051266772083678D+11,  &
             -1.2829912364088687306D+11,  6.7016641869173237784D+14,  &
             -8.0728726905150210443D+15, -1.7439661319197499338D+01 /)
REAL (dp), PARAMETER  :: QY2(7) =  &
           (/ 8.7903362168128450017D+02, 5.3924739209768057030D+05,  &
              2.4727219475672302327D+08, 8.6926121104209825246D+10,  &
              2.2598377924042897629D+13, 3.9272425569640309819D+15,  &
              3.4563724628846457519D+17 /)
!-------------------------------------------------------------------
!  Coefficients for Hart,s approximation,  |X| > 8.0
!-------------------------------------------------------------------
REAL (dp), PARAMETER  :: P0(6) =  &
           (/ 3.4806486443249270347D+03, 2.1170523380864944322D+04,  &
              4.1345386639580765797D+04, 2.2779090197304684302D+04,  &
              8.8961548424210455236D-01, 1.5376201909008354296D+02 /)

REAL (dp), PARAMETER  :: Q0(5) =  &
           (/ 3.5028735138235608207D+03, 2.1215350561880115730D+04,  &
              4.1370412495510416640D+04, 2.2779090197304684318D+04,  &
              1.5711159858080893649D+02 /)
REAL (dp), PARAMETER  :: P1(6) =  &
           (/ -2.2300261666214198472D+01, -1.1183429920482737611D+02,  &
              -1.8591953644342993800D+02, -8.9226600200800094098D+01,  &
              -8.8033303048680751817D-03, -1.2441026745835638459D+00 /)
REAL (dp), PARAMETER  :: Q1(5) =  &
           (/ 1.4887231232283756582D+03, 7.2642780169211018836D+03,  &
              1.1951131543434613647D+04, 5.7105024128512061905D+03,  &
              9.0593769594993125859D+01 /)
!-------------------------------------------------------------------
!  Check for error conditions
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
!  Calculate J0 for appropriate interval, preserving
!     accuracy near the zero of J0
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
!  Calculate Y0.  First find  RESJ = pi/2 ln(x/xn) J0(x),
!    where xn is a zero of Y0
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
!  Now calculate Y0 for appropriate interval, preserving
!     accuracy near the zero of Y0
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
!  Calculate J0 or Y0 for |ARG|  >  8.0
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
!---------- Last line of CALJY0 ----------
END SUBROUTINE caljy0

SUBROUTINE caljy1(arg, result, jint)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!---------------------------------------------------------------------

! This packet computes first-order Bessel functions of the first and
!   second kind (J1 and Y1), for real arguments X, where 0 < X <= XMAX
!   for Y1, and |X| <= XMAX for J1.  It contains two function-type
!   subprograms,  BESJ1  and  BESY1,  and one subroutine-type
!   subprogram,  CALJY1.  The calling statements for the primary
!   entries are:

!           Y = BESJ1(X)
!   and
!           Y = BESY1(X),

!   where the entry points correspond to the functions J1(X) and Y1(X),
!   respectively.  The routine  CALJY1  is intended for internal packet
!   use only, all computations within the packet being concentrated in
!   this one routine.  The function subprograms invoke  CALJY1  with
!   the statement
!           CALL CALJY1(ARG,RESULT,JINT),
!   where the parameter usage is as follows:

!      Function                  Parameters for CALJY1
!       call              ARG             RESULT          JINT

!     BESJ1(ARG)     |ARG| .LE. XMAX       J1(ARG)          0
!     BESY1(ARG)   0 .LT. ARG .LE. XMAX    Y1(ARG)          1

!   The main computation uses unpublished minimax rational
!   approximations for X .LE. 8.0, and an approximation from the
!   book  Computer Approximations  by Hart, et. al., Wiley and Sons,
!   New York, 1968, for arguments larger than 8.0   Part of this
!   transportable packet is patterned after the machine-dependent
!   FUNPACK program BESJ1(X), but cannot match that version for
!   efficiency or accuracy.  This version uses rational functions
!   that are theoretically accurate to at least 18 significant decimal
!   digits for X <= 8, and at least 18 decimal places for X > 8.  The
!   accuracy achieved depends on the arithmetic system, the compiler,
!   the intrinsic functions, and proper selection of the machine-
!   dependent constants.

!*******************************************************************

! The following machine-dependent constants must be declared in
!   DATA statements.  IEEE values are provided as a default.

!   XINF   = largest positive machine number
!   XMAX   = largest acceptable argument.  The functions AINT, SIN
!            and COS must perform properly for  ABS(X) .LE. XMAX.
!            We recommend that XMAX be a small integer multiple of
!            sqrt(1/eps), where eps is the smallest positive number
!            such that  1+eps > 1.
!   XSMALL = positive argument such that  1.0-(1/2)(X/2)**2 = 1.0
!            to machine precision for all  ABS(X) .LE. XSMALL.
!            We recommend that  XSMALL < sqrt(eps)/beta, where beta
!            is the floating-point radix (usually 2 or 16).

!     Approximate values for some important machines are

!                          eps      XMAX     XSMALL      XINF

!  CDC 7600      (S.P.)  7.11E-15  1.34E+08  2.98E-08  1.26E+322
!  CRAY-1        (S.P.)  7.11E-15  1.34E+08  2.98E-08  5.45E+2465
!  IBM PC (8087) (S.P.)  5.96E-08  8.19E+03  1.22E-04  3.40E+38
!  IBM PC (8087) (D.P.)  1.11D-16  2.68D+08  3.72D-09  1.79D+308
!  IBM 195       (D.P.)  2.22D-16  6.87D+09  9.09D-13  7.23D+75
!  UNIVAC 1108   (D.P.)  1.73D-18  4.30D+09  2.33D-10  8.98D+307
!  VAX 11/780    (D.P.)  1.39D-17  1.07D+09  9.31D-10  1.70D+38

!*******************************************************************

! Error Returns

!  The program returns the value zero for  X .GT. XMAX, and returns
!    -XINF when BESLY1 is called with a negative or zero argument.


! Intrinsic functions required are:

!     ABS, AINT, COS, LOG, SIN, SQRT


!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439, USA

!  Latest modification: March 13, 1992

!--------------------------------------------------------------------

REAL (dp), INTENT(IN)   :: arg
REAL (dp), INTENT(OUT)  :: result
INTEGER, INTENT(IN)     :: jint

INTEGER    :: i
REAL (dp)  :: ax, down, prod, resj, r0, r1, up, w, wsq, xden, xnum, xy, z, zsq
!-------------------------------------------------------------------
!  Mathematical constants
!-------------------------------------------------------------------
REAL (dp), PARAMETER  :: eight = 8.0_dp, four = 4.0_dp, half = 0.5_dp,  &
                         throv8 = 0.375_dp, p17 = 0.1716_dp,  &
                         pi2 = 6.3661977236758134308D-1, zero = 0.0_dp,  &
                         twopi = 6.2831853071795864769_dp, twopi1 = 6.28125_dp, &
                         TWOPI2 = 1.9353071795864769253D-03, two56 = 256.0_dp, &
                         rtpi2 = 7.9788456080286535588D-1
!-------------------------------------------------------------------
!  Machine-dependent constants
!-------------------------------------------------------------------
REAL (dp), PARAMETER  :: XMAX = 2.68D+08, XSMALL = 3.72D-09, XINF = 1.79D+308
!-------------------------------------------------------------------
!  Zeroes of Bessel functions
!-------------------------------------------------------------------
REAL (dp), PARAMETER  :: XJ0 = 3.8317059702075123156D+0,   &
           XJ1 = 7.0155866698156187535_dp, XY0 = 2.1971413260310170351_dp, &
           XY1 = 5.4296810407941351328_dp, XJ01 =  981.0_dp,  &
           XJ02 = -3.2527979248768438556D-04, XJ11 = 1796.0_dp,   &
           XJ12 = -3.8330184381246462950D-05, XY01 =  562.0_dp,   &
           XY02 =  1.8288260310170351490D-03, XY11 = 1390.0_dp,   &
           XY12 = -6.4592058648672279948D-06
!-------------------------------------------------------------------
!  Coefficients for rational approximation to ln(x/a)
!--------------------------------------------------------------------
REAL (dp), PARAMETER  :: PLG(4) =  &
           (/ -2.4562334077563243311D+01, 2.3642701335621505212D+02,  &
              -5.4989956895857911039D+02, 3.5687548468071500413D+02 /)
REAL (dp), PARAMETER  :: QLG(4) =  &
           (/ -3.5553900764052419184D+01, 1.9400230218539473193D+02,  &
              -3.3442903192607538956D+02, 1.7843774234035750207D+02 /)
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!  J1(X) / (X * (X**2 - XJ0**2)),  XSMALL  <  |X|  <=  4.0
!--------------------------------------------------------------------
REAL (dp), PARAMETER  :: PJ0(7) =  &
           (/ 9.8062904098958257677D+05, -1.1548696764841276794D+08,  &
              6.6781041261492395835D+09, -1.4258509801366645672D+11,  &
             -4.4615792982775076130D+03,  1.0650724020080236441D+01,  &
             -1.0767857011487300348D-02 /)
REAL (dp), PARAMETER  :: QJ0(5) =  &
           (/ 5.9117614494174794095D+05, 2.0228375140097033958D+08,  &
              4.2091902282580133541D+10, 4.1868604460820175290D+12,  &
              1.0742272239517380498D+03 /)
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!  J1(X) / (X * (X**2 - XJ1**2)),  4.0  <  |X|  <=  8.0
!-------------------------------------------------------------------
REAL (dp), PARAMETER  :: PJ1(8) =  &
           (/ 4.6179191852758252280D+00, -7.1329006872560947377D+03,  &
              4.5039658105749078904D+06, -1.4437717718363239107D+09,  &
              2.3569285397217157313D+11, -1.6324168293282543629D+13,  &
              1.1357022719979468624D+14,  1.0051899717115285432D+15 /)
REAL (dp), PARAMETER  :: QJ1(7) =  &
           (/ 1.1267125065029138050D+06, 6.4872502899596389593D+08,  &
              2.7622777286244082666D+11, 8.4899346165481429307D+13,  &
              1.7128800897135812012D+16, 1.7253905888447681194D+18,  &
              1.3886978985861357615D+03 /)
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y1(X) - 2 LN(X/XY0) J1(X)) / (X**2 - XY0**2),
!        XSMALL  <  |X|  <=  4.0
!--------------------------------------------------------------------
REAL (dp), PARAMETER  :: PY0(7) =  &
           (/ 2.2157953222280260820D+05, -5.9157479997408395984D+07,  &
              7.2144548214502560419D+09, -3.7595974497819597599D+11,  &
              5.4708611716525426053D+12,  4.0535726612579544093D+13,  &
             -3.1714424660046133456D+02 /)
REAL (dp), PARAMETER  :: QY0(6) =  &
           (/ 8.2079908168393867438D+02, 3.8136470753052572164D+05,  &
              1.2250435122182963220D+08, 2.7800352738690585613D+10,  &
              4.1272286200406461981D+12, 3.0737873921079286084D+14 /)
!--------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y1(X) - 2 LN(X/XY1) J1(X)) / (X**2 - XY1**2),
!        4.0  <  |X|  <=  8.0
!--------------------------------------------------------------------
REAL (dp), PARAMETER  :: PY1(9) =  &
           (/ 1.9153806858264202986D+06, -1.1957961912070617006D+09,  &
              3.7453673962438488783D+11, -5.9530713129741981618D+13,  &
              4.0686275289804744814D+15, -2.3638408497043134724D+16,  &
             -5.6808094574724204577D+18,  1.1514276357909013326D+19,  &
             -1.2337180442012953128D+03 /)
REAL (dp), PARAMETER  :: QY1(8) =  &
           (/ 1.2855164849321609336D+03, 1.0453748201934079734D+06,  &
              6.3550318087088919566D+08, 3.0221766852960403645D+11,  &
              1.1187010065856971027D+14, 3.0837179548112881950D+16,  &
              5.6968198822857178911D+18, 5.3321844313316185697D+20 /)
!-------------------------------------------------------------------
!  Coefficients for Hart,s approximation,  |X| > 8.0
!-------------------------------------------------------------------
REAL (dp), PARAMETER  :: P0(6) =  &
           (/ -1.0982405543459346727D+05, -1.5235293511811373833D+06,  &
              -6.6033732483649391093D+06, -9.9422465050776411957D+06,  &
              -4.4357578167941278571D+06, -1.6116166443246101165D+03 /)
REAL (dp), PARAMETER  :: Q0(6) =  &
           (/ -1.0726385991103820119D+05,-1.5118095066341608816D+06,  &
              -6.5853394797230870728D+06,-9.9341243899345856590D+06,  &
              -4.4357578167941278568D+06,-1.4550094401904961825D+03 /)
REAL (dp), PARAMETER  :: P1(6) =  &
           (/ 1.7063754290207680021D+03, 1.8494262873223866797D+04,  &
              6.6178836581270835179D+04, 8.5145160675335701966D+04,  &
              3.3220913409857223519D+04, 3.5265133846636032186D+01 /)
REAL (dp), PARAMETER  :: Q1(6) =  &
           (/ 3.7890229745772202641D+04, 4.0029443582266975117D+05,  &
              1.4194606696037208929D+06, 1.8194580422439972989D+06,  &
              7.0871281941028743574D+05, 8.6383677696049909675D+02 /)
!-------------------------------------------------------------------
!  Check for error conditions
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
!  Calculate J1 for appropriate interval, preserving accuracy near
!  the zero of J1
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
!  Calculate Y1.  First find  RESJ = pi/2 ln(x/xn) J1(x),
!    where xn is a zero of Y1
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
!  Now calculate Y1 for appropriate interval, preserving
!     accuracy near the zero of Y1
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
!  Calculate J1 or Y1 for |ARG|  >  8.0
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
!---------- Last card of CALJY1 ----------
END SUBROUTINE caljy1

end module calgo_715
