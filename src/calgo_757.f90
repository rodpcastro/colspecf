! Licensed under the ACM Software License Agreement
! Copyright © 1970–2012 Association for Computing Machinery (ACM)
! See ColSpecF LICENSE file for details.

module calgo_757
!* # CALGO 757
! Algorithm 757. 
!
! Procedures:
!
! - `strvh0`: Struve function \(\mathbf{H}_0(x)\)
! - `strvh1`: Struve function \(\mathbf{H}_1(x)\)
!
! Other CALGO 757 procedures **not** yet included in this module:
!
! - `abram0`, `abram1`, `abram2`: Abramowitz functions \(f_m(x) = 
!     \int_{0}^{\infty} t^m e^{-t^2 -\frac{x}{t}} \, dt, \, \text{for } m = 0, 1, 2\)
! - `airygi`: Modified Airy function \(\mathrm{Gi}(x)\)
! - `airyhi`: Modified Airy function \(\mathrm{Hi}(x)\)
! - `airint`: Integral of the Airy function \(\mathrm{Ai}(x)\) from \(0\) to \(x\)
! - `birint`: Integral of the Airy function \(\mathrm{Bi}(x)\) from \(0\) to \(x\)
! - `j0int`: Integral of the Bessel function \(J_0(x)\) from \(0\) to \(x\)
! - `y0int`: Integral of the Bessel function \(Y_0(x)\) from \(0\) to \(x\)
! - `i0int`: Integral of the modified Bessel function \(I_0(x)\) from \(0\) to \(x\)
! - `k0int`: Integral of the modified Bessel function \(K_0(x)\) from \(0\) to \(x\)
! - `debye1`, `debye2`, `debye3`, `debye4`: Debye functions \(D_n(x) = 
!     \frac{n}{x^n} \int_{0}^{x} \frac{t^n}{e^t-1} \, dt,
!     \, \text{for } n = 1, 2, 3, 4\)
! - `strvl0`: Modified Struve function \(\mathbf{L}_0(x)\)
! - `strvl1`: Modified Struve function \(\mathbf{L}_1(x)\)
! - `i0ml0`: Difference between the modified Bessel and Struve functions
!     \(I_0(x) - \mathbf{L}_0(x)\)
! - `i1ml1`: Difference between the modified Bessel and Struve functions
!     \(I_1(x) - \mathbf{L}_1(x)\)
! - `synch1`: Synchrotron radiation function
!     \(f_1(x) = x \int_{x}^{\infty} K_{5/3}(t) \, dt\)
! - `synch2`: Synchrotron radiation function
!     \(f_2(x) = x K_{2/3}(x)\)
! - `tran02`, `tran03`, ..., `tran09`: Transport integrals \(J_n(x) = 
!     \int_{0}^{x} \frac{t^n e^t}{(e^t-1)^2} \, dt,
!     \, \text{for } n = 2, 3, \ldots, 9\) 
! - `atnint`: Inverse-tangent integral \(\int_{0}^{x} \frac{\arctan t}{t} \, dt\)
! - `clausn`: Clausen's integral
!     \(-\int_{0}^{x} \ln{\lvert 2 \sin\frac{t}{2} \rvert} \, dt\)
! - `exp3`: \(\int_{0}^{x} e^{-t^3} \, dt\)
! - `goodst`: \(\int_{0}^{\infty} \frac{e^{-t^2}}{t + x} \, dt\)
! - `lobach`: Lobachevski's integral
!     \(-\int_{0}^{\infty} \log \lvert \cos t \rvert \, dt\)
! - `strom`: Stromgren's integral
!     \(\frac{15}{4 \pi^4} \int_{0}^{x} \frac{t^7 e^{2 t}}{(e^t - 1)^3} \, dt\)
!
! ## Author
! Allan J. MacLeod
!
! ## History
! - 1996-01-23 - Allan J. MacLeod
!     - Original code.
! - 2000-06-13 - Allan J. MacLeod
!     - F77 code distributed by ACM: <https://calgo.acm.org/757.zip>
! - 2001-10-10 - Alan Miller
!     - F90 code adaptation by Alan Miller:
!       <https://jblevins.org/mirror/amiller/toms757.zip>
! - 2025-07-27 - Rodrigo Castro (GitHub: rodpcastro)
!     - Retained only the subroutines `strvh0` and `strvh1`; additional subroutines
!       will be included as required.
!     - Replaced array constructor `(/.../)` by the less verbose `[...]`.
!     - Replaced `dp` (double precision) by `wp` (working precision).
!     - Replaced `EPSILON(0.0_dp)` by `eps_wp` (CSF working precision machine epsilon).
!     - Removed subroutine `errprn` that prints error messages, which are unnecessary
!       for `strvh0` and `strvh1`.
!
! ## References
! 1. Allan J. MacLeod. 1996. Algorithm 757: MISCFUN, a software package to compute
!    uncommon special functions. ACM Trans. Math. Softw. 22, 3 (Sept. 1996), 288–301.
!*   <https://doi.org/10.1145/232826.232846>

  use csf_kinds, only: wp
  use csf_numerror, only: eps_wp

  implicit none
  private
  public :: strvh0, strvh1

contains

  FUNCTION strvh0(xvalue) RESULT(fn_val)
    !! CALGO 757 Struve function \(\mathbf{H}_0(x)\).
    !
    ! DESCRIPTION:
    !     This function calculates the value of the Struve function
    !     of order 0, denoted H0(x), for the argument XVALUE, defined
    !
    !         STRVHO(x) = (2/pi) integral{0 to pi/2} sin(x cos(t)) dt
    !
    !     H0 also satisfies the second-order equation
    !
    !                 x*D(Df) + Df + x*f = 2x/pi
    !
    !     The code uses Chebyshev expansions whose coefficients are given to 20D.
    !
    !   ERROR RETURNS:
    !     As the asymptotic expansion of H0 involves the Bessel function
    !     of the second kind Y0, there is a problem for large x, since
    !     we cannot accurately calculate the value of Y0. An error message
    !     is printed and STRVH0 returns the value 0.0.
    !
    !   MACHINE-DEPENDENT CONSTANTS:
    !     NTERM1 - The no. of terms to be used in the array ARRH0. The
    !              recommended value is such that
    !                      ABS(ARRH0(NTERM1)) < EPS/100.
    !     NTERM2 - The no. of terms to be used in the array ARRH0A. The
    !              recommended value is such that
    !                      ABS(ARRH0A(NTERM2)) < EPS/100.
    !     NTERM3 - The no. of terms to be used in the array AY0ASP. The
    !              recommended value is such that
    !                      ABS(AY0ASP(NTERM3)) < EPS/100.
    !     NTERM4 - The no. of terms to be used in the array AY0ASQ. The
    !              recommended value is such that
    !                      ABS(AY0ASQ(NTERM4)) < EPS/100.
    !     XLOW - The value for which H0(x) = 2*x/pi to machine precision, if
    !            abs(x) < XLOW. The recommended value is
    !                      XLOW = 3 * SQRT(EPSNEG)
    !     XHIGH - The value above which we are unable to calculate Y0 with
    !             any reasonable accuracy. An error message is printed and
    !              STRVH0 returns the value 0.0. The recommended value is
    !                      XHIGH = 1/EPS.
    !
    !     For values of EPS and EPSNEG refer to the file MACHCON.TXT.

    REAL(wp), INTENT(IN) :: xvalue
    REAL(wp)             :: fn_val

    INTEGER  :: indsgn, nterm1, nterm2, nterm3, nterm4
    REAL(wp) :: h0as, t, x, xhigh, xlow, xmp4, xsq, y0p, y0q, y0val

    REAL(wp), PARAMETER :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp, &
                           eight = 8.0_wp, eleven = 11.0_wp, twenty = 20.0_wp, &
                           onehun = 100.0_wp, sixtp5 = 60.5_wp, &
                           two62 = 262.0_wp, thr2p5 = 302.5_wp, &
                           piby4 = 0.78539816339744830962_wp, &
                           rt2bpi = 0.79788456080286535588_wp, &
                           twobpi = 0.63661977236758134308_wp

    REAL(wp), PARAMETER :: arrh0(0:19) = [ &
      0.28696487399013225740_wp, -0.25405332681618352305_wp, &
      0.20774026739323894439_wp, -0.20364029560386585140_wp, &
      0.12888469086866186016_wp, -0.4825632815622261202e-1_wp, &
      0.1168629347569001242e-1_wp, -0.198118135642418416e-2_wp, &
      0.24899138512421286e-3_wp, -0.2418827913785950e-4_wp, &
      0.187437547993431e-5_wp, -0.11873346074362e-6_wp, &
      0.626984943346e-8_wp, -0.28045546793e-9_wp, &
      0.1076941205e-10_wp, -0.35904793e-12_wp, &
      0.1049447e-13_wp, -0.27119e-15_wp, &
      0.624e-17_wp, -0.13e-18_wp &
    ]
    REAL(wp), PARAMETER :: arrh0a(0:20) = [ &
      1.99291885751992305515_wp, -0.384232668701456887e-2_wp, &
     -0.32871993712353050e-3_wp, -0.2941181203703409e-4_wp, &
     -0.267315351987066e-5_wp, -0.24681031075013e-6_wp, &
     -0.2295014861143e-7_wp, -0.215682231833e-8_wp, &
     -0.20303506483e-9_wp, -0.1934575509e-10_wp, &
     -0.182773144e-11_wp, -0.17768424e-12_wp, &
     -0.1643296e-13_wp, -0.171569e-14_wp, &
     -0.13368e-15_wp, -0.2077e-16_wp, &
      0.2e-19_wp, -0.55e-18_wp, &
      0.10e-18_wp, -0.4e-19_wp, &
      0.1e-19_wp &
    ]
    REAL(wp), PARAMETER :: ay0asp(0:12) = [ &
      1.99944639402398271568_wp, -0.28650778647031958e-3_wp, &
     -0.1005072797437620e-4_wp, -0.35835941002463e-6_wp, &
     -0.1287965120531e-7_wp, -0.46609486636e-9_wp, &
     -0.1693769454e-10_wp, -0.61852269e-12_wp, &
     -0.2261841e-13_wp, -0.83268e-15_wp, &
     -0.3042e-16_wp, -0.115e-17_wp, &
     -0.4e-19_wp &
    ]
    REAL(wp), PARAMETER :: ay0asq(0:13) = [ &
      1.99542681386828604092_wp, -0.236013192867514472e-2_wp, &
     -0.7601538908502966e-4_wp, -0.256108871456343e-5_wp, &
     -0.8750292185106e-7_wp, -0.304304212159e-8_wp, &
     -0.10621428314e-9_wp, -0.377371479e-11_wp, &
     -0.13213687e-12_wp, -0.488621e-14_wp, &
     -0.15809e-15_wp, -0.762e-17_wp, &
     -0.3e-19_wp, -0.3e-19_wp &
    ]

    ! Start computation
    x = xvalue
    indsgn = 1
    IF (x < zero) THEN
      x = -x
      indsgn = -1
    END IF

    ! Compute the machine-dependent constants.
    h0as = eps_wp
    xhigh = one / (2*eps_wp)

    ! Error test
    IF (ABS(xvalue) > xhigh) THEN
      fn_val = zero
      RETURN
    END IF

    ! continue with machine constants
    t = h0as / onehun
    IF (x <= eleven) THEN
      DO  nterm1 = 19, 0, -1
        IF (ABS(arrh0(nterm1)) > t) EXIT
      END DO
      y0p = SQRT(h0as)
      xlow = y0p + y0p + y0p
    ELSE
      DO  nterm2 = 20, 0, -1
        IF (ABS(arrh0a(nterm2)) > t) EXIT
      END DO
      DO  nterm3 = 12, 0, -1
        IF (ABS(ay0asp(nterm3)) > t) EXIT
      END DO
      DO  nterm4 = 13, 0, -1
        IF (ABS(ay0asq(nterm4)) > t) EXIT
      END DO
    END IF

    ! Code for abs(x) <= 11
    IF (x <= eleven) THEN
      IF (x < xlow) THEN
        fn_val = twobpi * x
      ELSE
        t = ((x*x)/sixtp5-half) - half
        fn_val = twobpi * x * cheval(nterm1,arrh0,t)
      END IF
    ELSE
      ! Code for abs(x) > 11
      xsq = x * x
      t = (two62-xsq) / (twenty+xsq)
      y0p = cheval(nterm3,ay0asp,t)
      y0q = cheval(nterm4,ay0asq,t) / (eight*x)
      xmp4 = x - piby4
      y0val = y0p * SIN(xmp4) - y0q * COS(xmp4)
      y0val = y0val * rt2bpi / SQRT(x)
      t = (thr2p5-xsq) / (sixtp5+xsq)
      h0as = twobpi * cheval(nterm2,arrh0a,t) / x
      fn_val = y0val + h0as
    END IF
    IF (indsgn == -1) fn_val = -fn_val
    RETURN
  END FUNCTION strvh0


  FUNCTION strvh1(xvalue) RESULT(fn_val)
    !! CALGO 757 Struve function \(\mathbf{H}_1(x)\).
    !
    ! DESCRIPTION:
    !     This function calculates the value of the Struve function
    !     of order 1, denoted H1(x), for the argument XVALUE, defined as
    !                                                                  2
    !        STRVH1(x) = (2x/pi) integral{0 to pi/2} sin( x cos(t))*sin t dt
    !
    !     H1 also satisfies the second-order differential equation
    !
    !                    2   2                   2            2
    !                   x * D f  +  x * Df  +  (x - 1)f  =  2x / pi
    !
    !     The code uses Chebyshev expansions with the coefficients given to 20D.
    !
    ! ERROR RETURNS:
    !     As the asymptotic expansion of H1 involves the Bessel function
    !     of the second kind Y1, there is a problem for large x, since
    !     we cannot accurately calculate the value of Y1. An error message
    !     is printed and STRVH1 returns the value 0.0.
    !
    !   MACHINE-DEPENDENT CONSTANTS:
    !     NTERM1 - The no. of terms to be used in the array ARRH1. The
    !              recommended value is such that
    !                      ABS(ARRH1(NTERM1)) < EPS/100.
    !     NTERM2 - The no. of terms to be used in the array ARRH1A. The
    !              recommended value is such that
    !                      ABS(ARRH1A(NTERM2)) < EPS/100.
    !     NTERM3 - The no. of terms to be used in the array AY1ASP. The
    !              recommended value is such that
    !                      ABS(AY1ASP(NTERM3)) < EPS/100.
    !     NTERM4 - The no. of terms to be used in the array AY1ASQ. The
    !              recommended value is such that
    !                      ABS(AY1ASQ(NTERM4)) < EPS/100.
    !     XLOW1 - The value of x, below which H1(x) set to zero, if
    !             abs(x)<XLOW1. The recommended value is
    !                      XLOW1 = 1.5 * SQRT(XMIN)
    !     XLOW2 - The value for which H1(x) = 2*x*x/pi to machine precision, if
    !             abs(x) < XLOW2. The recommended value is
    !                      XLOW2 = SQRT(15*EPSNEG)
    !     XHIGH - The value above which we are unable to calculate Y1 with
    !             any reasonable accuracy. An error message is printed and
    !             STRVH1 returns the value 0.0.  The recommended value is
    !                      XHIGH = 1/EPS.
    !
    !     For values of EPS, EPSNEG and XMIN refer to the file MACHCON.TXT.

    REAL(wp), INTENT(IN) :: xvalue
    REAL(wp)             :: fn_val

    INTEGER  :: nterm1, nterm2, nterm3, nterm4
    REAL(wp) :: h1as, t, x, xhigh, xlow1, xlow2, xm3p4, xsq, y1p, y1q, y1val

    REAL(wp), PARAMETER :: zero = 0.0_wp, half = 0.5_wp, eight = 8.0_wp, &
                           nine = 9.0_wp, fiften = 15.0_wp, twenty = 20.0_wp, &
                           onehun = 100.0_wp, fortp5 = 40.5_wp, &
                           one82 = 182.0_wp, tw02p5 = 202.5_wp, &
                           rt2bpi = 0.79788456080286535588_wp, &
                           thpby4 = 2.35619449019234492885_wp, &
                           twobpi = 0.63661977236758134308_wp

    REAL(wp), PARAMETER :: arrh1(0:17) = [ &
      0.17319061083675439319_wp, -0.12606917591352672005_wp, &
      0.7908576160495357500e-1_wp, -0.3196493222321870820e-1_wp, &
      0.808040581404918834e-2_wp, -0.136000820693074148e-2_wp, &
      0.16227148619889471e-3_wp, -0.1442352451485929e-4_wp, &
      0.99219525734072e-6_wp, -0.5441628049180e-7_wp, &
      0.243631662563e-8_wp, -0.9077071338e-10_wp, &
      0.285926585e-11_wp, -0.7716975e-13_wp, &
      0.180489e-14_wp, -0.3694e-16_wp, &
      0.67e-18_wp, -0.1e-19_wp &
    ]
    REAL(wp), PARAMETER  :: arrh1a(0:21) = [ &
      2.01083504951473379407_wp, 0.592218610036099903e-2_wp, &
      0.55274322698414130e-3_wp, 0.5269873856311036e-4_wp, &
      0.506374522140969e-5_wp, 0.49028736420678e-6_wp, &
      0.4763540023525e-7_wp, 0.465258652283e-8_wp, &
      0.45465166081e-9_wp, 0.4472462193e-10_wp, &
      0.437308292e-11_wp, 0.43568368e-12_wp, &
      0.4182190e-13_wp, 0.441044e-14_wp, &
      0.36391e-15_wp, 0.5558e-16_wp, &
     -0.4e-19_wp, 0.163e-17_wp, &
     -0.34e-18_wp, 0.13e-18_wp, &
     -0.4e-19_wp, 0.1e-19_wp &
    ]
    REAL(wp), PARAMETER  :: ay1asp(0:14) = [ &
      2.00135240045889396402_wp, 0.71104241596461938e-3_wp, &
      0.3665977028232449e-4_wp, 0.191301568657728e-5_wp, &
      0.10046911389777e-6_wp, 0.530401742538e-8_wp, &
      0.28100886176e-9_wp, 0.1493886051e-10_wp, &
      0.79578420e-12_wp, 0.4252363e-13_wp, &
      0.227195e-14_wp, 0.12216e-15_wp, &
      0.650e-17_wp, 0.36e-18_wp, &
      0.2e-19_wp &
    ]
    REAL(wp), PARAMETER  :: ay1asq(0:15) = [ &
      5.99065109477888189116_wp, -0.489593262336579635e-2_wp, &
     -0.23238321307070626e-3_wp, -0.1144734723857679e-4_wp, &
     -0.57169926189106e-6_wp, -0.2895516716917e-7_wp, &
     -0.147513345636e-8_wp, -0.7596537378e-10_wp, &
     -0.390658184e-11_wp, -0.20464654e-12_wp, &
     -0.1042636e-13_wp, -0.57702e-15_wp, &
     -0.2550e-16_wp, -0.210e-17_wp, &
      0.2e-19_wp, -0.2e-19_wp &
    ]

    ! Start computation
    x = ABS(xvalue)

    ! Compute the machine-dependent constants.
    xhigh = (half+half) / (2*eps_wp)

    ! Error test
    IF (x > xhigh) THEN
      fn_val = zero
      RETURN
    END IF

    ! continue with machine constants
    h1as = eps_wp
    t = h1as / onehun
    IF (x <= nine) THEN
      DO  nterm1 = 17, 0, -1
        IF (ABS(arrh1(nterm1)) > t) EXIT
      END DO
      xlow1 = half * SQRT(TINY(zero))
      xlow1 = xlow1 + xlow1 + xlow1
      xlow2 = SQRT(fiften*h1as)
    ELSE
      DO  nterm2 = 21, 0, -1
        IF (ABS(arrh1a(nterm2)) > t) EXIT
      END DO
      DO  nterm3 = 14, 0, -1
        IF (ABS(ay1asp(nterm3)) > t) EXIT
      END DO
      DO  nterm4 = 15, 0, -1
        IF (ABS(ay1asq(nterm4)) > t) EXIT
      END DO
    END IF

    ! Code for abs(x) <= 9
    IF (x <= nine) THEN
      IF (x < xlow1) THEN
        fn_val = zero
      ELSE
        xsq = x * x
        IF (x < xlow2) THEN
          fn_val = twobpi * xsq
        ELSE
          t = (xsq/fortp5-half) - half
          fn_val = twobpi * xsq * cheval(nterm1,arrh1,t)
        END IF
      END IF
    ELSE
      ! Code for abs(x) > 9
      xsq = x * x
      t = (one82-xsq) / (twenty+xsq)
      y1p = cheval(nterm3,ay1asp,t)
      y1q = cheval(nterm4,ay1asq,t) / (eight*x)
      xm3p4 = x - thpby4
      y1val = y1p * SIN(xm3p4) + y1q * COS(xm3p4)
      y1val = y1val * rt2bpi / SQRT(x)
      t = (tw02p5-xsq) / (fortp5+xsq)
      h1as = twobpi * cheval(nterm2,arrh1a,t)
      fn_val = y1val + h1as
    END IF
    RETURN
  END FUNCTION strvh1


  FUNCTION cheval(n, a, t) RESULT(fn_val)
    ! This function evaluates a Chebyshev series, using the Clenshaw method
    ! with Reinsch modification, as analysed in the paper by Oliver.
    !
    ! INPUT PARAMETERS
    !     N - INTEGER - The no. of terms in the sequence
    !     A - REAL(wp) ARRAY, dimension 0 to N - The coefficients of
    !         the Chebyshev series
    !     T - REAL(wp) - The value at which the series is to be evaluated
    !
    ! REFERENCES
    !      "An error analysis of the modified Clenshaw method for
    !       evaluating Chebyshev and Fourier series" J. Oliver,
    !       J.I.M.A., vol. 20, 1977, pp379-391

    INTEGER, INTENT(IN)  :: n
    REAL(wp), INTENT(IN) :: a(0:n)
    REAL(wp), INTENT(IN) :: t
    REAL(wp)             :: fn_val
    INTEGER  :: i
    REAL(wp) :: d1, d2, tt, u0, u1, u2
    REAL(wp), PARAMETER :: zero = 0.0_wp, half = 0.5_wp, test = 0.6_wp, two = 2.0_wp

    u1 = zero

    ! If ABS(T) < 0.6 use the standard Clenshaw method
    IF (ABS(t) < test) THEN
      u0 = zero
      tt = t + t
      DO  i = n, 0, -1
        u2 = u1
        u1 = u0
        u0 = tt * u1 + a(i) - u2
      END DO
      fn_val = (u0-u2) / two
    ELSE
      ! If ABS(T) >= 0.6 use the Reinsch modification
      d1 = zero
      
      ! T >= 0.6 code
      IF (t > zero) THEN
        tt = (t-half) - half
        tt = tt + tt
        DO  i = n, 0, -1
          d2 = d1
          u2 = u1
          d1 = tt * u2 + a(i) + d2
          u1 = d1 + u2
        END DO
        fn_val = (d1+d2) / two
      ELSE
        ! T <= -0.6 code
        tt = (t+half) + half
        tt = tt + tt
        DO  i = n, 0, -1
          d2 = d1
          u2 = u1
          d1 = tt * u2 + a(i) - d2
          u1 = d1 - u2
        END DO
        fn_val = (d1-d2) / two
      END IF
    END IF
    RETURN
  END FUNCTION cheval

end module calgo_757
