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
! Untested procedures:  
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
!     - Replaced array constructor `(/.../)` by the less verbose `[...]`.
!     - Replaced `dp` (double precision) by `wp` (working precision).
!
! ## References
! 1. Allan J. MacLeod. 1996. Algorithm 757: MISCFUN, a software package to compute
!    uncommon special functions. ACM Trans. Math. Softw. 22, 3 (Sept. 1996), 288–301.
!*   <https://doi.org/10.1145/232826.232846>

  implicit none
  private
  public :: abram0, abram1, abram2, &
            airygi, airyhi, airint, birint, &
            j0int, y0int, i0int, k0int, &
            debye1, debye2, debye3, debye4, &
            strvh0, strvh1, strvl0, strvl1, i0ml0, i1ml1, &
            synch1, synch2, &
            tran02, tran03, tran04, tran05, tran06, tran07, tran08, tran09, &
            atnint, clausn, exp3, goodst, lobach, strom

  use csf_kinds, only: wp

contains

  FUNCTION abram0(xvalue) RESULT(fn_val)
    !! CALGO 757 Abramowitz function \(f_0(x) = 
    !! \int_{0}^{\infty} e^{-t^2 -\frac{x}{t}} \, dt\)
    !
    ! DESCRIPTION:
    !    This function calculates the Abramowitz function of order 0, defined as
    !
    !     ABRAM0(x) = integral{ 0 to infinity } exp( -t*t - x/t ) dt
    !
    !     The code uses Chebyshev expansions with the coefficients
    !     given to an accuracy of 20 decimal places.
    !
    ! ERROR RETURNS:
    !    If XVALUE < 0.0, the function prints a message and returns the value 0.0.
    !
    ! MACHINE-DEPENDENT CONSTANTS:
    !    NTERMF - INTEGER - No. of terms needed for the array AB0F.
    !             Recommended value such that
    !                   ABS( AB0F(NTERMF) ) < EPS/100
    !    NTERMG - INTEGER - No. of terms needed for array AB0G.
    !             Recommended value such that
    !                   ABS( AB0G(NTERMG) ) < EPS/100
    !    NTERMH - INTEGER - No. of terms needed for array AB0H.
    !             Recommended value such that
    !                   ABS( AB0H(NTERMH) ) < EPS/100
    !    NTERMA - INTEGER - No. of terms needed for array AB0AS.
    !             Recommended value such that
    !                   ABS( AB0AS(NTERMA) ) < EPS/100
    !   XLOW1 - REAL(wp) - The value below which
    !            ABRAM0 = root(pi)/2 + X ( ln X - GVAL0 )
    !           Recommended value is SQRT(2*EPSNEG)
    !   LNXMIN - REAL(wp) - The value of ln XMIN. Used to prevent
    !            exponential underflow for large X.
    !
    !   For values of EPS, EPSNEG, XMIN refer to the file MACHCON.TXT.
    !
    !   The machine-dependent constants are computed internally by
    !   using the D1MACH subroutine.
    !
    ! INTRINSIC FUNCTIONS USED:
    !   LOG, EXP, SQRT
    !
    ! OTHER MISCFUN SUBROUTINES USED:
    !        CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN) :: xvalue
    REAL(wp)             :: fn_val
    INTEGER  :: nterma, ntermf, ntermg, ntermh
    REAL(wp) :: asln, asval, fval, gval, hval, lnxmin, t, v, x, xlow1
    CHARACTER (LEN=33) :: errmsg = 'FUNCTION CALLED WITH ARGUMENT < 0'
    CHARACTER (LEN=6)  :: fnname = 'ABRAM0'

    REAL(wp), PARAMETER :: ab0f(0:8) = [ &
     -0.68121927093549469816_wp, -0.78867919816149252495_wp, &
      0.5121581776818819543e-1_wp, -0.71092352894541296e-3_wp, &
      0.368681808504287e-5_wp, -0.917832337237e-8_wp, 0.1270202563e-10_wp, &
     -0.1076888e-13_wp, 0.599e-17_wp]
    REAL(wp), PARAMETER :: ab0g(0:8) = [ &
     -0.60506039430868273190_wp, -0.41950398163201779803_wp, &
      0.1703265125190370333e-1_wp, -0.16938917842491397e-3_wp, &
      0.67638089519710e-6_wp, -0.135723636255e-8_wp, 0.156297065e-11_wp, &
     -0.112887e-14_wp, 0.55e-18_wp]
    REAL(wp), PARAMETER :: ab0h(0:8) = [ &
      1.38202655230574989705_wp, -0.30097929073974904355_wp, &
      0.794288809364887241e-2_wp, -0.6431910276847563e-4_wp, &
      0.22549830684374e-6_wp, -0.41220966195e-9_wp, 0.44185282e-12_wp, &
     -0.30123e-15_wp, 0.14e-18_wp]
    REAL(wp), PARAMETER :: ab0as(0:27) = [ &
      1.97755499723693067407_wp, -0.1046024792004819485e-1_wp, &
      0.69680790253625366e-3_wp, -0.5898298299996599e-4_wp, &
      0.577164455305320e-5_wp, -0.61523013365756e-6_wp, 0.6785396884767e-7_wp, &
     -0.723062537907e-8_wp, 0.63306627365e-9_wp, -0.989453793e-11_wp, &
     -0.1681980530e-10_wp, 0.673799551e-11_wp, -0.200997939e-11_wp, &
      0.54055903e-12_wp, -0.13816679e-12_wp, 0.3422205e-13_wp, -0.826686e-14_wp, &
      0.194566e-14_wp, -0.44268e-15_wp, 0.9562e-16_wp, -0.1883e-16_wp, 0.301e-17_wp, &
     -0.19e-18_wp, -0.14e-18_wp, 0.11e-18_wp, -0.4e-19_wp, 0.2e-19_wp, -0.1e-19_wp]

    REAL(wp), PARAMETER :: zero = 0.0_wp, half = 0.5_wp, two = 2.0_wp, &
                           three = 3.0_wp, six = 6.0_wp, onehun = 100.0_wp, &
                           rt3bpi = 0.97720502380583984317_wp, &
                           rtpib2 = 0.88622692545275801365_wp, &
                           gval0 = 0.13417650264770070909_wp, &
                           onerpi = 0.56418958354775628695_wp

    ! Start computation
    x = xvalue

    ! Error test
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    t = 2*EPSILON(zero) / onehun
    IF (x <= two) THEN
      DO  ntermf = 8, 0, -1
        IF (ABS(ab0f(ntermf)) > t) EXIT
      END DO
      DO  ntermg = 8, 0, -1
        IF (ABS(ab0g(ntermg)) > t) EXIT
      END DO
      DO  ntermh = 8, 0, -1
        IF (ABS(ab0h(ntermh)) > t) EXIT
      END DO
      xlow1 = SQRT(two*EPSILON(zero))
    ELSE
      DO  nterma = 27, 0, -1
        IF (ABS(ab0as(nterma)) > t) EXIT
      END DO
      lnxmin = LOG(TINY(zero))
    END IF

    ! Code for 0 <= XVALUE <= 2
    IF (x <= two) THEN
      IF (x == zero) THEN
        fn_val = rtpib2
        RETURN
      END IF
      IF (x < xlow1) THEN
        fn_val = rtpib2 + x * (LOG(x)-gval0)
        RETURN
      ELSE
        t = (x*x/two-half) - half
        fval = cheval(ntermf,ab0f,t)
        gval = cheval(ntermg,ab0g,t)
        hval = cheval(ntermh,ab0h,t)
        fn_val = fval / onerpi + x * (LOG(x)*hval-gval)
        RETURN
      END IF
    ELSE
      ! Code for XVALUE > 2
      v = three * ((x/two)**(two/three))
      t = (six/v-half) - half
      asval = cheval(nterma,ab0as,t)
      asln = LOG(asval/rt3bpi) - v
      IF (asln < lnxmin) THEN
        fn_val = zero
      ELSE
        fn_val = EXP(asln)
      END IF
    END IF
    RETURN
  END FUNCTION abram0


  FUNCTION abram1(xvalue) RESULT(fn_val)
    !! CALGO 757 Abramowitz function \(f_1(x) = 
    !! \int_{0}^{\infty} t e^{-t^2 -\frac{x}{t}} \, dt\)
    !
    !   DESCRIPTION:
    !      This function calculates the Abramowitz function of order 1, defined as
    !
    !       ABRAM1(x) = integral{ 0 to infinity } t * exp( -t*t - x/t ) dt
    !
    !       The code uses Chebyshev expansions with the coefficients
    !       given to an accuracy of 20 decimal places.
    !
    !   ERROR RETURNS:
    !      If XVALUE < 0.0, the function prints a message and returns the value 0.0.
    !
    !   MACHINE-DEPENDENT CONSTANTS:
    !      NTERMF - INTEGER - No. of terms needed for the array AB1F.
    !               Recommended value such that
    !                     ABS( AB1F(NTERMF) ) < EPS/100
    !      NTERMG - INTEGER - No. of terms needed for array AB1G.
    !               Recommended value such that
    !                     ABS( AB1G(NTERMG) ) < EPS/100
    !      NTERMH - INTEGER - No. of terms needed for array AB1H.
    !               Recommended value such that
    !                     ABS( AB1H(NTERMH) ) < EPS/100
    !      NTERMA - INTEGER - No. of terms needed for array AB1AS.
    !               Recommended value such that
    !                     ABS( AB1AS(NTERMA) ) < EPS/100
    !      XLOW - REAL(wp) - The value below which
    !                ABRAM1(x) = 0.5 to machine precision.
    !             The recommended value is EPSNEG/2
    !      XLOW1 - REAL(wp) - The value below which
    !                ABRAM1(x) = (1 - x ( sqrt(pi) + xln(x) ) / 2
    !              Recommended value is SQRT(2*EPSNEG)
    !      LNXMIN - REAL(wp) - The value of ln XMIN. Used to prevent
    !              exponential underflow for large X.
    !
    !      For values of EPS, EPSNEG, XMIN refer to the file MACHCON.TXT
    !
    !      The machine-dependent constants are computed internally by using
    !      the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !     LOG, EXP, SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN) :: xvalue
    REAL(wp)             :: fn_val

    INTEGER  :: nterma, ntermf, ntermg, ntermh
    REAL(wp) :: asln, asval, fval, gval, hval, lnxmin, t, v, x, xlow, xlow1
    CHARACTER (LEN=33) :: errmsg = 'FUNCTION CALLED WITH ARGUMENT < 0'
    CHARACTER (LEN=6)  :: fnname = 'ABRAM1'

    REAL(wp), PARAMETER :: ab1f(0:9) = [ &
      1.47285192577978807369_wp, 0.10903497570168956257_wp, &
     -0.12430675360056569753_wp, 0.306197946853493315e-2_wp, &
     -0.2218410323076511e-4_wp, 0.6989978834451e-7_wp, -0.11597076444e-9_wp, &
      0.11389776e-12_wp, -0.7173e-16_wp, 0.3e-19_wp]
    REAL(wp), PARAMETER :: ab1g(0:8) = [ &
      0.39791277949054503528_wp, -0.29045285226454720849_wp,  &
      0.1048784695465363504e-1_wp, -0.10249869522691336e-3_wp, &
      0.41150279399110e-6_wp, -0.83652638940e-9_wp, 0.97862595e-12_wp, &
     -0.71868e-15_wp, 0.35e-18_wp]
    REAL(wp), PARAMETER :: ab1h(0:8) = [ &
      0.84150292152274947030_wp, -0.7790050698774143395e-1_wp, &
      0.133992455878390993e-2_wp, -0.808503907152788e-5_wp, 0.2261858281728e-7_wp, &
     -0.3441395838e-10_wp, 0.3159858e-13_wp, -0.1884e-16_wp, 0.1e-19_wp]
    REAL(wp), PARAMETER :: ab1as(0:27) = [ &
      2.13013643429065549448_wp, 0.6371526795218539933e-1_wp, &
     -0.129334917477510647e-2_wp, 0.5678328753228265e-4_wp, -0.279434939177646e-5_wp, &
      0.5600214736787e-7_wp, 0.2392009242798e-7_wp, -0.750984865009e-8_wp, &
      0.173015330776e-8_wp, -0.36648877955e-9_wp, 0.7520758307e-10_wp, &
     -0.1517990208e-10_wp, 0.301713710e-11_wp, -0.58596718e-12_wp, &
      0.10914455e-12_wp, -0.1870536e-13_wp, 0.262542e-14_wp, -0.14627e-15_wp, &
     -0.9500e-16_wp, 0.5873e-16_wp, -0.2420e-16_wp, 0.868e-17_wp, -0.290e-17_wp, &
      0.93e-18_wp, -0.29e-18_wp, 0.9e-19_wp, -0.3e-19_wp, 0.1e-19_wp]
    REAL(wp), PARAMETER :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp, &
                           two = 2.0_wp, three = 3.0_wp, six = 6.0_wp, &
                           onehun = 100.0_wp, &
                           rt3bpi = 0.97720502380583984317_wp, &
                           onerpi = 0.56418958354775628695_wp

    ! Start calculation
    x = xvalue

    ! Error test
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    t = 2*EPSILON(zero) / onehun
    IF (x <= two) THEN
      DO  ntermf = 9, 0, -1
        IF (ABS(ab1f(ntermf)) > t) EXIT
      END DO
      DO  ntermg = 8, 0, -1
        IF (ABS(ab1g(ntermg)) > t) EXIT
      END DO
      DO  ntermh = 8, 0, -1
        IF (ABS(ab1h(ntermh)) > t) EXIT
      END DO
      t = EPSILON(zero)
      xlow1 = SQRT(two*t)
      xlow = t / two
    ELSE
      DO  nterma = 27, 0, -1
        IF (ABS(ab1as(nterma)) > t) EXIT
      END DO
      lnxmin = LOG(TINY(zero))
    END IF

    ! Code for 0 <= XVALUE <= 2
    IF (x <= two) THEN
      IF (x == zero) THEN
        fn_val = half
        RETURN
      END IF
      IF (x < xlow1) THEN
        IF (x < xlow) THEN
          fn_val = half
        ELSE
          fn_val = (one-x/onerpi-x*x*LOG(x)) * half
        END IF
        RETURN
      ELSE
        t = (x*x/two-half) - half
        fval = cheval(ntermf,ab1f,t)
        gval = cheval(ntermg,ab1g,t)
        hval = cheval(ntermh,ab1h,t)
        fn_val = fval - x * (gval/onerpi+x*LOG(x)*hval)
        RETURN
      END IF
    ELSE
      ! Code for XVALUE > 2
      v = three * ((x/two)**(two/three))
      t = (six/v-half) - half
      asval = cheval(nterma,ab1as,t)
      asln = LOG(asval*SQRT(v/three)/rt3bpi) - v
      IF (asln < lnxmin) THEN
        fn_val = zero
      ELSE
        fn_val = EXP(asln)
      END IF
    END IF
    RETURN
  END FUNCTION abram1


  FUNCTION abram2(xvalue) RESULT(fn_val)
    !! CALGO 757 Abramowitz function \(f_2(x) = 
    !! \int_{0}^{\infty} t^2 e^{-t^2 -\frac{x}{t}} \, dt\)
    !
    !   DESCRIPTION:
    !      This function calculates the Abramowitz function of order 2, defined as
    !
    !       ABRAM2(x) = integral{ 0 to infinity } (t**2) * exp( -t*t - x/t ) dt
    !
    !      The code uses Chebyshev expansions with the coefficients
    !      given to an accuracy of 20 decimal places.
    !
    !   ERROR RETURNS:
    !      If XVALUE < 0.0, the function prints a message and returns the value 0.0.
    !
    !   MACHINE-DEPENDENT CONSTANTS:
    !      NTERMF - INTEGER - No. of terms needed for the array AB2F.
    !               Recommended value such that
    !                     ABS( AB2F(NTERMF) ) < EPS/100
    !      NTERMG - INTEGER - No. of terms needed for array AB2G.
    !               Recommended value such that
    !                     ABS( AB2G(NTERMG) ) < EPS/100
    !      NTERMH - INTEGER - No. of terms needed for array AB2H.
    !               Recommended value such that
    !                     ABS( AB2H(NTERMH) ) < EPS/100
    !      NTERMA - INTEGER - No. of terms needed for array AB2AS.
    !               Recommended value such that
    !                     ABS( AB2AS(NTERMA) ) < EPS/100
    !      XLOW - REAL(wp) - The value below which
    !               ABRAM2 = root(pi)/4 to machine precision.
    !             The recommended value is EPSNEG
    !      XLOW1 - REAL(wp) - The value below which
    !                ABRAM2 = root(pi)/4 - x/2 + x**3ln(x)/6
    !              Recommended value is SQRT(2*EPSNEG)
    !      LNXMIN - REAL(wp) - The value of ln XMIN. Used to prevent
    !               exponential underflow for large X.
    !
    !     For values of EPS, EPSNEG, XMIN refer to the file MACHCON.TXT
    !
    !     The machine-dependent constants are computed internally by
    !     using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !       LOG, EXP
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !       CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN) :: xvalue
    REAL(wp)             :: fn_val

    INTEGER  :: nterma, ntermf, ntermg, ntermh
    REAL(wp) :: asln, asval, fval, gval, hval, lnxmin, t, v, x, xlow, xlow1
    CHARACTER (LEN=33) :: errmsg = 'FUNCTION CALLED WITH ARGUMENT < 0'
    CHARACTER (LEN= 6) :: fnname = 'ABRAM2'

    REAL(wp), PARAMETER :: ab2f(0:9) = [ &
      1.03612162804243713846_wp, 0.19371246626794570012_wp, &
     -0.7258758839233007378e-1_wp, 0.174790590864327399e-2_wp, &
     -0.1281223233756549e-4_wp, 0.4115018153651e-7_wp, -0.6971047256e-10_wp, &
      0.6990183e-13_wp, -0.4492e-16_wp, 0.2e-19_wp]
    REAL(wp), PARAMETER :: ab2g(0:8) = [ &
      1.46290157198630741150_wp, 0.20189466883154014317_wp, &
     -0.2908292087997129022e-1_wp, 0.47061049035270050e-3_wp, &
     -0.257922080359333e-5_wp, 0.656133712946e-8_wp, -0.914110203e-11_wp, &
      0.774276e-14_wp, -0.429e-17_wp]
    REAL(wp), PARAMETER :: ab2h(0:7) = [ &
      0.30117225010910488881_wp, -0.1588667818317623783e-1_wp, &
      0.19295936935584526e-3_wp, -0.90199587849300e-6_wp, 0.206105041837e-8_wp, &
     -0.265111806e-11_wp, 0.210864e-14_wp, -0.111e-17_wp]
    REAL(wp), PARAMETER :: ab2as(0:26) = [ &
      2.46492325304334856893_wp, 0.23142797422248905432_wp, &
     -0.94068173010085773e-3_wp, 0.8290270038089733e-4_wp, -0.883894704245866e-5_wp, &
      0.106638543567985e-5_wp, -0.13991128538529e-6_wp, 0.1939793208445e-7_wp, &
     -0.277049938375e-8_wp, 0.39590687186e-9_wp, -0.5408354342e-10_wp, &
      0.635546076e-11_wp, -0.38461613e-12_wp, -0.11696067e-12_wp, 0.6896671e-13_wp, &
     -0.2503113e-13_wp, 0.785586e-14_wp, -0.230334e-14_wp, 0.64914e-15_wp, &
     -0.17797e-15_wp, 0.4766e-16_wp, -0.1246e-16_wp, 0.316e-17_wp, -0.77e-18_wp, &
      0.18e-18_wp, -0.4e-19_wp, 0.1e-19_wp]
    REAL(wp), PARAMETER :: zero = 0.0_wp, half = 0.5_wp, two = 2.0_wp, &
                           three = 3.0_wp, six = 6.0_wp, onehun = 100.0_wp, &
                           rt3bpi = 0.97720502380583984317_wp, &
                           rtpib4 = 0.44311346272637900682_wp, &
                           onerpi = 0.56418958354775628695_wp

    ! Start calculation
    x = xvalue

    ! Error test
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    t = 2*EPSILON(zero) / onehun
    IF (x <= two) THEN
      DO  ntermf = 9, 0, -1
        IF (ABS(ab2f(ntermf)) > t) EXIT
      END DO
      DO  ntermg = 8, 0, -1
        IF (ABS(ab2g(ntermg)) > t) EXIT
      END DO
      DO  ntermh = 7, 0, -1
        IF (ABS(ab2h(ntermh)) > t) EXIT
      END DO
      xlow = EPSILON(zero)
      xlow1 = SQRT(two*xlow)
    ELSE
      DO  nterma = 26, 0, -1
        IF (ABS(ab2as(nterma)) > t) EXIT
      END DO
      lnxmin = LOG(TINY(zero))
    END IF

    ! Code for 0 <= XVALUE <= 2
    IF (x <= two) THEN
      IF (x == zero) THEN
        fn_val = rtpib4
        RETURN
      END IF
      IF (x < xlow1) THEN
        IF (x < xlow) THEN
          fn_val = rtpib4
        ELSE
          fn_val = rtpib4 - half * x + x * x * x * LOG(x) / six
        END IF
        RETURN
      ELSE
        t = (x*x/two-half) - half
        fval = cheval(ntermf,ab2f,t)
        gval = cheval(ntermg,ab2g,t)
        hval = cheval(ntermh,ab2h,t)
        fn_val = fval / onerpi + x * (x*x*LOG(x)*hval-gval)
        RETURN
      END IF
    ELSE
      ! Code for XVALUE > 2
      v = three * ((x/two)**(two/three))
      t = (six/v-half) - half
      asval = cheval(nterma,ab2as,t)
      asln = LOG(asval/rt3bpi) + LOG(v/three) - v
      IF (asln < lnxmin) THEN
        fn_val = zero
      ELSE
        fn_val = EXP(asln)
      END IF
    END IF
    RETURN
  END FUNCTION abram2


  FUNCTION airygi(xvalue) RESULT(fn_val)
    !! CALGO 757 Modified Airy function \(\mathrm{Gi}(x)\) 
    !
    !   DESCRIPTION:
    !      This subroutine computes the modified Airy function Gi(x), defined as
    !
    !        AIRYGI(x) = [ Integral{0 to infinity} sin(x*t+t^3/3) dt ] / pi
    !
    !      The approximation uses Chebyshev expansions with the coefficients
    !      given to 20 decimal places.
    !
    !   ERROR RETURNS:
    !      If x < -XHIGH1*XHIGH1 (see below for definition of XHIGH1), then
    !      the trig. functions needed for the asymptotic expansion of Bi(x)
    !      cannot be computed to any accuracy. An error message is printed
    !      and the code returns the value 0.0.
    !
    !   MACHINE-DEPENDENT CONSTANTS:
    !      NTERM1 - INTEGER - The no. of terms to be used from the array
    !                         ARGIP1. The recommended value is such that
    !                                ABS(ARGIP1(NTERM1)) < EPS/100
    !                         subject to 1 <= NTERM1 <= 30.
    !      NTERM2 - INTEGER - The no. of terms to be used from the array
    !                         ARGIP2. The recommended value is such that
    !                                ABS(ARGIP2(NTERM2)) < EPS/100
    !                         subject to 1 <= NTERM2 <= 29.
    !      NTERM3 - INTEGER - The no. of terms to be used from the array
    !                         ARGIN1. The recommended value is such that
    !                                ABS(ARGIN1(NTERM3)) < EPS/100
    !                         subject to 1 <= NTERM3 <= 42.
    !      NTERM4 - INTEGER - The no. of terms to be used from the array
    !                         ARBIN1. The recommended value is such that
    !                                ABS(ARBIN1(NTERM4)) < EPS/100
    !                         subject to 1 <= NTERM4 <= 10.
    !      NTERM5 - INTEGER - The no. of terms to be used from the array
    !                         ARBIN2. The recommended value is such that
    !                                ABS(ARBIN2(NTERM5)) < EPS/100
    !                         subject to 1 <= NTERM5 <= 11.
    !      NTERM6 - INTEGER - The no. of terms to be used from the array
    !                         ARGH2. The recommended value is such that
    !                                ABS(ARHIN1(NTERM6)) < EPS/100
    !                         subject to 1 <= NTERM6 <= 15.
    !      XLOW1 - REAL(wp) - The value such that, if -XLOW1 < x < XLOW1,
    !                     then AIRYGI = Gi(0) to machine precision.
    !                     The recommended value is   EPS.
    !      XHIGH1 - REAL(wp) - The value such that, if x > XHIGH1, then
    !                      AIRYGI = 1/(Pi*x) to machine precision.
    !                      Also used for error test - see above.
    !                      The recommended value is
    !                          cube root( 2/EPS ).
    !      XHIGH2 - REAL(wp) - The value above which AIRYGI = 0.0.
    !                      The recommended value is
    !                          1/(Pi*XMIN).
    !      XHIGH3 - REAL(wp) - The value such that, if x < XHIGH3,
    !                      then the Chebyshev expansions for the
    !                      asymptotic form of Bi(x) are not needed.
    !                      The recommended value is
    !                          -8 * cube root( 2/EPSNEG ).
    !
    !      For values of EPS, EPSNEG, and XMIN refer to the file
    !      MACHCON.TXT.
    !
    !     The machine-dependent constants are computed internally by
    !     using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !                             COS , SIN , SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: nterm1, nterm2, nterm3, nterm4, nterm5, nterm6
    REAL(wp)  :: arg, bi, cheb1, cheb2, cosz, sinz, t, temp, x, xcube,  &
                  xhigh1, xhigh2, xhigh3, xlow1, xminus, z, zeta
    CHARACTER (LEN=46) :: errmsg = 'ARGUMENT TOO NEGATIVE FOR ACCURATE COMPUTATION'
    CHARACTER (LEN= 6) :: fnname = 'AIRYGI'

    REAL(wp), PARAMETER  :: argip1(0:30) = [ &
      0.26585770795022745082_wp, -0.10500333097501922907_wp, &
      0.841347475328454492e-2_wp, 0.2021067387813439541e-1_wp, &
     -0.1559576113863552234e-1_wp, 0.564342939043256481e-2_wp, &
     -0.59776844826655809e-3_wp, -0.42833850264867728e-3_wp, &
      0.22605662380909027e-3_wp, -0.3608332945592260e-4_wp, &
     -0.785518988788901e-5_wp, 0.473252480746370e-5_wp, -0.59743513977694e-6_wp, &
     -0.15917609165602e-6_wp, 0.6336129065570e-7_wp, -0.276090232648e-8_wp, &
     -0.256064154085e-8_wp, 0.47798676856e-9_wp, 0.4488131863e-10_wp, &
     -0.2346508882e-10_wp, 0.76839085e-12_wp, 0.73227985e-12_wp, -0.8513687e-13_wp, &
     -0.1630201e-13_wp, 0.356769e-14_wp, 0.25001e-15_wp, -0.10859e-15_wp, &
     -0.158e-17_wp, 0.275e-17_wp, -0.5e-19_wp, -0.6e-19_wp]
    REAL(wp), PARAMETER  :: argip2(0:29) = [ &
      2.00473712275801486391_wp, 0.294184139364406724e-2_wp, &
      0.71369249006340167e-3_wp, 0.17526563430502267e-3_wp, 0.4359182094029882e-4_wp, &
      0.1092626947604307e-4_wp, 0.272382418399029e-5_wp, 0.66230900947687e-6_wp, &
      0.15425323370315e-6_wp, 0.3418465242306e-7_wp, 0.728157724894e-8_wp, &
      0.151588525452e-8_wp, 0.30940048039e-9_wp, 0.6149672614e-10_wp, &
      0.1202877045e-10_wp, 0.233690586e-11_wp, 0.43778068e-12_wp, 0.7996447e-13_wp, &
      0.1494075e-13_wp, 0.246790e-14_wp, 0.37672e-15_wp, 0.7701e-16_wp, &
      0.354e-17_wp, -0.49e-18_wp, 0.62e-18_wp, -0.40e-18_wp, -0.1e-19_wp, 0.2e-19_wp, &
     -0.3e-19_wp, 0.1e-19_wp]
    REAL(wp), PARAMETER  :: argin1(0:42) = [ &
     -0.20118965056732089130_wp, -0.7244175303324530499e-1_wp, &
      0.4505018923894780120e-1_wp, -0.24221371122078791099_wp, &
      0.2717884964361678294e-1_wp, -0.5729321004818179697e-1_wp, &
     -0.18382107860337763587_wp, 0.7751546082149475511e-1_wp, &
      0.18386564733927560387_wp, 0.2921504250185567173e-1_wp, &
     -0.6142294846788018811e-1_wp, -0.2999312505794616238e-1_wp, &
      0.585937118327706636e-2_wp, 0.822221658497402529e-2_wp, &
      0.132579817166846893e-2_wp, -0.96248310766565126e-3_wp, &
     -0.45065515998211807e-3_wp, 0.772423474325474e-5_wp, &
      0.5481874134758052e-4_wp, 0.1245898039742876e-4_wp, &
     -0.246196891092083e-5_wp, -0.169154183545285e-5_wp, &
     -0.16769153169442e-6_wp, 0.9636509337672e-7_wp, &
      0.3253314928030e-7_wp, 0.5091804231e-10_wp, -0.209180453553e-8_wp, &
     -0.41237387870e-9_wp, 0.4163338253e-10_wp, 0.3032532117e-10_wp, &
      0.340580529e-11_wp, -0.88444592e-12_wp, -0.31639612e-12_wp, &
     -0.1505076e-13_wp, 0.1104148e-13_wp, 0.246508e-14_wp, -0.3107e-16_wp, &
     -0.9851e-16_wp, -0.1453e-16_wp, 0.118e-17_wp, &
      0.67e-18_wp, 0.6e-19_wp, -0.1e-19_wp]
    REAL(wp), PARAMETER  :: arbin1(0:10) = [ &
      1.99983763583586155980_wp, -0.8104660923669418e-4_wp, &
      0.13475665984689e-6_wp, -0.70855847143e-9_wp, 0.748184187e-11_wp, &
     -0.12902774e-12_wp, 0.322504e-14_wp, -0.10809e-15_wp, 0.460e-17_wp, &
     -0.24e-18_wp, 0.1e-19_wp]
    REAL(wp), PARAMETER  :: arbin2(0:11) = [ &
      0.13872356453879120276_wp, -0.8239286225558228e-4_wp, &
      0.26720919509866e-6_wp, -0.207423685368e-8_wp, 0.2873392593e-10_wp, &
     -0.60873521e-12_wp, 0.1792489e-13_wp, -0.68760e-15_wp, 0.3280e-16_wp, &
     -0.188e-17_wp, 0.13e-18_wp, -0.1e-19_wp]
    REAL(wp), PARAMETER  :: arhin1(0:15) = [ &
      1.99647720399779650525_wp, -0.187563779407173213e-2_wp, &
     -0.12186470897787339e-3_wp, -0.814021609659287e-5_wp, -0.55050925953537e-6_wp, &
     -0.3763008043303e-7_wp, -0.258858362365e-8_wp, -0.17931829265e-9_wp, &
     -0.1245916873e-10_wp, -0.87171247e-12_wp, -0.6084943e-13_wp, -0.431178e-14_wp, &
     -0.29787e-15_wp, -0.2210e-16_wp, -0.136e-17_wp, -0.14e-18_wp]

    REAL(wp), PARAMETER  :: zero = 0.0_wp, one = 1.0_wp, three = 3.0_wp,  &
                             four = 4.0_wp, five = 5.0_wp, seven = 7.0_wp, &
                             minate = -8.0_wp, nine = 9.0_wp, twent8 = 28.0_wp, &
                             seven2 = 72.0_wp, onehun = 100.0_wp, one76 = 176.0_wp, &
                             five14 = 514.0_wp, one024 = 1024.0_wp,  &
                             twelhu = 1200.0_wp,  &
                             gizero = 0.20497554248200024505_wp,  &
                             onebpi = 0.31830988618379067154_wp,  &
                             piby4 = 0.78539816339744830962_wp,   &
                             rtpiin = 0.56418958354775628695_wp

    ! Start computation
    x = xvalue

    ! Compute the machine-dependent constants.
    z = EPSILON(zero)
    xlow1 = z
    arg = 2*EPSILON(zero)
    xhigh1 = one / arg
    xhigh1 = (xhigh1+xhigh1) ** (one/three)

    ! Error test
    IF (x < -xhigh1*xhigh1) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! continue with machine-dependent constants
    t = arg / onehun
    IF (x >= zero) THEN
      DO  nterm1 = 30, 0, -1
        IF (ABS(argip1(nterm1)) > t) EXIT
      END DO
      DO  nterm2 = 29, 0, -1
        IF (ABS(argip2(nterm2)) > t) EXIT
      END DO
      temp = four * piby4
      xhigh2 = one / (temp*TINY(zero))
    ELSE
      DO  nterm3 = 42, 0, -1
        IF (ABS(argin1(nterm3)) > t) EXIT
      END DO
      DO  nterm4 = 10, 0, -1
        IF (ABS(arbin1(nterm4)) > t) EXIT
      END DO
      DO  nterm5 = 11, 0, -1
        IF (ABS(arbin2(nterm5)) > t) EXIT
      END DO
      DO  nterm6 = 15, 0, -1
        IF (ABS(arhin1(nterm6)) > t) EXIT
      END DO
      temp = one / z
      xhigh3 = minate * (temp+temp) ** (one/three)
    END IF

    ! Code for x >= 0.0
    IF (x >= zero) THEN
      IF (x <= seven) THEN
        IF (x < xlow1) THEN
          fn_val = gizero
        ELSE
          t = (nine*x-twent8) / (x+twent8)
          fn_val = cheval(nterm1,argip1,t)
        END IF
      ELSE
        IF (x > xhigh1) THEN
          IF (x > xhigh2) THEN
            fn_val = zero
          ELSE
            fn_val = onebpi / x
          END IF
        ELSE
          xcube = x * x * x
          t = (twelhu-xcube) / (five14+xcube)
          fn_val = onebpi * cheval(nterm2,argip2,t) / x
        END IF
      END IF
    ELSE
      ! Code for x < 0.0
      IF (x >= minate) THEN
        IF (x > -xlow1) THEN
          fn_val = gizero
        ELSE
          t = -(x+four) / four
          fn_val = cheval(nterm3,argin1,t)
        END IF
      ELSE
        xminus = -x
        t = xminus * SQRT(xminus)
        zeta = (t+t) / three
        temp = rtpiin / SQRT(SQRT(xminus))
        cosz = COS(zeta+piby4)
        sinz = SIN(zeta+piby4) / zeta
        xcube = x * x * x
        IF (x > xhigh3) THEN
          t = -(one024/(xcube)+one)
          cheb1 = cheval(nterm4,arbin1,t)
          cheb2 = cheval(nterm5,arbin2,t)
          bi = (cosz*cheb1+sinz*cheb2) * temp
        ELSE
          bi = (cosz+sinz*five/seven2) * temp
        END IF
        t = (xcube+twelhu) / (one76-xcube)
        fn_val = bi + cheval(nterm6,arhin1,t) * onebpi / x
      END IF
    END IF
    RETURN
  END FUNCTION airygi


  FUNCTION airyhi(xvalue) RESULT(fn_val)
    !! CALGO 757 Modified Airy function \(\mathrm{Hi}(x)\) 
    !
    !   DESCRIPTION:
    !      This subroutine computes the modified Airy function Hi(x), defined as
    !
    !         AIRYHI(x) = [ Integral{0 to infinity} exp(x*t-t^3/3) dt ] / pi
    !
    !      The approximation uses Chebyshev expansions with the coefficients
    !      given to 20 decimal places.
    !
    !   ERROR RETURNS:
    !      If x > XHIGH1 (see below for definition of XHIGH1), then
    !      the asymptotic expansion of Hi(x) will cause an overflow.
    !      An error message is printed and the code returns the largest
    !      floating-pt number as the result.
    !
    !   MACHINE-DEPENDENT CONSTANTS:
    !      NTERM1 - INTEGER - The no. of terms to be used from the array
    !                         ARHIP. The recommended value is such that
    !                                ABS(ARHIP(NTERM1)) < EPS/100
    !                         subject to 1 <= NTERM1 <= 31.
    !      NTERM2 - INTEGER - The no. of terms to be used from the array
    !                         ARBIP. The recommended value is such that
    !                                ABS(ARBIP(NTERM2)) < EPS/100
    !                         subject to 1 <= NTERM2 <= 23.
    !      NTERM3 - INTEGER - The no. of terms to be used from the array
    !                         ARGIP. The recommended value is such that
    !                                ABS(ARGIP1(NTERM3)) < EPS/100
    !                         subject to 1 <= NTERM3 <= 29.
    !      NTERM4 - INTEGER - The no. of terms to be used from the array
    !                         ARHIN1. The recommended value is such that
    !                                ABS(ARHIN1(NTERM4)) < EPS/100
    !                         subject to 1 <= NTERM4 <= 21.
    !      NTERM5 - INTEGER - The no. of terms to be used from the array
    !                         ARHIN2. The recommended value is such that
    !                                ABS(ARHIN2(NTERM5)) < EPS/100
    !                         subject to 1 <= NTERM5 <= 15.
    !      XLOW1 - REAL(wp) - The value such that, if -XLOW1 < x < XLOW1,
    !                     then AIRYGI = Hi(0) to machine precision.
    !                     The recommended value is   EPS.
    !      XHIGH1 - REAL(wp) - The value such that, if x > XHIGH1, then
    !                      overflow might occur. The recommended value is
    !                      computed as follows:
    !                           compute Z = 1.5*LOG(XMAX)
    !                        XHIGH1 = ( Z + LOG(Z)/4 + LOG(PI)/2 )**(2/3)
    !      XNEG1 - REAL(wp) - The value below which AIRYHI = 0.0.
    !                     The recommended value is
    !                          -1/(Pi*XMIN).
    !      XNEG2 - REAL(wp) - The value such that, if x < XNEG2, then
    !                      AIRYHI = -1/(Pi*x) to machine precision.
    !                      The recommended value is
    !                          -cube root( 2/EPS ).
    !      XMAX - REAL(wp) - The largest possible floating-pt. number.
    !                    This is the value given to the function
    !                    if x > XHIGH1.
    !
    !      For values of EPS, EPSNEG, XMIN  and XMAX refer to the file
    !      MACHCON.TXT.
    !
    !     The machine-dependent constants are computed internally by
    !     using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !       EXP, LOG, SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !       CHEVAL, ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: nterm1, nterm2, nterm3, nterm4, nterm5
    REAL(wp)  :: bi, gi, t, temp, x, xcube, xhigh1,  &
                  xlow1, xmax, xneg1, xneg2, z, zeta
    CHARACTER (LEN=30)  :: errmsg = 'ARGUMENT TO FUNCTION TOO LARGE'
    CHARACTER (LEN= 6)  :: fnname = 'AIRYHI'

    REAL(wp), PARAMETER  :: arhip(0:31) = [ &
      1.24013562561762831114_wp, 0.64856341973926535804_wp, &
      0.55236252592114903246_wp, 0.20975122073857566794_wp, &
      0.12025669118052373568_wp, 0.3768224931095393785e-1_wp, &
      0.1651088671548071651e-1_wp, 0.455922755211570993e-2_wp, &
      0.161828480477635013e-2_wp, 0.40841282508126663e-3_wp, &
      0.12196479721394051e-3_wp, 0.2865064098657610e-4_wp, 0.742221556424344e-5_wp, &
      0.163536231932831e-5_wp, 0.37713908188749e-6_wp, 0.7815800336008e-7_wp, &
      0.1638447121370e-7_wp, 0.319857665992e-8_wp, 0.61933905307e-9_wp, &
      0.11411161191e-9_wp, 0.2064923454e-10_wp, &
      0.360018664e-11_wp, 0.61401849e-12_wp, 0.10162125e-12_wp, 0.1643701e-13_wp, &
      0.259084e-14_wp, 0.39931e-15_wp, 0.6014e-16_wp, 0.886e-17_wp, 0.128e-17_wp, &
      0.18e-18_wp, 0.3e-19_wp]
    REAL(wp), PARAMETER  :: arbip(0:23) = [ &
      2.00582138209759064905_wp, 0.294478449170441549e-2_wp,  &
      0.3489754514775355e-4_wp, 0.83389733374343e-6_wp, 0.3136215471813e-7_wp,  &
      0.167865306015e-8_wp, 0.12217934059e-9_wp, 0.1191584139e-10_wp, &
      0.154142553e-11_wp, 0.24844455e-12_wp, 0.4213012e-13_wp, 0.505293e-14_wp, &
     -0.60032e-15_wp, -0.65474e-15_wp, -0.22364e-15_wp, -0.3015e-16_wp, &
      0.959e-17_wp, 0.616e-17_wp, 0.97e-18_wp, -0.37e-18_wp, &
     -0.21e-18_wp, -0.1e-19_wp, 0.2e-19_wp, 0.1e-19_wp]
    REAL(wp), PARAMETER  :: argip1(0:29) = [ &
      2.00473712275801486391_wp, 0.294184139364406724e-2_wp,
      0.71369249006340167e-3_wp, 0.17526563430502267e-3_wp, 0.4359182094029882e-4_wp, &
      0.1092626947604307e-4_wp, 0.272382418399029e-5_wp, 0.66230900947687e-6_wp, &
      0.15425323370315e-6_wp, 0.3418465242306e-7_wp, 0.728157724894e-8_wp, &
      0.151588525452e-8_wp, 0.30940048039e-9_wp,  0.6149672614e-10_wp, &
      0.1202877045e-10_wp, 0.233690586e-11_wp, 0.43778068e-12_wp, 0.7996447e-13_wp, &
      0.1494075e-13_wp, 0.246790e-14_wp, 0.37672e-15_wp, 0.7701e-16_wp, 0.354e-17_wp, &
     -0.49e-18_wp, 0.62e-18_wp, -0.40e-18_wp, -0.1e-19_wp, 0.2e-19_wp, -0.3e-19_wp, &
      0.1e-19_wp]
    REAL(wp), PARAMETER  :: arhin1(0:21) = [ &
      0.31481017206423404116_wp, -0.16414499216588964341_wp,  &
      0.6176651597730913071e-1_wp, -0.1971881185935933028e-1_wp,   &
      0.536902830023331343e-2_wp, -0.124977068439663038e-2_wp,
      0.24835515596994933e-3_wp, -0.4187024096746630e-4_wp, 0.590945437979124e-5_wp, &
     -0.68063541184345e-6_wp, 0.6072897629164e-7_wp, -0.367130349242e-8_wp,
      0.7078017552e-10_wp, 0.1187894334e-10_wp, -0.120898723e-11_wp, &
      0.1189656e-13_wp, 0.594128e-14_wp, -0.32257e-15_wp, -0.2290e-16_wp, &
      0.253e-17_wp, 0.9e-19_wp, -0.2e-19_wp]
    REAL(wp), PARAMETER  :: arhin2(0:15) = [ &
      1.99647720399779650525_wp, -0.187563779407173213e-2_wp,  &
     -0.12186470897787339e-3_wp, -0.814021609659287e-5_wp, -0.55050925953537e-6_wp,  &
     -0.3763008043303e-7_wp, -0.258858362365e-8_wp, -0.17931829265e-9_wp,  &
     -0.1245916873e-10_wp, -0.87171247e-12_wp, -0.6084943e-13_wp, -0.431178e-14_wp,  &
     -0.29787e-15_wp, -0.2210e-16_wp, -0.136e-17_wp, -0.14e-18_wp]
    REAL(wp), PARAMETER  :: zero = 0.0_wp, one = 1.0_wp, two = 2.0_wp,  &
                             three = 3.0_wp, four = 4.0_wp, seven = 7.0_wp,  &
                             minate = -8.0_wp, twelve = 12.0_wp, one76 = 176.0_wp, &
                             thre43 = 343.0_wp, five14 = 514.0_wp,  &
                             twelhu = 1200.0_wp, onehun = 100.0_wp,  &
                             hizero = 0.40995108496400049010_wp,  &
                             lnrtpi = 0.57236494292470008707_wp,  &
                             onebpi = 0.31830988618379067154_wp

    ! Start computation
    x = xvalue

    ! Compute the machine-dependent constants.
    xmax = HUGE(zero)
    temp = three * LOG(xmax) / two
    zeta = (temp + LOG(temp)/four - LOG(onebpi)/two)
    xhigh1 = zeta ** (two/three)

    ! Error test
    IF (x > xhigh1) THEN
      CALL errprn(fnname,errmsg)
      fn_val = xmax
      RETURN
    END IF

    ! continue with machine-dependent constants
    z = EPSILON(zero)
    xlow1 = z
    t = z / onehun
    IF (x >= zero) THEN
      DO  nterm1 = 31, 0, -1
        IF (ABS(arhip(nterm1)) > t) EXIT
      END DO
      DO  nterm2 = 23, 0, -1
        IF (ABS(arbip(nterm2)) > t) EXIT
      END DO
      DO  nterm3 = 29, 0, -1
        IF (ABS(argip1(nterm3)) > t) EXIT
      END DO
    ELSE
      DO  nterm4 = 21, 0, -1
        IF (ABS(arhin1(nterm4)) > t) EXIT
      END DO
      DO  nterm5 = 15, 0, -1
        IF (ABS(arhin2(nterm5)) > t) EXIT
      END DO
      temp = one / onebpi
      xneg1 = -one / (temp*TINY(zero))
      xneg2 = -((two/z)**(one/three))
    END IF

    ! Code for x >= 0.0
    IF (x >= zero) THEN
      IF (x <= seven) THEN
        IF (x < xlow1) THEN
          fn_val = hizero
        ELSE
          t = (x+x) / seven - one
          temp = (x+x+x) / two
          fn_val = EXP(temp) * cheval(nterm1,arhip,t)
        END IF
      ELSE
        xcube = x * x * x
        temp = SQRT(xcube)
        zeta = (temp+temp) / three
        t = two * (SQRT(thre43/xcube)) - one
        temp = cheval(nterm2,arbip,t)
        temp = zeta + LOG(temp) - LOG(x) / four - lnrtpi
        bi = EXP(temp)
        t = (twelhu-xcube) / (xcube+five14)
        gi = cheval(nterm3,argip1,t) * onebpi / x
        fn_val = bi - gi
      END IF
    ELSE
      ! Code for x < 0.0
      IF (x >= minate) THEN
        IF (x > -xlow1) THEN
          fn_val = hizero
        ELSE
          t = (four*x+twelve) / (x-twelve)
          fn_val = cheval(nterm4,arhin1,t)
        END IF
      ELSE
        IF (x < xneg1) THEN
          fn_val = zero
        ELSE
          IF (x < xneg2) THEN
            temp = one
          ELSE
            xcube = x * x * x
            t = (xcube+twelhu) / (one76-xcube)
            temp = cheval(nterm5,arhin2,t)
          END IF
          fn_val = -temp * onebpi / x
        END IF
      END IF
    END IF
    RETURN
  END FUNCTION airyhi


  FUNCTION airint(xvalue) RESULT(fn_val)
    !! CALGO 757 Integral of the Airy function \(\mathrm{Ai}(x)\) from \(0\) to \(x\) 
    ! 
    !   DESCRIPTION:
    !      This function calculates the integral of the Airy function Ai, defined as
    !
    !         AIRINT(x) = {integral 0 to x} Ai(t) dt
    !
    !      The program uses Chebyshev expansions, the coefficients of which
    !      are given to 20 decimal places.
    !
    !   ERROR RETURNS:
    !      If the argument is too large and negative, it is impossible
    !      to accurately compute the necessary SIN and COS functions.
    !      An error message is printed, and the program returns the
    !      value -2/3 (the value at -infinity).
    !
    !   MACHINE-DEPENDENT CONSTANTS:
    !      NTERM1 - INTEGER - The no. of terms to be used from the array
    !                          AAINT1. The recommended value is such that
    !                             ABS(AAINT1(NTERM1)) < EPS/100,
    !                          subject to 1 <= NTERM1 <= 25.
    !      NTERM2 - INTEGER - The no. of terms to be used from the array
    !                          AAINT2. The recommended value is such that
    !                             ABS(AAINT2(NTERM2)) < EPS/100,
    !                          subject to 1 <= NTERM2 <= 21.
    !      NTERM3 - INTEGER - The no. of terms to be used from the array
    !                          AAINT3. The recommended value is such that
    !                             ABS(AAINT3(NTERM3)) < EPS/100,
    !                          subject to 1 <= NTERM3 <= 40.
    !      NTERM4 - INTEGER - The no. of terms to be used from the array
    !                          AAINT4. The recommended value is such that
    !                             ABS(AAINT4(NTERM4)) < EPS/100,
    !                          subject to 1 <= NTERM4 <= 17.
    !      NTERM5 - INTEGER - The no. of terms to be used from the array
    !                          AAINT5. The recommended value is such that
    !                             ABS(AAINT5(NTERM5)) < EPS/100,
    !                          subject to 1 <= NTERM5 <= 17.
    !      XLOW1 - REAL(wp) - The value such that, if |x| < XLOW1,
    !                          AIRINT(x) = x * Ai(0)
    !                     to machine precision. The recommended value is
    !                          2 * EPSNEG.
    !      XHIGH1 - REAL(wp) - The value such that, if x > XHIGH1,
    !                          AIRINT(x) = 1/3,
    !                      to machine precision. The recommended value is
    !                          (-1.5*LOG(EPSNEG)) ** (2/3).
    !      XNEG1 - REAL(wp) - The value such that, if x < XNEG1,
    !                     the trigonometric functions in the asymptotic
    !                     expansion cannot be calculated accurately.
    !                     The recommended value is
    !                          -(1/((EPS)**2/3))
    !
    !      For values of EPS and EPSNEG, refer to the file MACHCON.TXT.
    !
    !     The machine-dependent constants are computed internally by
    !     using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !       COS, EXP, SIN, SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !       CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: nterm1, nterm2, nterm3, nterm4, nterm5
    REAL(wp)  :: arg, gval, hval, t, temp, x, xhigh1, xlow1, xneg1, z
    CHARACTER (LEN=46)  :: errmsg = 'AIRINT'
    CHARACTER (LEN= 6)  :: fnname = 'FUNCTION TOO NEGATIVE FOR ACCURATE COMPUTATION'

    REAL(wp), PARAMETER  :: aaint1(0:25) = [ &
      0.37713517694683695526_wp, -0.13318868432407947431_wp,  &
      0.3152497374782884809e-1_wp, -0.318543076436574077e-2_wp,  &
     -0.87398764698621915e-3_wp, 0.46699497655396971e-3_wp, &
     -0.9544936738983692e-4_wp, 0.542705687156716e-5_wp, 0.239496406252188e-5_wp, &
     -0.75690270205649e-6_wp, 0.9050138584518e-7_wp, 0.320529456043e-8_wp, &
     -0.303825536444e-8_wp, 0.48900118596e-9_wp, -0.1839820572e-10_wp, &
     -0.711247519e-11_wp, 0.151774419e-11_wp, -0.10801922e-12_wp, -0.963542e-14_wp, &
      0.313425e-14_wp, -0.29446e-15_wp, -0.477e-17_wp, &
      0.461e-17_wp, -0.53e-18_wp, 0.1e-19_wp, 0.1e-19_wp]
    REAL(wp), PARAMETER  :: aaint2(0:21) = [ &
      1.92002524081984009769_wp, -0.4220049417256287021e-1_wp, &
     -0.239457722965939223e-2_wp, -0.19564070483352971e-3_wp, &
     -0.1547252891056112e-4_wp, -0.140490186137889e-5_wp, -0.12128014271367e-6_wp, &
     -0.1179186050192e-7_wp, -0.104315578788e-8_wp, -0.10908209293e-9_wp, &
     -0.929633045e-11_wp, -0.110946520e-11_wp, -0.7816483e-13_wp, -0.1319661e-13_wp, &
     -0.36823e-15_wp, -0.21505e-15_wp, 0.1238e-16_wp, -0.557e-17_wp, 0.84e-18_wp, &
     -0.21e-18_wp, 0.4e-19_wp, -0.1e-19_wp]
    REAL(wp), PARAMETER  :: aaint3(0:40) = [ &
        0.47985893264791052053_wp, -0.19272375126169608863_wp, &
        0.2051154129525428189e-1_wp, 0.6332000070732488786e-1_wp, &
       -0.5093322261845754082e-1_wp, 0.1284424078661663016e-1_wp, &
        0.2760137088989479413e-1_wp, -0.1547066673866649507e-1_wp, &
       -0.1496864655389316026e-1_wp, 0.336617614173574541e-2_wp, &
        0.530851163518892985e-2_wp, 0.41371226458555081e-3_wp, &
       -0.102490579926726266e-2_wp, -0.32508221672025853e-3_wp, &
        0.8608660957169213e-4_wp, 0.6671367298120775e-4_wp, 0.449205999318095e-5_wp,  &
       -0.670427230958249e-5_wp, -0.196636570085009e-5_wp, 0.22229677407226e-6_wp,  &
        0.22332222949137e-6_wp, 0.2803313766457e-7_wp, -0.1155651663619e-7_wp,  &
       -0.433069821736e-8_wp, -0.6227777938e-10_wp, 0.26432664903e-9_wp, 0.5333881114e-10_wp, &
       -0.522957269e-11_wp, -0.382229283e-11_wp, -0.40958233e-12_wp, 0.11515622e-12_wp,  &
        0.3875766e-13_wp, 0.140283e-14_wp, -0.141526e-14_wp, -0.28746e-15_wp, 0.923e-17_wp,  &
        0.1224e-16_wp, 0.157e-17_wp, -0.19e-18_wp, -0.8e-19_wp, -0.1e-19_wp]

    REAL(wp), PARAMETER  :: aaint4(0:17) = [ &
        1.99653305828522730048_wp, -0.187541177605417759e-2_wp,  &
       -0.15377536280305750e-3_wp, -0.1283112967682349e-4_wp, -0.108128481964162e-5_wp, &
       -0.9182131174057e-7_wp, -0.784160590960e-8_wp, -0.67292453878e-9_wp,  &
       -0.5796325198e-10_wp, -0.501040991e-11_wp, -0.43420222e-12_wp, -0.3774305e-13_wp,  &
       -0.328473e-14_wp, -0.28700e-15_wp, -0.2502e-16_wp, -0.220e-17_wp, -0.19e-18_wp, -0.2e-19_wp]
    REAL(wp), PARAMETER  :: aaint5(0:17) = [ &
        1.13024602034465716133_wp, -0.464718064639872334e-2_wp,  &
       -0.35137413382693203e-3_wp, -0.2768117872545185e-4_wp, -0.222057452558107e-5_wp,  &
       -0.18089142365974e-6_wp, -0.1487613383373e-7_wp, -0.123515388168e-8_wp,  &
       -0.10310104257e-9_wp, -0.867493013e-11_wp, -0.73080054e-12_wp, -0.6223561e-13_wp,   &
       -0.525128e-14_wp, -0.45677e-15_wp, -0.3748e-16_wp, -0.356e-17_wp, -0.23e-18_wp, -0.4e-19_wp]

    REAL(wp), PARAMETER  :: zero = 0.0_wp, one = 1.0_wp, two = 2.0_wp,  &
                             three = 3.0_wp, four = 4.0_wp, eight = 8.0_wp,  &
                             nine = 9.0_wp, forty1 = 41.0_wp, onehun = 100.0_wp, &
                             ninhun = 900.0_wp, fr996 = 4996.0_wp,  &
                             piby4 = 0.78539816339744830962_wp,  &
                             pitim6 = 18.84955592153875943078_wp,  &
                             rt2b3p = 0.46065886596178063902_wp,  &
                             airzer = 0.35502805388781723926_wp

    ! Start computation
    x = xvalue

    ! Compute the machine-dependent constants.
    z = EPSILON(zero)
    xlow1 = two * z
    arg = 2*EPSILON(zero)
    xneg1 = -one / (arg**(two/three))

    ! Error test
    IF (x < xneg1) THEN
      CALL errprn(fnname,errmsg)
      fn_val = -two / three
      RETURN
    END IF

    ! continue with machine-dependent constants
    t = arg / onehun
    IF (x >= zero) THEN
      DO  nterm1 = 25, 0, -1
        IF (ABS(aaint1(nterm1)) > t) EXIT
      END DO
      DO  nterm2 = 21, 0, -1
        IF (ABS(aaint2(nterm2)) > t) EXIT
      END DO
      xhigh1 = (-three*LOG(z)/two) ** (two/three)
    ELSE
      DO  nterm3 = 40, 0, -1
        IF (ABS(aaint3(nterm3)) > t) EXIT
      END DO
      DO  nterm4 = 17, 0, -1
        IF (ABS(aaint4(nterm4)) > t) EXIT
      END DO
      DO  nterm5 = 17, 0, -1
        IF (ABS(aaint5(nterm5)) > t) EXIT
      END DO
    END IF

    ! Code for x >= 0
    IF (x >= zero) THEN
      IF (x <= four) THEN
        IF (x < xlow1) THEN
          fn_val = airzer * x
        ELSE
          t = x / two - one
          fn_val = cheval(nterm1,aaint1,t) * x
        END IF
      ELSE
        IF (x > xhigh1) THEN
          temp = zero
        ELSE
          z = (x+x) * SQRT(x) / three
          temp = three * z
          t = (forty1-temp) / (nine+temp)
          temp = EXP(-z) * cheval(nterm2,aaint2,t) / SQRT(pitim6*z)
        END IF
        fn_val = one / three - temp
      END IF
    ELSE
      ! Code for x < 0
      IF (x >= -eight) THEN
        IF (x > -xlow1) THEN
          fn_val = airzer * x
        ELSE
          t = -x / four - one
          fn_val = x * cheval(nterm3,aaint3,t)
        END IF
      ELSE
        z = -(x+x) * SQRT(-x) / three
        arg = z + piby4
        temp = nine * z * z
        t = (fr996-temp) / (ninhun+temp)
        gval = cheval(nterm4,aaint4,t)
        hval = cheval(nterm5,aaint5,t)
        temp = gval * COS(arg) + hval * SIN(arg) / z
        fn_val = rt2b3p * temp / SQRT(z) - two / three
      END IF
    END IF
    RETURN
  END FUNCTION airint


  FUNCTION birint(xvalue) RESULT(fn_val)
    !! CALGO 757 Integral of the Airy function \(\mathrm{Bi}(x)\) from \(0\) to \(x\) 
    !
    !   DESCRIPTION:
    !      This function calculates the integral of the Airy function Bi, defined
    !
    !          BIRINT(x) = integral{0 to x} Bi(t) dt
    !
    !      The program uses Chebyshev expansions, the coefficients of which
    !      are given to 20 decimal places.
    !
    !   ERROR RETURNS:
    !      If the function is too large and positive the correct
    !      value would overflow. An error message is printed and the
    !      program returns the value XMAX.
    !
    !      If the argument is too large and negative, it is impossible
    !      to accurately compute the necessary SIN and COS functions,
    !      for the asymptotic expansion.
    !      An error message is printed, and the program returns the
    !      value 0 (the value at -infinity).
    !
    !   MACHINE-DEPENDENT CONSTANTS:
    !      NTERM1 - INTEGER - The no. of terms to be used from the array
    !                          ABINT1. The recommended value is such that
    !                             ABS(ABINT1(NTERM1)) < EPS/100,
    !                          subject to 1 <= NTERM1 <= 36.
    !      NTERM2 - INTEGER - The no. of terms to be used from the array
    !                          ABINT2. The recommended value is such that
    !                             ABS(ABINT2(NTERM2)) < EPS/100,
    !                          subject to 1 <= NTERM2 <= 37.
    !      NTERM3 - INTEGER - The no. of terms to be used from the array
    !                          ABINT3. The recommended value is such that
    !                             ABS(ABINT3(NTERM3)) < EPS/100,
    !                          subject to 1 <= NTERM3 <= 37.
    !      NTERM4 - INTEGER - The no. of terms to be used from the array
    !                          ABINT4. The recommended value is such that
    !                             ABS(ABINT4(NTERM4)) < EPS/100,
    !                          subject to 1 <= NTERM4 <= 20.
    !      NTERM5 - INTEGER - The no. of terms to be used from the array
    !                          ABINT5. The recommended value is such that
    !                             ABS(ABINT5(NTERM5)) < EPS/100,
    !                          subject to 1 <= NTERM5 <= 20.
    !      XLOW1 - REAL(wp) - The value such that, if |x| < XLOW1,
    !                          BIRINT(x) = x * Bi(0)
    !                     to machine precision. The recommended value is
    !                          2 * EPSNEG.
    !      XHIGH1 - REAL(wp) - The value such that, if x > XHIGH1,
    !                      the function value would overflow.
    !                      The recommended value is computed as
    !                          z = ln(XMAX) + 0.5ln(ln(XMAX)),
    !                          XHIGH1 = (3z/2)^(2/3)
    !      XNEG1 - REAL(wp) - The value such that, if x < XNEG1,
    !                     the trigonometric functions in the asymptotic
    !                     expansion cannot be calculated accurately.
    !                     The recommended value is
    !                          -(1/((EPS)**2/3))
    !      XMAX - REAL(wp) - The value of the largest positive floating-pt
    !                    number. Used in giving a value to the function
    !                    if x > XHIGH1.
    !
    !     For values of EPS, EPSNEG, and XMAX see the file MACHCON.TXT.
    !
    !     The machine-dependent constants are computed internally by
    !     using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !       COS, EXP, LOG, SIN, SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !       CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: nterm1, nterm2, nterm3, nterm4, nterm5
    REAL(wp)  :: arg, f1, f2, t, temp, x, xlow1, xhigh1, xmax, xneg1, z
    CHARACTER (LEN=31) :: ermsg2 = 'ARGUMENT TOO LARGE AND NEGATIVE'
    CHARACTER (LEN=31) :: ermsg1 = 'ARGUMENT TOO LARGE AND POSITIVE'
    CHARACTER (LEN= 6) :: fnname = 'BIRINT'

    REAL(wp), PARAMETER  :: abint1(0:36) = [ &
        0.38683352445038543350_wp, -0.8823213550888908821e-1_wp,  &
        0.21463937440355429239_wp, -0.4205347375891315126e-1_wp,  &
        0.5932422547496086771e-1_wp, -0.840787081124270210e-2_wp,  &
        0.871824772778487955e-2_wp, -0.12191600199613455e-3_wp, 0.44024821786023234e-3_wp, &
        0.27894686666386678e-3_wp, -0.7052804689785537e-4_wp, 0.5901080066770100e-4_wp,  &
       -0.1370862587982142e-4_wp, 0.505962573749073e-5_wp, -0.51598837766735e-6_wp,  &
        0.397511312349e-8_wp, 0.9524985978055e-7_wp, -0.3681435887321e-7_wp,  &
        0.1248391688136e-7_wp, -0.249097619137e-8_wp, 0.31775245551e-9_wp,  &
        0.5434365270e-10_wp, -0.4024566915e-10_wp, 0.1393855527e-10_wp, -0.303817509e-11_wp,  &
        0.40809511e-12_wp, 0.1634116e-13_wp, -0.2683809e-13_wp, 0.896641e-14_wp,  &
       -0.183089e-14_wp, 0.21333e-15_wp, 0.1108e-16_wp, -0.1276e-16_wp, 0.363e-17_wp, -0.62e-18_wp, &
        0.5e-19_wp, 0.1e-19_wp]

    REAL(wp), PARAMETER  :: abint2(0:37) = [ &
        2.04122078602516135181_wp, 0.2124133918621221230e-1_wp,  &
        0.66617599766706276e-3_wp, 0.3842047982808254e-4_wp, 0.362310366020439e-5_wp,  &
        0.50351990115074e-6_wp, 0.7961648702253e-7_wp, 0.717808442336e-8_wp,  &
       -0.267770159104e-8_wp, -0.168489514699e-8_wp, -0.36811757255e-9_wp, &
        0.4757128727e-10_wp, 0.5263621945e-10_wp, 0.778973500e-11_wp, -0.460546143e-11_wp,  &
       -0.183433736e-11_wp, 0.32191249e-12_wp, 0.29352060e-12_wp, -0.1657935e-13_wp,  &
       -0.4483808e-13_wp, 0.27907e-15_wp, 0.711921e-14_wp, -0.1042e-16_wp, -0.119591e-14_wp,  &
        0.4606e-16_wp, 0.20884e-15_wp, -0.2416e-16_wp, -0.3638e-16_wp, 0.863e-17_wp, 0.591e-17_wp,  &
       -0.256e-17_wp, -0.77e-18_wp, 0.66e-18_wp, 0.3e-19_wp, -0.15e-18_wp, 0.2e-19_wp, 0.3e-19_wp,  &
       -0.1e-19_wp]

    REAL(wp), PARAMETER  :: abint3(0:37) = [ &
        0.31076961598640349251_wp, -0.27528845887452542718_wp,  &
        0.17355965706136543928_wp, -0.5544017909492843130e-1_wp,  &
       -0.2251265478295950941e-1_wp, 0.4107347447812521894e-1_wp,  &
        0.984761275464262480e-2_wp, -0.1555618141666041932e-1_wp,  &
       -0.560871870730279234e-2_wp, 0.246017783322230475e-2_wp, 0.165740392292336978e-2_wp, &
       -0.3277587501435402e-4_wp, -0.24434680860514925e-3_wp, -0.5035305196152321e-4_wp,  &
        0.1630264722247854e-4_wp, 0.851914057780934e-5_wp, 0.29790363004664e-6_wp,  &
       -0.64389707896401e-6_wp, -0.15046988145803e-6_wp, 0.1587013535823e-7_wp,  &
        0.1276766299622e-7_wp, 0.140578534199e-8_wp, -0.46564739741e-9_wp,  &
       -0.15682748791e-9_wp, -0.403893560e-11_wp, 0.666708192e-11_wp, 0.128869380e-11_wp,  &
       -0.6968663e-13_wp, -0.6254319e-13_wp, -0.718392e-14_wp, 0.115296e-14_wp, 0.42276e-15_wp,  &
        0.2493e-16_wp, -0.971e-17_wp, -0.216e-17_wp, -0.2e-19_wp, 0.6e-19_wp, 0.1e-19_wp]

    REAL(wp), PARAMETER  :: abint4(0:20) = [ &
        1.99507959313352047614_wp, -0.273736375970692738e-2_wp,  &
       -0.30897113081285850e-3_wp, -0.3550101982798577e-4_wp, -0.412179271520133e-5_wp,  &
       -0.48235892316833e-6_wp, -0.5678730727927e-7_wp, -0.671874810365e-8_wp,  &
       -0.79811649857e-9_wp, -0.9514271478e-10_wp, -0.1137468966e-10_wp, -0.136359969e-11_wp, &
       -0.16381418e-12_wp, -0.1972575e-13_wp, -0.237844e-14_wp, -0.28752e-15_wp, -0.3475e-16_wp, &
       -0.422e-17_wp, -0.51e-18_wp, -0.6e-19_wp, -0.1e-19_wp]

    REAL(wp), PARAMETER  :: abint5(0:20) = [ &
        1.12672081961782566017_wp, -0.671405567525561198e-2_wp,  &
       -0.69812918017832969e-3_wp, -0.7561689886425276e-4_wp, -0.834985574510207e-5_wp,  &
       -0.93630298232480e-6_wp, -0.10608556296250e-6_wp, -0.1213128916741e-7_wp,  &
       -0.139631129765e-8_wp, -0.16178918054e-9_wp, -0.1882307907e-10_wp, -0.220272985e-11_wp, &
       -0.25816189e-12_wp, -0.3047964e-13_wp, -0.358370e-14_wp, -0.42831e-15_wp, -0.4993e-16_wp,  &
       -0.617e-17_wp, -0.68e-18_wp, -0.10e-18_wp, -0.1e-19_wp]

    REAL(wp), PARAMETER  :: zero = 0.0_wp, one = 1.0_wp, onept5 = 1.5_wp,  &
                             three = 3.0_wp, four = 4.0_wp, seven = 7.0_wp,  &
                             eight = 8.0_wp, nine = 9.0_wp, sixten = 16.0_wp,  &
                             onehun = 100.0_wp, ninhun = 900.0_wp,  &
                             thr644 = 3644.0_wp,  &
                             piby4 = 0.78539816339744830962_wp,  &
                             rt2b3p = 0.46065886596178063902_wp, &
                             birzer = 0.61492662744600073515_wp

    ! Start computation
    x = xvalue

    ! Compute the machine-dependent constants.
    t = EPSILON(zero)
    f2 = one + one
    xneg1 = -one / (t**(f2/three))
    xmax = HUGE(zero)
    f1 = LOG(xmax)
    temp = f1 + LOG(f1) / f2
    xhigh1 = (three*temp/f2) ** (f2/three)

    ! Error test
    IF (x > xhigh1) THEN
      CALL errprn(fnname,ermsg1)
      fn_val = xmax
      RETURN
    END IF
    IF (x < xneg1) THEN
      CALL errprn(fnname,ermsg2)
      fn_val = zero
      RETURN
    END IF

    ! continue with machine-dependent constants
    xlow1 = f2 * t
    t = t / onehun
    IF (x >= zero) THEN
      DO  nterm1 = 36, 0, -1
        IF (ABS(abint1(nterm1)) > t) EXIT
      END DO
      DO  nterm2 = 37, 0, -1
        IF (ABS(abint2(nterm2)) > t) EXIT
      END DO
    ELSE
      DO  nterm3 = 37, 0, -1
        IF (ABS(abint3(nterm3)) > t) EXIT
      END DO
      DO  nterm4 = 20, 0, -1
        IF (ABS(abint4(nterm4)) > t) EXIT
      END DO
      DO  nterm5 = 20, 0, -1
        IF (ABS(abint5(nterm5)) > t) EXIT
      END DO
    END IF

    ! Code for x >= 0.0
    IF (x >= zero) THEN
      IF (x < xlow1) THEN
        fn_val = birzer * x
      ELSE
        IF (x <= eight) THEN
          t = x / four - one
          fn_val = x * EXP(onept5*x) * cheval(nterm1,abint1,t)
        ELSE
          t = sixten * SQRT(eight/x) / x - one
          z = (x+x) * SQRT(x) / three
          temp = rt2b3p * cheval(nterm2,abint2,t) / SQRT(z)
          temp = z + LOG(temp)
          fn_val = EXP(temp)
        END IF
      END IF
    ELSE
      ! Code for x < 0.0
      IF (x >= -seven) THEN
        IF (x > -xlow1) THEN
          fn_val = birzer * x
        ELSE
          t = -(x+x) / seven - one
          fn_val = x * cheval(nterm3,abint3,t)
        END IF
      ELSE
        z = -(x+x) * SQRT(-x) / three
        arg = z + piby4
        temp = nine * z * z
        t = (thr644-temp) / (ninhun+temp)
        f1 = cheval(nterm4,abint4,t) * SIN(arg)
        f2 = cheval(nterm5,abint5,t) * COS(arg) / z
        fn_val = (f2-f1) * rt2b3p / SQRT(z)
      END IF
    END IF
    RETURN
  END FUNCTION birint


  FUNCTION j0int(xvalue) RESULT(fn_val)
    !! CALGO 757 Integral of the Bessel function \(J_0(x)\) from \(0\) to \(x\) 
    !
    !   DESCRIPTION:
    !      This function calculates the integral of the Bessel
    !      function J0, defined as
    !
    !        J0INT(x) = {integral 0 to x} J0(t) dt
    !
    !      The code uses Chebyshev expansions whose coefficients are
    !      given to 20 decimal places.
    !
    !   ERROR RETURNS:
    !      If the value of |x| is too large, it is impossible to
    !      accurately compute the trigonometric functions used. An
    !      error message is printed, and the function returns the value 1.0.
    !
    !   MACHINE-DEPENDENT CONSTANTS:
    !      NTERM1 - The no. of terms to be used from the array
    !                ARJ01. The recommended value is such that
    !                   ABS(ARJ01(NTERM1)) < EPS/100, provided that
    !      NTERM2 - The no. of terms to be used from the array
    !                ARJ0A1. The recommended value is such that
    !                   ABS(ARJ0A1(NTERM2)) < EPS/100, provided that
    !      NTERM3 - The no. of terms to be used from the array
    !                ARJ0A2. The recommended value is such that
    !                   ABS(ARJ0A2(NTERM3)) < EPS/100, provided that
    !      XLOW - The value of |x| below which J0INT(x) = x to
    !             machine-precision. The recommended value is
    !                 sqrt(12*EPSNEG)
    !      XHIGH - The value of |x| above which it is impossible
    !              to calculate (x-pi/4) accurately. The recommended
    !              value is      1/EPSNEG
    !
    !      For values of EPS and EPSNEG for various machine/compiler
    !      combinations refer to the file MACHCON.TXT.
    !
    !      The machine-dependent constants are computed internally by
    !      using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !      COS , SIN , SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: ind, nterm1, nterm2, nterm3
    REAL(wp)  :: pib41, t, temp, x, xhigh, xlow, xmpi4
    CHARACTER (LEN=26)  :: errmsg = 'ARGUMENT TOO LARGE IN SIZE'
    CHARACTER (LEN= 6)  :: fnname = 'J0INT '

    REAL(wp), PARAMETER  :: zero = 0.0_wp, one = 1.0_wp, twelve = 12.0_wp,  &
                             sixten = 16.0_wp, onehun = 100.0_wp, one28 = 128.0_wp, &
                             five12 = 512.0_wp,  &
                             rt2bpi = 0.79788456080286535588_wp,  &
                             pib411 = 201.0_wp, pib412 = 256.0_wp,  &
                             pib42 = 0.24191339744830961566e-3_wp

    REAL(wp), PARAMETER  :: arj01(0:23) = [ &
        0.38179279321690173518_wp, -0.21275636350505321870_wp,  &
        0.16754213407215794187_wp, -0.12853209772196398954_wp,  &
        0.10114405455778847013_wp, -0.9100795343201568859e-1_wp,  &
        0.6401345264656873103e-1_wp, -0.3066963029926754312e-1_wp,  &
        0.1030836525325064201e-1_wp, -0.255670650399956918e-2_wp,  &
        0.48832755805798304e-3_wp, -0.7424935126036077e-4_wp, 0.922260563730861e-5_wp,  &
       -0.95522828307083e-6_wp, 0.8388355845986e-7_wp, -0.633184488858e-8_wp,  &
        0.41560504221e-9_wp, -0.2395529307e-10_wp, 0.122286885e-11_wp, -0.5569711e-13_wp,  &
        0.227820e-14_wp, -0.8417e-16_wp, 0.282e-17_wp, -0.9e-19_wp]

    REAL(wp), PARAMETER  :: arj0a1(0:21) = [ &
        1.24030133037518970827_wp, -0.478125353632280693e-2_wp,  &
        0.6613148891706678e-4_wp, -0.186042740486349e-5_wp, 0.8362735565080e-7_wp,  &
       -0.525857036731e-8_wp, 0.42606363251e-9_wp, -0.4211761024e-10_wp, 0.488946426e-11_wp,  &
       -0.64834929e-12_wp, 0.9617234e-13_wp, -0.1570367e-13_wp, 0.278712e-14_wp, -0.53222e-15_wp, &
        0.10844e-15_wp, -0.2342e-16_wp, 0.533e-17_wp, -0.127e-17_wp, 0.32e-18_wp, -0.8e-19_wp,  &
        0.2e-19_wp, -0.1e-19_wp]

    REAL(wp), PARAMETER  :: arj0a2(0:18) = [ &
        1.99616096301341675339_wp, -0.190379819246668161e-2_wp,  &
        0.1539710927044226e-4_wp, -0.31145088328103e-6_wp, 0.1110850971321e-7_wp,  &
       -0.58666787123e-9_wp, 0.4139926949e-10_wp, -0.365398763e-11_wp, 0.38557568e-12_wp,  &
       -0.4709800e-13_wp, 0.650220e-14_wp, -0.99624e-15_wp, 0.16700e-15_wp, -0.3028e-16_wp,  &
        0.589e-17_wp, -0.122e-17_wp, 0.27e-18_wp, -0.6e-19_wp, 0.1e-19_wp]

    ! Start computation
    x = xvalue
    ind = 1
    IF (x < zero) THEN
      x = -x
      ind = -1
    END IF

    ! Compute the machine-dependent constants.
    temp = EPSILON(zero)
    xhigh = one / temp

    ! Error test
    IF (x > xhigh) THEN
      CALL errprn(fnname,errmsg)
      fn_val = one
      IF (ind == -1) fn_val = -fn_val
      RETURN
    END IF

    ! continue with constants
    t = temp / onehun
    IF (x <= sixten) THEN
      DO  nterm1 = 23, 0, -1
        IF (ABS(arj01(nterm1)) > t) EXIT
      END DO
      xlow = SQRT(twelve*temp)
    ELSE
      DO  nterm2 = 21, 0, -1
        IF (ABS(arj0a1(nterm2)) > t) EXIT
      END DO
      DO  nterm3 = 18, 0, -1
        IF (ABS(arj0a2(nterm3)) > t) EXIT
      END DO
    END IF

    ! Code for 0 <= |x| <= 16
    IF (x <= sixten) THEN
      IF (x < xlow) THEN
        fn_val = x
      ELSE
        t = x * x / one28 - one
        fn_val = x * cheval(nterm1,arj01,t)
      END IF
    ELSE
      ! Code for |x| > 16
      t = five12 / (x*x) - one
      pib41 = pib411 / pib412
      xmpi4 = (x-pib41) - pib42
      temp = COS(xmpi4) * cheval(nterm2,arj0a1,t) / x
      temp = temp - SIN(xmpi4) * cheval(nterm3,arj0a2,t)
      fn_val = one - rt2bpi * temp / SQRT(x)
    END IF
    IF (ind == -1) fn_val = -fn_val
    RETURN
  END FUNCTION j0int


  FUNCTION y0int(xvalue) RESULT(fn_val)
    !! CALGO 757 Integral of the Bessel function \(Y_0(x)\) from \(0\) to \(x\) 
    !
    ! DESCRIPTION:
    !    This function calculates the integral of the Bessel
    !    function Y0, defined as
    !
    !      Y0INT(x) = {integral 0 to x} Y0(t) dt
    !
    !    The code uses Chebyshev expansions whose coefficients are
    !    given to 20 decimal places.
    !
    ! ERROR RETURNS:
    !    If x < 0.0, the function is undefined. An error message
    !    is printed and the function returns the value 0.0.
    !    If the value of x is too large, it is impossible to
    !    accurately compute the trigonometric functions used. An
    !    error message is printed, and the function returns the
    !    value 1.0.
    !
    ! MACHINE-DEPENDENT CONSTANTS:
    !    NTERM1 - The no. of terms to be used from the array
    !             ARJ01. The recommended value is such that
    !                   ABS(ARJ01(NTERM1)) < EPS/100
    !    NTERM2 - The no. of terms to be used from the array
    !             ARY01. The recommended value is such that
    !                   ABS(ARY01(NTERM2)) < EPS/100
    !    NTERM3 - The no. of terms to be used from the array
    !             ARY0A1. The recommended value is such that
    !                   ABS(ARY0A1(NTERM3)) < EPS/100
    !    NTERM4 - The no. of terms to be used from the array
    !             ARY0A2. The recommended value is such that
    !                   ABS(ARY0A2(NTERM4)) < EPS/100
    !    XLOW - The value of x below which
    !                   Y0INT(x) = x*(ln(x) - 0.11593)*2/pi
    !           to machine-precision. The recommended value is
    !                   sqrt(9*EPSNEG)
    !    XHIGH - The value of x above which it is impossible
    !            to calculate (x-pi/4) accurately. The recommended
    !            value is 1/EPSNEG
    !    For values of EPS and EPSNEG, refer to the file MACHCON.TXT
    !    The machine-dependent constants are computed internally by
    !    using the D1MACH subroutine.
    !
    ! INTRINSIC FUNCTIONS USED:
    !   COS, LOG, SIN, SQRT
    !
    ! OTHER MISCFUN SUBROUTINES USED:
    !   CHEVAL, ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: nterm1, nterm2, nterm3, nterm4
    REAL(wp)  :: pib41, t, temp, x, xhigh, xlow, xmpi4
    CHARACTER (LEN=18)  :: ermsg2 = 'ARGUMENT TOO LARGE'
    CHARACTER (LEN=14)  :: ermsg1 = 'ARGUMENT < 0.0'
    CHARACTER (LEN= 6)  :: fnname = 'Y0INT '

    REAL(wp), PARAMETER  :: zero = 0.0_wp, one = 1.0_wp, nine = 9.0_wp,  &
                             sixten = 16.0_wp, onehun = 100.0_wp,   &
                             one28 = 128.0_wp, five12 = 512.0_wp,   &
                             rt2bpi = 0.79788456080286535588_wp,    &
                             pib411 = 201.0_wp, pib412 = 256.0_wp,  &
                             pib42 = 0.24191339744830961566e-3_wp,     &
                             twobpi = 0.63661977236758134308_wp,    &
                             gal2m1 =  -1.11593151565841244881_wp,  &
                             gamln2 =  -0.11593151565841244881_wp

    REAL(wp), PARAMETER  :: arj01(0:23) = [ &
        0.38179279321690173518_wp, -0.21275636350505321870_wp,  &
        0.16754213407215794187_wp, -0.12853209772196398954_wp,  &
        0.10114405455778847013_wp, -0.9100795343201568859e-1_wp,   &
        0.6401345264656873103e-1_wp, -0.3066963029926754312e-1_wp,    &
        0.1030836525325064201e-1_wp, -0.255670650399956918e-2_wp,     &
        0.48832755805798304e-3_wp, -0.7424935126036077e-4_wp,         &
        0.922260563730861e-5_wp, -0.95522828307083e-6_wp, 0.8388355845986e-7_wp, &
       -0.633184488858e-8_wp, 0.41560504221e-9_wp, -0.2395529307e-10_wp, 0.122286885e-11_wp,  &
       -0.5569711e-13_wp, 0.227820e-14_wp, -0.8417e-16_wp, 0.282e-17_wp, -0.9e-19_wp]

    REAL(wp), PARAMETER  :: ary01(0:24) = [ &
        0.54492696302724365490_wp, -0.14957323588684782157_wp,  &
        0.11085634486254842337_wp, -0.9495330018683777109e-1_wp,   &
        0.6820817786991456963e-1_wp, -0.10324653383368200408_wp,   &
        0.10625703287534425491_wp, -0.6258367679961681990e-1_wp,   &
        0.2385645760338293285e-1_wp, -0.644864913015404481e-2_wp,     &
        0.131287082891002331e-2_wp, -0.20988088174989640e-3_wp,       &
        0.2716042484138347e-4_wp, -0.291199114014694e-5_wp, 0.26344333093795e-6_wp, &
       -0.2041172069780e-7_wp, 0.137124781317e-8_wp, -0.8070680792e-10_wp, 0.419883057e-11_wp, &
       -0.19459104e-12_wp, 0.808782e-14_wp, -0.30329e-15_wp, 0.1032e-16_wp, -0.32e-18_wp, &
        0.1e-19_wp]

    REAL(wp), PARAMETER  :: ary0a1(0:21) = [ &
        1.24030133037518970827_wp, -0.478125353632280693e-2_wp,  &
        0.6613148891706678e-4_wp, -0.186042740486349e-5_wp, &
        0.8362735565080e-7_wp, -0.525857036731e-8_wp, 0.42606363251e-9_wp, &
       -0.4211761024e-10_wp, 0.488946426e-11_wp, -0.64834929e-12_wp, 0.9617234e-13_wp, &
       -0.1570367e-13_wp, 0.278712e-14_wp, -0.53222e-15_wp, 0.10844e-15_wp, -0.2342e-16_wp,  &
        0.533e-17_wp, -0.127e-17_wp, 0.32e-18_wp, -0.8e-19_wp, 0.2e-19_wp, -0.1e-19_wp]

    REAL(wp), PARAMETER  :: ary0a2(0:18) = [ &
        1.99616096301341675339_wp, -0.190379819246668161e-2_wp, &
        0.1539710927044226e-4_wp, -0.31145088328103e-6_wp, 0.1110850971321e-7_wp,  &
       -0.58666787123e-9_wp, 0.4139926949e-10_wp, -0.365398763e-11_wp, 0.38557568e-12_wp,  &
       -0.4709800e-13_wp, 0.650220e-14_wp, -0.99624e-15_wp, 0.16700e-15_wp, -0.3028e-16_wp,  &
        0.589e-17_wp, -0.122e-17_wp, 0.27e-18_wp, -0.6e-19_wp, 0.1e-19_wp]

    ! Start computation
    x = xvalue

    ! First error test
    IF (x < zero) THEN
      CALL errprn(fnname,ermsg1)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    temp = EPSILON(zero)
    xhigh = one / temp

    ! Second error test
    IF (x > xhigh) THEN
      CALL errprn(fnname,ermsg2)
      fn_val = zero
      RETURN
    END IF

    ! continue with machine constants
    t = temp / onehun
    IF (x <= sixten) THEN
      DO  nterm1 = 23, 0, -1
        IF (ABS(arj01(nterm1)) > t) EXIT
      END DO
      DO  nterm2 = 24, 0, -1
        IF (ABS(ary01(nterm2)) > t) EXIT
      END DO
      xlow = SQRT(nine*temp)
    ELSE
      DO  nterm3 = 21, 0, -1
        IF (ABS(ary0a1(nterm3)) > t) EXIT
      END DO
      DO  nterm4 = 18, 0, -1
        IF (ABS(ary0a2(nterm4)) > t) EXIT
      END DO
    END IF

    ! Code for 0 <= x <= 16
    IF (x <= sixten) THEN
      IF (x < xlow) THEN
        IF (x == zero) THEN
          fn_val = zero
        ELSE
          fn_val = (LOG(x)+gal2m1) * twobpi * x
        END IF
      ELSE
        t = x * x / one28 - one
        temp = (LOG(x)+gamln2) * cheval(nterm1,arj01,t)
        temp = temp - cheval(nterm2,ary01,t)
        fn_val = twobpi * x * temp
      END IF
    ELSE
      ! Code for x > 16
      t = five12 / (x*x) - one
      pib41 = pib411 / pib412
      xmpi4 = (x-pib41) - pib42
      temp = SIN(xmpi4) * cheval(nterm3,ary0a1,t) / x
      temp = temp + COS(xmpi4) * cheval(nterm4,ary0a2,t)
      fn_val = -rt2bpi * temp / SQRT(x)
    END IF
    RETURN
  END FUNCTION y0int


  FUNCTION i0int(xvalue) RESULT(fn_val)
    !! CALGO 757 Integral of the modified Bessel function
    !! \(I_0(x)\) from \(0\) to \(x\) 
    !
    !   DESCRIPTION:
    !      This program computes the integral of the modified Bessel
    !      function I0(x) using the definition
    !
    !         I0INT(x) = {integral 0 to x} I0(t) dt
    !
    !      The program uses Chebyshev expansions, the coefficients of
    !      which are given to 20 decimal places.
    !
    !   ERROR RETURNS:
    !      If |XVALUE| larger than a certain limit, the value of
    !      I0INT would cause an overflow. If such a situation occurs
    !      the programs prints an error message, and returns the
    !      value sign(XVALUE)*XMAX, where XMAX is the largest
    !      acceptable floating-pt. value.
    !
    !   MACHINE-DEPENDENT CONSTANTS:
    !      NTERM1 - The no. of terms to be used from the array ARI01.
    !                The recommended value is such that
    !                    ABS(ARI01(NTERM1)) < EPS/100
    !      NTERM2 - The no. of terms to be used from the array ARI0A.
    !                The recommended value is such that
    !                    ABS(ARI0A(NTERM2)) < EPS/100
    !      XLOW - The value below which I0INT(x) = x, to machine precision.
    !             The recommended value is
    !                  sqrt(12*EPS).
    !      XHIGH - The value above which overflow will occur. The
    !              recommended value is
    !                  ln(XMAX) + 0.5*ln(ln(XMAX)) + ln(2).
    !
    !      For values of EPS and XMAX refer to the file MACHCON.TXT.
    !
    !      The machine-dependent constants are computed internally by
    !      using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !      EXP , LOG , SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: ind, nterm1, nterm2
    REAL(wp)  :: t, temp, x, xhigh, xlow
    CHARACTER (LEN=26) :: errmsg = 'SIZE OF ARGUMENT TOO LARGE'
    CHARACTER (LEN= 6) :: fnname = 'I0INT '

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, three = 3.0_wp,  &
                             ateen = 18.0_wp, thirt6 = 36.0_wp, onehun = 100.0_wp, &
                             lnr2pi = 0.91893853320467274178_wp

    REAL(wp), PARAMETER  :: ari01(0:28) = [ &
        0.41227906926781516801_wp, -0.34336345150081519562_wp,  &
        0.22667588715751242585_wp, -0.12608164718742260032_wp,  &
        0.6012484628777990271e-1_wp, -0.2480120462913358248e-1_wp,  &
        0.892773389565563897e-2_wp, -0.283253729936696605e-2_wp, 0.79891339041712994e-3_wp, &
       -0.20053933660964890e-3_wp, 0.4416816783014313e-4_wp, -0.822377042246068e-5_wp,  &
        0.120059794219015e-5_wp, -0.11350865004889e-6_wp, 0.69606014466e-9_wp,  &
        0.180622772836e-8_wp, -0.26039481370e-9_wp, -0.166188103e-11_wp, 0.510500232e-11_wp,  &
       -0.41515879e-12_wp, -0.7368138e-13_wp, 0.1279323e-13_wp, 0.103247e-14_wp, -0.30379e-15_wp, &
       -0.1789e-16_wp, 0.673e-17_wp, 0.44e-18_wp, -0.14e-18_wp, -0.1e-19_wp]

    REAL(wp), PARAMETER  :: ari0a(0:33) = [ &
        2.03739654571143287070_wp, 0.1917631647503310248e-1_wp,  &
        0.49923334519288147e-3_wp, 0.2263187103659815e-4_wp, 0.158682108285561e-5_wp,  &
        0.16507855636318e-6_wp, 0.2385058373640e-7_wp, 0.392985182304e-8_wp,  &
        0.46042714199e-9_wp, -0.7072558172e-10_wp, -0.6747183961e-10_wp, -0.2026962001e-10_wp, &
       -0.87320338e-12_wp, 0.175520014e-11_wp, 0.60383944e-12_wp, -0.3977983e-13_wp,  &
       -0.8049048e-13_wp, -0.1158955e-13_wp, 0.827318e-14_wp, 0.282290e-14_wp, -0.77667e-15_wp,  &
       -0.48731e-15_wp, 0.7279e-16_wp, 0.7873e-16_wp, -0.785e-17_wp, -0.1281e-16_wp, 0.121e-17_wp,  &
        0.214e-17_wp, -0.27e-18_wp, -0.36e-18_wp, 0.7e-19_wp, 0.6e-19_wp, -0.2e-19_wp, -0.1e-19_wp]

    ! Start computation
    ind = 1
    x = xvalue
    IF (xvalue < zero) THEN
      ind = -1
      x = -x
    END IF

    ! Compute the machine-dependent constants.
    t = LOG(HUGE(zero))
    xhigh = t + LOG(t) * half - LOG(half)

    ! Error test
    IF (x > xhigh) THEN
      CALL errprn(fnname,errmsg)
      fn_val = EXP(xhigh-lnr2pi-half*LOG(xhigh))
      IF (ind == -1) fn_val = -fn_val
      RETURN
    END IF

    ! Continue with machine-constants
    temp = EPSILON(zero)
    t = temp / onehun
    IF (x <= ateen) THEN
      DO  nterm1 = 28, 0, -1
        IF (ABS(ari01(nterm1)) > t) EXIT
      END DO
      xlow = SQRT(thirt6*temp/three)
    ELSE
      DO  nterm2 = 33, 0, -1
        IF (ABS(ari0a(nterm2)) > t) EXIT
      END DO
    END IF

    ! Code for 0 <= |x| <= 18
    IF (x <= ateen) THEN
      IF (x < xlow) THEN
        fn_val = x
      ELSE
        t = (three*x-ateen) / (x+ateen)
        fn_val = x * EXP(x) * cheval(nterm1,ari01,t)
      END IF
    ELSE
      ! Code for |x| > 18
      t = (thirt6/x-half) - half
      temp = x - half * LOG(x) - lnr2pi + LOG(cheval(nterm2,ari0a,t))
      fn_val = EXP(temp)
    END IF
    IF (ind == -1) fn_val = -fn_val
    RETURN
  END FUNCTION i0int


  FUNCTION k0int(xvalue) RESULT(fn_val)
    !! CALGO 757 Integral of the modified Bessel function
    !! \(K_0(x)\) from \(0\) to \(x\) 
    !
    !   DESCRIPTION:
    !      This function calculates the integral of the modified Bessel function
    !      defined by
    !
    !         K0INT(x) = {integral 0 to x} K0(t) dt
    !
    !      The code uses Chebyshev expansions, whose coefficients are
    !      given to 20 decimal places.
    !
    !   ERROR RETURNS:
    !      If XVALUE < 0.0, the function is undefined. An error message is
    !      printed and the function returns the value 0.0.
    !
    !   MACHINE-DEPENDENT CONSTANTS:
    !      NTERM1 - The no. of terms to be used in the array AK0IN1. The
    !                recommended value is such that
    !                   ABS(AK0IN1(NTERM1)) < EPS/100,
    !      NTERM2 - The no. of terms to be used in the array AK0IN2. The
    !                recommended value is such that
    !                   ABS(AK0IN2(NTERM2)) < EPS/100,
    !      NTERM3 - The no. of terms to be used in the array AK0INA. The
    !                recommended value is such that
    !                   ABS(AK0INA(NTERM3)) < EPS/100,
    !      XLOW - The value below which K0INT = x * ( const - ln(x) ) to
    !             machine precision. The recommended value is
    !                   sqrt (18*EPSNEG).
    !      XHIGH - The value above which K0INT = pi/2 to machine precision.
    !              The recommended value is
    !                   - log (2*EPSNEG)
    !
    !      For values of EPS and EPSNEG refer to the file MACHCON.TXT.
    !
    !      The machine-dependent constants are computed internally by
    !      using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !      EXP , LOG , SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: nterm1, nterm2, nterm3
    REAL(wp)  :: fval, t, temp, x, xhigh, xlow
    CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
    CHARACTER (LEN= 8)  :: fnname = 'K0INT '

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, six = 6.0_wp,  &
                             twelve = 12.0_wp, eightn = 18.0_wp, onehun = 100.0_wp, &
                             const1 = 1.11593151565841244881_wp,   &
                             const2 =  -0.11593151565841244881_wp, &
                             piby2  = 1.57079632679489661923_wp,   &
                             rt2bpi = 0.79788456080286535588_wp

    REAL(wp), PARAMETER  :: ak0in1(0:15) = [ &
        16.79702714464710959477_wp, 9.79134687676889407070_wp,  &
         2.80501316044337939300_wp, 0.45615620531888502068_wp,  &
         0.4716224457074760784e-1_wp, 0.335265148269698289e-2_wp,  &
         0.17335181193874727e-3_wp, 0.679951889364702e-5_wp, 0.20900268359924e-6_wp,  &
         0.516603846976e-8_wp, 0.10485708331e-9_wp, 0.177829320e-11_wp, 0.2556844e-13_wp,  &
         0.31557e-15_wp, 0.338e-17_wp, 0.3e-19_wp]
    REAL(wp), PARAMETER  :: ak0in2(0:15) = [ &
        10.76266558227809174077_wp, 5.62333479849997511550_wp,  &
         1.43543664879290867158_wp, 0.21250410143743896043_wp,  &
         0.2036537393100009554e-1_wp, 0.136023584095623632e-2_wp,  &
         0.6675388699209093e-4_wp, 0.250430035707337e-5_wp, 0.7406423741728e-7_wp,  &
         0.176974704314e-8_wp, 0.3485775254e-10_wp, 0.57544785e-12_wp, 0.807481e-14_wp,  &
         0.9747e-16_wp, 0.102e-17_wp, 0.1e-19_wp]
    REAL(wp), PARAMETER  :: ak0ina(0:27) = [ &
         1.91172065445060453895_wp, -0.4183064565769581085e-1_wp,  &
         0.213352508068147486e-2_wp, -0.15859497284504181e-3_wp, 0.1497624699858351e-4_wp, &
        -0.167955955322241e-5_wp, 0.21495472478804e-6_wp, -0.3058356654790e-7_wp, &
         0.474946413343e-8_wp, -0.79424660432e-9_wp, 0.14156555325e-9_wp,  &
        -0.2667825359e-10_wp, 0.528149717e-11_wp, -0.109263199e-11_wp, 0.23518838e-12_wp,  &
        -0.5247991e-13_wp, 0.1210191e-13_wp, -0.287632e-14_wp, 0.70297e-15_wp, -0.17631e-15_wp, &
         0.4530e-16_wp, -0.1190e-16_wp, 0.319e-17_wp, -0.87e-18_wp, 0.24e-18_wp, -0.7e-19_wp,  &
         0.2e-19_wp, -0.1e-19_wp]

    ! Start computation
    x = xvalue

    ! Error test
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    temp = EPSILON(zero)
    t = temp / onehun
    IF (x <= six) THEN
      DO  nterm1 = 15, 0, -1
        IF (ABS(ak0in1(nterm1)) > t) EXIT
      END DO
      DO  nterm2 = 15, 0, -1
        IF (ABS(ak0in2(nterm2)) > t) EXIT
      END DO
      xlow = SQRT(eightn*temp)
    ELSE
      DO  nterm3 = 27, 0, -1
        IF (ABS(ak0ina(nterm3)) > t) EXIT
      END DO
      xhigh = -LOG(temp+temp)
    END IF

    ! Code for 0 <= XVALUE <= 6
    IF (x <= six) THEN
      IF (x < xlow) THEN
        fval = x
        IF (x > zero) THEN
          fval = fval * (const1-LOG(x))
        END IF
        fn_val = fval
      ELSE
        t = ((x*x)/eightn-half) - half
        fval = (const2+LOG(x)) * cheval(nterm2,ak0in2,t)
        fn_val = x * (cheval(nterm1,ak0in1,t)-fval)
      END IF
      
    ! Code for x > 6
    ELSE
      fval = piby2
      IF (x < xhigh) THEN
        t = (twelve/x-half) - half
        temp = EXP(-x) * cheval(nterm3,ak0ina,t)
        fval = fval - temp / (SQRT(x)*rt2bpi)
      END IF
      fn_val = fval
    END IF
    RETURN
  END FUNCTION k0int


  FUNCTION debye1(xvalue) RESULT(fn_val)
    !! CALGO 757 Debye function \(D_1(x) = 
    !! \frac{1}{x} \int_{0}^{x} \frac{t}{e^t-1} \, dt\)
    !
    !   DEFINITION:
    !      This program calculates the Debye function of order 1, defined as
    !
    !            DEBYE1(x) = [Integral {0 to x} t/(exp(t)-1) dt] / x
    !
    !      The code uses Chebyshev series whose coefficients
    !      are given to 20 decimal places.
    !
    !   ERROR RETURNS:
    !      If XVALUE < 0.0 an error message is printed and the
    !      function returns the value 0.0
    !
    !   MACHINE-DEPENDENT PARAMETERS:
    !      NTERMS - INTEGER - The no. of elements of the array ADEB1.
    !                         The recommended value is such that
    !                             ABS(ADEB1(NTERMS)) < EPS/100 , with
    !                                   1 <= NTERMS <= 18
    !      XLOW - REAL(wp) - The value below which
    !                    DEBYE1 = 1 - x/4 + x*x/36 to machine precision.
    !                    The recommended value is
    !                        SQRT(8*EPSNEG)
    !      XUPPER - REAL(wp) - The value above which
    !                      DEBYE1 = (pi*pi/(6*x)) - exp(-x)(x+1)/x.
    !                      The recommended value is
    !                          -LOG(2*EPS)
    !      XLIM - REAL(wp) - The value above which DEBYE1 = pi*pi/(6*x)
    !                    The recommended value is
    !                          -LOG(XMIN)
    !
    !      For values of EPS, EPSNEG, and XMIN see the file MACHCON.TXT
    !
    !      The machine-dependent constants are computed internally by
    !      using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !      AINT , EXP , INT , LOG , SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: i, nexp, nterms
    REAL(wp)  :: expmx, rk, sum, t, x, xk, xlim, xlow, xupper
    CHARACTER (LEN=17)  :: errmsg = 'ARGUMENT NEGATIVE'
    CHARACTER (LEN= 6)  :: fnname = 'DEBYE1'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, quart = 0.25_wp, half = 0.5_wp,  &
                             one = 1.0_wp, four = 4.0_wp, eight = 8.0_wp,  &
                             nine = 9.0_wp, thirt6 = 36.0_wp, onehun = 100.0_wp, &
                             debinf = 0.60792710185402662866_wp

    REAL(wp), PARAMETER  :: adeb1(0:18) = [ &
        2.40065971903814101941_wp, 0.19372130421893600885_wp,  &
       -0.623291245548957703e-2_wp, 0.35111747702064800e-3_wp, -0.2282224667012310e-4_wp, &
        0.158054678750300e-5_wp, -0.11353781970719e-6_wp, 0.835833611875e-8_wp,  &
       -0.62644247872e-9_wp, 0.4760334890e-10_wp, -0.365741540e-11_wp, 0.28354310e-12_wp,  &
       -0.2214729e-13_wp, 0.174092e-14_wp, -0.13759e-15_wp, 0.1093e-16_wp, -0.87e-18_wp,  &
        0.7e-19_wp, -0.1e-19_wp]

    ! Start computation
    x = xvalue

    ! Check XVALUE >= 0.0
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    t = EPSILON(zero)
    xlow = SQRT(t*eight)
    xupper = -LOG(t+t)
    xlim = -LOG(TINY(zero))
    t = t / onehun
    DO  nterms = 18, 0, -1
      IF (ABS(adeb1(nterms)) > t) GO TO 20
    END DO

    ! Code for x <= 4.0
  20 IF (x <= four) THEN
      IF (x < xlow) THEN
        fn_val = ((x-nine)*x+thirt6) / thirt6
      ELSE
        t = ((x*x/eight)-half) - half
        fn_val = cheval(nterms,adeb1,t) - quart * x
      END IF
    ELSE
      ! Code for x > 4.0
      fn_val = one / (x*debinf)
      IF (x < xlim) THEN
        expmx = EXP(-x)
        IF (x > xupper) THEN
          fn_val = fn_val - expmx * (one+one/x)
        ELSE
          sum = zero
          rk = AINT(xlim/x)
          nexp = INT(rk)
          xk = rk * x
          DO  i = nexp, 1, -1
            t = (one+one/xk) / rk
            sum = sum * expmx + t
            rk = rk - one
            xk = xk - x
          END DO
          fn_val = fn_val - sum * expmx
        END IF
      END IF
    END IF
    RETURN
  END FUNCTION debye1


  FUNCTION debye2(xvalue) RESULT(fn_val)
    !! CALGO 757 Debye function \(D_2(x) = 
    !! \frac{2}{x^2} \int_{0}^{x} \frac{t^2}{e^t-1} \, dt\)
    !
    !   DEFINITION:
    !      This program calculates the Debye function of order 1, defined as
    !
    !            DEBYE2(x) = 2*[Integral {0 to x} t*t/(exp(t)-1) dt] / (x*x)
    !
    !      The code uses Chebyshev series whose coefficients
    !      are given to 20 decimal places.
    !
    !   ERROR RETURNS:
    !      If XVALUE < 0.0 an error message is printed and the
    !      function returns the value 0.0
    !
    !   MACHINE-DEPENDENT PARAMETERS:
    !      NTERMS - INTEGER - The no. of elements of the array ADEB2.
    !                         The recommended value is such that
    !                             ABS(ADEB2(NTERMS)) < EPS/100,
    !                         subject to 1 <= NTERMS <= 18.
    !      XLOW - REAL(wp) - The value below which
    !                    DEBYE2 = 1 - x/3 + x*x/24 to machine precision.
    !                    The recommended value is
    !                        SQRT(8*EPSNEG)
    !      XUPPER - REAL(wp) - The value above which
    !                      DEBYE2 = (4*zeta(3)/x^2) - 2*exp(-x)(x^2+2x+1)/x^2.
    !                      The recommended value is
    !                          -LOG(2*EPS)
    !      XLIM1 - REAL(wp) - The value above which DEBYE2 = 4*zeta(3)/x^2
    !                     The recommended value is
    !                          -LOG(XMIN)
    !      XLIM2 - REAL(wp) - The value above which DEBYE2 = 0.0 to machine
    !                     precision. The recommended value is
    !                           SQRT(4.8/XMIN)
    !
    !      For values of EPS, EPSNEG, and XMIN see the file MACHCON.TXT
    !
    !      The machine-dependent constants are computed internally by
    !      using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !      AINT , EXP , INT , LOG , SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: i, nexp, nterms
    REAL(wp)  :: expmx, rk, sum, t, x, xk, xlim1, xlim2, xlow, xupper
    CHARACTER (LEN=17)  :: errmsg = 'ARGUMENT NEGATIVE'
    CHARACTER (LEN= 6)  :: fnname = 'DEBYE2'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                             two = 2.0_wp, three = 3.0_wp, four = 4.0_wp,  &
                             eight = 8.0_wp, twent4 = 24.0_wp, onehun = 100.0_wp, &
                             debinf = 4.80822761263837714160_wp
    REAL(wp), PARAMETER  :: adeb2(0:18) = [ &
        2.59438102325707702826_wp, 0.28633572045307198337_wp,  &
       -0.1020626561580467129e-1_wp, 0.60491097753468435e-3_wp, -0.4052576589502104e-4_wp, &
        0.286338263288107e-5_wp, -0.20863943030651e-6_wp, 0.1552378758264e-7_wp,  &
       -0.117312800866e-8_wp, 0.8973585888e-10_wp, -0.693176137e-11_wp, 0.53980568e-12_wp,  &
       -0.4232405e-13_wp, 0.333778e-14_wp, -0.26455e-15_wp, 0.2106e-16_wp, -0.168e-17_wp,  &
        0.13e-18_wp, -0.1e-19_wp]

    ! Start computation
    x = xvalue

    ! Check XVALUE >= 0.0
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    t = TINY(zero)
    xlim1 = -LOG(t)
    xlim2 = SQRT(debinf) / SQRT(t)
    t = EPSILON(zero)
    xlow = SQRT(t*eight)
    xupper = -LOG(t+t)
    t = t / onehun
    DO  nterms = 18, 0, -1
      IF (ABS(adeb2(nterms)) > t) GO TO 20
    END DO

    ! Code for x <= 4.0
  20 IF (x <= four) THEN
      IF (x < xlow) THEN
        fn_val = ((x-eight)*x+twent4) / twent4
      ELSE
        t = ((x*x/eight)-half) - half
        fn_val = cheval(nterms,adeb2,t) - x / three
      END IF
    ELSE
      ! Code for x > 4.0
      IF (x > xlim2) THEN
        fn_val = zero
      ELSE
        fn_val = debinf / (x*x)
        IF (x < xlim1) THEN
          expmx = EXP(-x)
          IF (x > xupper) THEN
            sum = ((x+two)*x+two) / (x*x)
          ELSE
            sum = zero
            rk = AINT(xlim1/x)
            nexp = INT(rk)
            xk = rk * x
            DO  i = nexp, 1, -1
              t = (one+two/xk+two/(xk*xk)) / rk
              sum = sum * expmx + t
              rk = rk - one
              xk = xk - x
            END DO
          END IF
          fn_val = fn_val - two * sum * expmx
        END IF
      END IF
    END IF
    RETURN
  END FUNCTION debye2


  FUNCTION debye3(xvalue) RESULT(fn_val)
    !! CALGO 757 Debye function \(D_3(x) = 
    !! \frac{3}{x^3} \int_{0}^{x} \frac{t^3}{e^t-1} \, dt\)
    !
    !   DEFINITION:
    !      This program calculates the Debye function of order 3, defined as
    !
    !            DEBYE3(x) = 3*[Integral {0 to x} t^3/(exp(t)-1) dt] / (x^3)
    !
    !      The code uses Chebyshev series whose coefficients
    !      are given to 20 decimal places.
    !
    !   ERROR RETURNS:
    !      If XVALUE < 0.0 an error message is printed and the
    !      function returns the value 0.0
    !
    !   MACHINE-DEPENDENT PARAMETERS:
    !      NTERMS - INTEGER - The no. of elements of the array ADEB3.
    !                         The recommended value is such that
    !                             ABS(ADEB3(NTERMS)) < EPS/100,
    !                         subject to 1 <= NTERMS <= 18
    !      XLOW - REAL(wp) - The value below which
    !                    DEBYE3 = 1 - 3x/8 + x*x/20 to machine precision.
    !                    The recommended value is
    !                        SQRT(8*EPSNEG)
    !      XUPPER - REAL(wp) - The value above which
    !               DEBYE3 = (18*zeta(4)/x^3) - 3*exp(-x)(x^3+3x^2+6x+6)/x^3.
    !                      The recommended value is
    !                          -LOG(2*EPS)
    !      XLIM1 - REAL(wp) - The value above which DEBYE3 = 18*zeta(4)/x^3
    !                     The recommended value is
    !                          -LOG(XMIN)
    !      XLIM2 - REAL(wp) - The value above which DEBYE3 = 0.0 to machine
    !                     precision. The recommended value is
    !                          CUBE ROOT(19/XMIN)
    !
    !      For values of EPS, EPSNEG, and XMIN see the file MACHCON.TXT
    !
    !      The machine-dependent constants are computed internally by
    !      using the D1MACH subroutine.
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH
    !
    !   INTRINSIC FUNCTIONS USED:
    !      AINT , EXP , INT , LOG , SQRT

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: i, nexp, nterms
    REAL(wp)  :: expmx, rk, sum, t, x, xk, xki, xlim1, xlim2, xlow, xupper
    CHARACTER (LEN=17)  :: errmsg = 'ARGUMENT NEGATIVE'
    CHARACTER (LEN= 6)  :: fnname = 'DEBYE3'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, pt375 = 0.375_wp, half = 0.5_wp,  &
                             one = 1.0_wp, three = 3.0_wp, four = 4.0_wp,  &
                             six = 6.0_wp, sevp5 = 7.5_wp, eight = 8.0_wp,  &
                             twenty = 20.0_wp, onehun = 100.0_wp,  &
                             debinf = 0.51329911273421675946e-1_wp
    REAL(wp), PARAMETER  :: adeb3(0:18) = [ &
        2.70773706832744094526_wp, 0.34006813521109175100_wp,  &
       -0.1294515018444086863e-1_wp, 0.79637553801738164e-3_wp, -0.5463600095908238e-4_wp, &
        0.392430195988049e-5_wp, -0.28940328235386e-6_wp, 0.2173176139625e-7_wp,  &
       -0.165420999498e-8_wp, 0.12727961892e-9_wp, -0.987963459e-11_wp, 0.77250740e-12_wp,  &
       -0.6077972e-13_wp, 0.480759e-14_wp, -0.38204e-15_wp, 0.3048e-16_wp, -0.244e-17_wp,  &
        0.20e-18_wp, -0.2e-19_wp]

    ! Start computation
    x = xvalue

    ! Error test
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    t = TINY(zero)
    xlim1 = -LOG(t)
    xk = one / three
    xki = (one/debinf) ** xk
    rk = t ** xk
    xlim2 = xki / rk
    t = EPSILON(zero)
    xlow = SQRT(t*eight)
    xupper = -LOG(t+t)
    t = t / onehun
    DO  nterms = 18, 0, -1
      IF (ABS(adeb3(nterms)) > t) GO TO 20
    END DO

    ! Code for x <= 4.0
  20 IF (x <= four) THEN
      IF (x < xlow) THEN
        fn_val = ((x-sevp5)*x+twenty) / twenty
      ELSE
        t = ((x*x/eight)-half) - half
        fn_val = cheval(nterms,adeb3,t) - pt375 * x
      END IF
    ELSE
      ! Code for x > 4.0
      IF (x > xlim2) THEN
        fn_val = zero
      ELSE
        fn_val = one / (debinf*x*x*x)
        IF (x < xlim1) THEN
          expmx = EXP(-x)
          IF (x > xupper) THEN
            sum = (((x+three)*x+six)*x+six) / (x*x*x)
          ELSE
            sum = zero
            rk = AINT(xlim1/x)
            nexp = INT(rk)
            xk = rk * x
            DO  i = nexp, 1, -1
              xki = one / xk
              t = (((six*xki+six)*xki+three)*xki+one) / rk
              sum = sum * expmx + t
              rk = rk - one
              xk = xk - x
            END DO
          END IF
          fn_val = fn_val - three * sum * expmx
        END IF
      END IF
    END IF
    RETURN
  END FUNCTION debye3


  FUNCTION debye4(xvalue) RESULT(fn_val)
    !! CALGO 757 Debye function \(D_4(x) = 
    !! \frac{4}{x^4} \int_{0}^{x} \frac{t^4}{e^t-1} \, dt\)
    !
    !   DEFINITION:
    !      This program calculates the Debye function of order 4, defined as
    !
    !            DEBYE4(x) = 4*[Integral {0 to x} t^4/(exp(t)-1) dt] / (x^4)
    !
    !      The code uses Chebyshev series whose coefficients
    !      are given to 20 decimal places.
    !
    !   ERROR RETURNS:
    !      If XVALUE < 0.0 an error message is printed and the
    !      function returns the value 0.0
    !
    !   MACHINE-DEPENDENT PARAMETERS:
    !      NTERMS - INTEGER - The no. of elements of the array ADEB4.
    !                         The recommended value is such that
    !                             ABS(ADEB4(NTERMS)) < EPS/100,
    !                         subject to 1 <= NTERMS <= 18
    !      XLOW - REAL(wp) - The value below which
    !                    DEBYE4 = 1 - 4x/10 + x*x/18 to machine precision.
    !                    The recommended value is
    !                        SQRT(8*EPSNEG)
    !      XUPPER - REAL(wp) - The value above which
    !               DEBYE4=(96*zeta(5)/x^4)-4*exp(-x)(x^4+4x^2+12x^2+24x+24)/x^4.
    !                      The recommended value is
    !                          -LOG(2*EPS)
    !      XLIM1 - REAL(wp) - The value above which DEBYE4 = 96*zeta(5)/x^4
    !                     The recommended value is
    !                          -LOG(XMIN)
    !      XLIM2 - REAL(wp) - The value above which DEBYE4 = 0.0 to machine
    !                     precision. The recommended value is
    !                          FOURTH ROOT(99/XMIN)
    !
    !      For values of EPS, EPSNEG, and XMIN see the file MACHCON.TXT
    !
    !      The machine-dependent constants are computed internally by
    !      using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !      AINT , EXP , INT , LOG , SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: i, nexp, nterms
    REAL(wp)  :: expmx, rk, sum, t, x, xk, xki, xlim1, xlim2, xlow, xupper
    CHARACTER (LEN=17)  :: errmsg = 'ARGUMENT NEGATIVE'
    CHARACTER (LEN= 6)  :: fnname = 'DEBYE4'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                             twopt5 = 2.5_wp, four = 4.0_wp, five = 5.0_wp,  &
                             eight = 8.0_wp, twelve = 12.0_wp, eightn = 18.0_wp,  &
                             twent4 = 24.0_wp, forty5 = 45.0_wp, onehun = 100.0_wp, &
                             debinf = 99.54506449376351292781_wp
    REAL(wp), PARAMETER  :: adeb4(0:18) = [ &
        2.78186941502052346008_wp, 0.37497678352689286364_wp,  &
       -0.1494090739903158326e-1_wp, 0.94567981143704274e-3_wp,  &
       -0.6613291613893255e-4_wp, 0.481563298214449e-5_wp, -0.35880839587593e-6_wp,  &
        0.2716011874160e-7_wp, -0.208070991223e-8_wp, 0.16093838692e-9_wp,  &
       -0.1254709791e-10_wp, 0.98472647e-12_wp, -0.7772369e-13_wp, 0.616483e-14_wp,  &
       -0.49107e-15_wp, 0.3927e-16_wp, -0.315e-17_wp, 0.25e-18_wp, -0.2e-19_wp]

    ! Start computation
    x = xvalue

    ! Check XVALUE >= 0.0
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    t = TINY(zero)
    xlim1 = -LOG(t)
    rk = one / four
    xk = debinf ** rk
    xki = t ** rk
    xlim2 = xk / xki
    t = EPSILON(zero)
    xlow = SQRT(t*eight)
    xupper = -LOG(t+t)
    t = t / onehun
    DO  nterms = 18, 0, -1
      IF (ABS(adeb4(nterms)) > t) GO TO 20
    END DO

    ! Code for x <= 4.0
  20 IF (x <= four) THEN
      IF (x < xlow) THEN
        fn_val = ((twopt5*x-eightn)*x+forty5) / forty5
      ELSE
        t = ((x*x/eight)-half) - half
        fn_val = cheval(nterms,adeb4,t) - (x+x) / five
      END IF
    ELSE
      ! Code for x > 4.0
      IF (x > xlim2) THEN
        fn_val = zero
      ELSE
        t = x * x
        fn_val = (debinf/t) / t
        IF (x < xlim1) THEN
          expmx = EXP(-x)
          IF (x > xupper) THEN
            sum = ((((x+four)*x+twelve)*x+twent4)*x+twent4) / (x*x*x*x )
          ELSE
            sum = zero
            rk = INT(xlim1/x)
            nexp = INT(rk)
            xk = rk * x
            DO  i = nexp, 1, -1
              xki = one / xk
              t = ((((twent4*xki+twent4)*xki+twelve)*xki+four)*xki+one ) / rk
              sum = sum * expmx + t
              rk = rk - one
              xk = xk - x
            END DO
          END IF
          fn_val = fn_val - four * sum * expmx
        END IF
      END IF
    END IF
    RETURN
  END FUNCTION debye4


  FUNCTION strvh0(xvalue) RESULT(fn_val)
    !! CALGO 757 Struve function \(\mathbf{H}_0(x)\)
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
    !
    !     The machine-dependent constants are computed internally by
    !     using the D1MACH subroutine.
    !
    ! INTRINSIC FUNCTIONS USED:
    !     ABS, COS, SIN, SQRT.
    !
    ! OTHER MISCFUN SUBROUTINES USED:
    !     CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: indsgn, nterm1, nterm2, nterm3, nterm4
    REAL(wp)  :: h0as, t, x, xhigh, xlow, xmp4, xsq, y0p, y0q, y0val
    CHARACTER (LEN=26)  :: errmsg = 'ARGUMENT TOO LARGE IN SIZE'
    CHARACTER (LEN= 6)  :: fnname = 'STRVH0'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                             eight = 8.0_wp, eleven = 11.0_wp, twenty = 20.0_wp,  &
                             onehun = 100.0_wp, sixtp5 = 60.5_wp,  &
                             two62 = 262.0_wp, thr2p5 = 302.5_wp,  &
                             piby4 = 0.78539816339744830962_wp,   &
                             rt2bpi = 0.79788456080286535588_wp,  &
                             twobpi = 0.63661977236758134308_wp

    REAL(wp), PARAMETER  :: arrh0(0:19) = [ &
        0.28696487399013225740_wp, -0.25405332681618352305_wp,  &
        0.20774026739323894439_wp, -0.20364029560386585140_wp,  &
        0.12888469086866186016_wp, -0.4825632815622261202e-1_wp,  &
        0.1168629347569001242e-1_wp, -0.198118135642418416e-2_wp,  &
        0.24899138512421286e-3_wp, -0.2418827913785950e-4_wp, 0.187437547993431e-5_wp,  &
       -0.11873346074362e-6_wp, 0.626984943346e-8_wp, -0.28045546793e-9_wp,  &
        0.1076941205e-10_wp, -0.35904793e-12_wp, 0.1049447e-13_wp, -0.27119e-15_wp,  &
        0.624e-17_wp, -0.13e-18_wp]

    REAL(wp), PARAMETER  :: arrh0a(0:20) = [ &
        1.99291885751992305515_wp, -0.384232668701456887e-2_wp,  &
       -0.32871993712353050e-3_wp, -0.2941181203703409e-4_wp, -0.267315351987066e-5_wp,  &
       -0.24681031075013e-6_wp, -0.2295014861143e-7_wp, -0.215682231833e-8_wp,  &
       -0.20303506483e-9_wp, -0.1934575509e-10_wp, -0.182773144e-11_wp, -0.17768424e-12_wp,  &
       -0.1643296e-13_wp, -0.171569e-14_wp, -0.13368e-15_wp, -0.2077e-16_wp, 0.2e-19_wp,  &
       -0.55e-18_wp, 0.10e-18_wp, -0.4e-19_wp, 0.1e-19_wp]

    REAL(wp), PARAMETER  :: ay0asp(0:12) = [ &
        1.99944639402398271568_wp, -0.28650778647031958e-3_wp,  &
       -0.1005072797437620e-4_wp, -0.35835941002463e-6_wp, -0.1287965120531e-7_wp,  &
       -0.46609486636e-9_wp, -0.1693769454e-10_wp, -0.61852269e-12_wp, -0.2261841e-13_wp,  &
       -0.83268e-15_wp, -0.3042e-16_wp, -0.115e-17_wp, -0.4e-19_wp]
    REAL(wp), PARAMETER  :: ay0asq(0:13) = [ &
        1.99542681386828604092_wp, -0.236013192867514472e-2_wp,  &
       -0.7601538908502966e-4_wp, -0.256108871456343e-5_wp, -0.8750292185106e-7_wp,  &
       -0.304304212159e-8_wp, -0.10621428314e-9_wp, -0.377371479e-11_wp, -0.13213687e-12_wp,  &
       -0.488621e-14_wp, -0.15809e-15_wp, -0.762e-17_wp, -0.3e-19_wp, -0.3e-19_wp]

    ! Start computation
    x = xvalue
    indsgn = 1
    IF (x < zero) THEN
      x = -x
      indsgn = -1
    END IF

    ! Compute the machine-dependent constants.
    h0as = EPSILON(zero)
    xhigh = one / (2*EPSILON(zero))

    ! Error test
    IF (ABS(xvalue) > xhigh) THEN
      CALL errprn(fnname,errmsg)
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
    !! CALGO 757 Struve function \(\mathbf{H}_1(x)\)
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
    !
    !     The machine-dependent constants are computed internally by
    !     using the D1MACH subroutine.
    !
    ! INTRINSIC FUNCTIONS USED:
    !     ABS, COS, SIN, SQRT.
    !
    ! OTHER MISCFUN SUBROUTINES USED:
    !     CHEVAL, ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: nterm1, nterm2, nterm3, nterm4
    REAL(wp)  :: h1as, t, x, xhigh, xlow1, xlow2, xm3p4, xsq, y1p, y1q, y1val
    CHARACTER (LEN=26)  :: errmsg = 'ARGUMENT TOO LARGE IN SIZE'
    CHARACTER (LEN= 6)  :: fnname = 'STRVH1'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, eight = 8.0_wp,  &
                             nine = 9.0_wp, fiften = 15.0_wp, twenty = 20.0_wp, &
                             onehun = 100.0_wp, fortp5 = 40.5_wp,  &
                             one82 = 182.0_wp, tw02p5 = 202.5_wp,  &
                             rt2bpi = 0.79788456080286535588_wp,  &
                             thpby4 = 2.35619449019234492885_wp,  &
                             twobpi = 0.63661977236758134308_wp

    REAL(wp), PARAMETER  :: arrh1(0:17) = [ &
        0.17319061083675439319_wp, -0.12606917591352672005_wp,  &
        0.7908576160495357500e-1_wp, -0.3196493222321870820e-1_wp,  &
        0.808040581404918834e-2_wp, -0.136000820693074148e-2_wp,  &
        0.16227148619889471e-3_wp, -0.1442352451485929e-4_wp,  &
        0.99219525734072e-6_wp, -0.5441628049180e-7_wp, 0.243631662563e-8_wp,  &
        -0.9077071338e-10_wp, 0.285926585e-11_wp, -0.7716975e-13_wp,  &
        0.180489e-14_wp, -0.3694e-16_wp, 0.67e-18_wp, -0.1e-19_wp]

    REAL(wp), PARAMETER  :: arrh1a(0:21) = [ &
        2.01083504951473379407_wp, 0.592218610036099903e-2_wp,  &
        0.55274322698414130e-3_wp, 0.5269873856311036e-4_wp, 0.506374522140969e-5_wp,  &
        0.49028736420678e-6_wp, 0.4763540023525e-7_wp, 0.465258652283e-8_wp,  &
        0.45465166081e-9_wp, 0.4472462193e-10_wp, 0.437308292e-11_wp, 0.43568368e-12_wp,  &
        0.4182190e-13_wp, 0.441044e-14_wp, 0.36391e-15_wp, 0.5558e-16_wp, -0.4e-19_wp,  &
        0.163e-17_wp, -0.34e-18_wp, 0.13e-18_wp, -0.4e-19_wp, 0.1e-19_wp]

    REAL(wp), PARAMETER  :: ay1asp(0:14) = [ &
        2.00135240045889396402_wp, 0.71104241596461938e-3_wp,  &
        0.3665977028232449e-4_wp, 0.191301568657728e-5_wp,  &
        0.10046911389777e-6_wp, 0.530401742538e-8_wp, 0.28100886176e-9_wp,  &
        0.1493886051e-10_wp, 0.79578420e-12_wp, 0.4252363e-13_wp, 0.227195e-14_wp,  &
        0.12216e-15_wp, 0.650e-17_wp, 0.36e-18_wp, 0.2e-19_wp]
    REAL(wp), PARAMETER  :: ay1asq(0:15) = [ &
        5.99065109477888189116_wp, -0.489593262336579635e-2_wp,  &
       -0.23238321307070626e-3_wp, -0.1144734723857679e-4_wp, -0.57169926189106e-6_wp,  &
       -0.2895516716917e-7_wp, -0.147513345636e-8_wp, -0.7596537378e-10_wp,  &
       -0.390658184e-11_wp, -0.20464654e-12_wp, -0.1042636e-13_wp, -0.57702e-15_wp,  &
       -0.2550e-16_wp, -0.210e-17_wp, 0.2e-19_wp, -0.2e-19_wp]

    ! Start computation
    x = ABS(xvalue)

    ! Compute the machine-dependent constants.
    xhigh = (half+half) / (2*EPSILON(zero))

    ! Error test
    IF (x > xhigh) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! continue with machine constants
    h1as = EPSILON(zero)
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


  FUNCTION strvl0(xvalue) RESULT(fn_val)
    !! CALGO 757 Struve function \(\mathbf{L}_0(x)\)
    !
    !   DESCRIPTION:
    !      This function calculates the modified Struve function of
    !      order 0, denoted L0(x), defined as the solution of the
    !      second-order equation
    !                  x*D(Df) + Df - x*f  =  2x/pi
    !
    !   ERROR RETURNS:
    !      If the value of |XVALUE| is too large, the result
    !      would cause an floating-pt overflow. An error message
    !      is printed and the function returns the value of
    !      sign(XVALUE)*XMAX where XMAX is the largest possible
    !      floating-pt argument.
    !
    !   MACHINE-DEPENDENT PARAMETERS:
    !      NTERM1 - INTEGER - The no. of terms for the array ARL0.
    !                         The recommended value is such that
    !                             ABS(ARL0(NTERM1)) < EPS/100
    !      NTERM2 - INTEGER - The no. of terms for the array ARL0AS.
    !                         The recommended value is such that
    !                             ABS(ARL0AS(NTERM2)) < EPS/100
    !      NTERM3 - INTEGER - The no. of terms for the array AI0ML0.
    !                         The recommended value is such that
    !                             ABS(AI0ML0(NTERM3)) < EPS/100
    !      XLOW - REAL(wp) - The value of x below which L0(x) = 2*x/pi
    !                    to machine precision. The recommended value is
    !                             3*SQRT(EPS)
    !      XHIGH1 - REAL(wp) - The value beyond which the Chebyshev series
    !                      in the asymptotic expansion of I0 - L0 gives
    !                      1.0 to machine precision. The recommended value
    !                      is   SQRT( 30/EPSNEG )
    !      XHIGH2 - REAL(wp) - The value beyond which the Chebyshev series
    !                      in the asymptotic expansion of I0 gives 1.0
    !                      to machine precision. The recommended value
    !                      is   28 / EPSNEG
    !      XMAX - REAL(wp) - The value of XMAX, where XMAX is the
    !                    largest possible floating-pt argument.
    !                    This is used to prevent overflow.
    !      For values of EPS, EPSNEG and XMAX the user should refer
    !      to the file MACHCON.TXT
    !      The machine-dependent constants are computed internally by
    !      using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !      EXP , LOG , SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: indsgn, nterm1, nterm2, nterm3
    REAL(wp)  :: ch1, ch2, t, test, x, xhigh1, xhigh2, xlow, xmax, xsq
    CHARACTER (LEN=24)  :: errmsg = 'ARGUMENT CAUSES OVERFLOW'
    CHARACTER (LEN= 6)  :: fnname = 'STRVL0'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, one = 1.0_wp, two = 2.0_wp,  &
                             four = 4.0_wp, sixten = 16.0_wp, twent4 = 24.0_wp,  &
                             twent8 = 28.0_wp, onehun = 100.0_wp,  &
                             two88 = 288.0_wp, atehun = 800.0_wp,  &
                             lnr2pi = 0.91893853320467274178_wp,  &
                             twobpi = 0.63661977236758134308_wp

    REAL(wp), PARAMETER  :: arl0(0:27) = [ &
        0.42127458349979924863_wp, -0.33859536391220612188_wp,  &
        0.21898994812710716064_wp, -0.12349482820713185712_wp,  &
        0.6214209793866958440e-1_wp, -0.2817806028109547545e-1_wp,  &
        0.1157419676638091209e-1_wp, -0.431658574306921179e-2_wp,  &
        0.146142349907298329e-2_wp, -0.44794211805461478e-3_wp, 0.12364746105943761e-3_wp, &
       -0.3049028334797044e-4_wp, 0.663941401521146e-5_wp, -0.125538357703889e-5_wp,  &
        0.20073446451228e-6_wp, -0.2588260170637e-7_wp, 0.241143742758e-8_wp,  &
       -0.10159674352e-9_wp, -0.1202430736e-10_wp, 0.262906137e-11_wp, -0.15313190e-12_wp,  &
       -0.1574760e-13_wp, 0.315635e-14_wp, -0.4096e-16_wp, -0.3620e-16_wp, 0.239e-17_wp,  &
        0.36e-18_wp, -0.4e-19_wp]

    REAL(wp), PARAMETER  :: arl0as(0:15) = [ &
        2.00861308235605888600_wp, 0.403737966500438470e-2_wp,  &
       -0.25199480286580267e-3_wp, 0.1605736682811176e-4_wp, -0.103692182473444e-5_wp,  &
        0.6765578876305e-7_wp, -0.444999906756e-8_wp, 0.29468889228e-9_wp,  &
       -0.1962180522e-10_wp, 0.131330306e-11_wp, -0.8819190e-13_wp, 0.595376e-14_wp,  &
       -0.40389e-15_wp, 0.2651e-16_wp, -0.208e-17_wp, 0.11e-18_wp]

    REAL(wp), PARAMETER  :: ai0ml0(0:23) = [ &
        2.00326510241160643125_wp, 0.195206851576492081e-2_wp,  &
        0.38239523569908328e-3_wp, 0.7534280817054436e-4_wp, 0.1495957655897078e-4_wp,  &
        0.299940531210557e-5_wp, 0.60769604822459e-6_wp, 0.12399495544506e-6_wp,  &
        0.2523262552649e-7_wp, 0.504634857332e-8_wp, 0.97913236230e-9_wp, 0.18389115241e-9_wp, &
        0.3376309278e-10_wp, 0.611179703e-11_wp, 0.108472972e-11_wp, 0.18861271e-12_wp,  &
        0.3280345e-13_wp, 0.565647e-14_wp, 0.93300e-15_wp, 0.15881e-15_wp, 0.2791e-16_wp,  &
        0.389e-17_wp, 0.70e-18_wp, 0.16e-18_wp]

    ! Start computation
    x = xvalue
    indsgn = 1
    IF (x < zero) THEN
      x = -x
      indsgn = -1
    END IF

    ! Compute the machine-dependent constants.
    test = EPSILON(zero)
    t = test / onehun
    IF (x <= sixten) THEN
      DO  nterm1 = 27, 0, -1
        IF (ABS(arl0(nterm1)) > t) EXIT
      END DO
      xlow = (one+two) * SQRT(test)
    ELSE
      DO  nterm2 = 15, 0, -1
        IF (ABS(arl0as(nterm2)) > t) EXIT
      END DO
      DO  nterm3 = 23, 0, -1
        IF (ABS(ai0ml0(nterm3)) > t) EXIT
      END DO
      xmax = HUGE(zero)
      xhigh1 = SQRT((twent8+two)/test)
      xhigh2 = twent8 / test
    END IF

    ! Code for |xvalue| <= 16
    IF (x <= sixten) THEN
      IF (x < xlow) THEN
        fn_val = twobpi * x
      ELSE
        t = (four*x-twent4) / (x+twent4)
        fn_val = twobpi * x * cheval(nterm1,arl0,t) * EXP(x)
      END IF
    ELSE
      ! Code for |xvalue| > 16
      IF (x > xhigh2) THEN
        ch1 = one
      ELSE
        t = (x-twent8) / (four-x)
        ch1 = cheval(nterm2,arl0as,t)
      END IF
      IF (x > xhigh1) THEN
        ch2 = one
      ELSE
        xsq = x * x
        t = (atehun-xsq) / (two88+xsq)
        ch2 = cheval(nterm3,ai0ml0,t)
      END IF
      test = LOG(ch1) - lnr2pi - LOG(x) / two + x
      IF (test > LOG(xmax)) THEN
        CALL errprn(fnname,errmsg)
        fn_val = xmax
      ELSE
        fn_val = EXP(test) - twobpi * ch2 / x
      END IF
    END IF
    IF (indsgn == -1) fn_val = -fn_val
    RETURN
  END FUNCTION strvl0


  FUNCTION strvl1(xvalue) RESULT(fn_val)
    !! CALGO 757 Struve function \(\mathbf{L}_1(x)\)
    !
    !   DESCRIPTION:
    !      This function calculates the modified Struve function of
    !      order 1, denoted L1(x), defined as the solution of
    !
    !               x*x*D(Df) + x*Df - (x*x+1)f = 2*x*x/pi
    !
    !   ERROR RETURNS:
    !      If the value of |XVALUE| is too large, the result
    !      would cause an floating-pt overflow. An error message
    !      is printed and the function returns the value of
    !      sign(XVALUE)*XMAX where XMAX is the largest possible
    !      floating-pt argument.
    !
    !   MACHINE-DEPENDENT PARAMETERS:
    !      NTERM1 - INTEGER - The no. of terms for the array ARL1.
    !                         The recommended value is such that
    !                             ABS(ARL1(NTERM1)) < EPS/100
    !      NTERM2 - INTEGER - The no. of terms for the array ARL1AS.
    !                         The recommended value is such that
    !                             ABS(ARL1AS(NTERM2)) < EPS/100
    !      NTERM3 - INTEGER - The no. of terms for the array AI1ML1.
    !                         The recommended value is such that
    !                             ABS(AI1ML1(NTERM3)) < EPS/100
    !      XLOW1 - REAL(wp) - The value of x below which L1(x) = 2*x*x/(3*pi)
    !                     to machine precision. The recommended value is
    !                              SQRT(15*EPS)
    !      XLOW2 - REAL(wp) - The value of x below which L1(x) set to 0.0.
    !                     This is used to prevent underflow. The
    !                     recommended value is
    !                              SQRT(5*XMIN)
    !      XHIGH1 - REAL(wp) - The value of |x| above which the Chebyshev
    !                      series in the asymptotic expansion of I1
    !                      equals 1.0 to machine precision. The
    !                      recommended value is  SQRT( 30 / EPSNEG ).
    !      XHIGH2 - REAL(wp) - The value of |x| above which the Chebyshev
    !                      series in the asymptotic expansion of I1 - L1
    !                      equals 1.0 to machine precision. The recommended
    !                      value is   30 / EPSNEG.
    !      XMAX - REAL(wp) - The value of XMAX, where XMAX is the
    !                    largest possible floating-pt argument.
    !                    This is used to prevent overflow.
    !      For values of EPS, EPSNEG, XMIN, and XMAX the user should refer
    !      to the file MACHCON.TXT
    !      The machine-dependent constants are computed internally by
    !      using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !      EXP , LOG , SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: nterm1, nterm2, nterm3
    REAL(wp)  :: ch1, ch2, t, test, x, xhigh1, xhigh2, xlow1, xlow2, xmax, xsq
    CHARACTER (LEN=24)  :: errmsg = 'ARGUMENT CAUSES OVERFLOW'
    CHARACTER (LEN= 6)  :: fnname = 'STRVL1'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, one = 1.0_wp, two = 2.0_wp,  &
                             four = 4.0_wp, sixten = 16.0_wp, twent4 = 24.0_wp,  &
                             thirty = 30.0_wp, onehun = 100.0_wp, two88 = 288.0_wp, &
                             atehun = 800.0_wp,  &
                             lnr2pi = 0.91893853320467274178_wp, &
                             pi3by2 = 4.71238898038468985769_wp, &
                             twobpi = 0.63661977236758134308_wp

    REAL(wp), PARAMETER  :: arl1(0:26) = [ &
        0.38996027351229538208_wp, -0.33658096101975749366_wp,  &
        0.23012467912501645616_wp, -0.13121594007960832327_wp,  &
        0.6425922289912846518e-1_wp, -0.2750032950616635833e-1_wp,  &
        0.1040234148637208871e-1_wp, -0.350532294936388080e-2_wp,  &
        0.105748498421439717e-2_wp, -0.28609426403666558e-3_wp, 0.6925708785942208e-4_wp,  &
       -0.1489693951122717e-4_wp, 0.281035582597128e-5_wp, -0.45503879297776e-6_wp,  &
        0.6090171561770e-7_wp, -0.623543724808e-8_wp, 0.38430012067e-9_wp, 0.790543916e-11_wp, &
       -0.489824083e-11_wp, 0.46356884e-12_wp, 0.684205e-14_wp, -0.569748e-14_wp,  &
        0.35324e-15_wp, 0.4244e-16_wp, -0.644e-17_wp, -0.21e-18_wp, 0.9e-19_wp]

    REAL(wp), PARAMETER  :: arl1as(0:16) = [ &
        1.97540378441652356868_wp, -0.1195130555088294181e-1_wp,  &
        0.33639485269196046e-3_wp, -0.1009115655481549e-4_wp, 0.30638951321998e-6_wp,  &
       -0.953704370396e-8_wp, 0.29524735558e-9_wp, -0.951078318e-11_wp, 0.28203667e-12_wp,  &
       -0.1134175e-13_wp, 0.147e-17_wp, -0.6232e-16_wp, -0.751e-17_wp, -0.17e-18_wp, 0.51e-18_wp,  &
        0.23e-18_wp, 0.5e-19_wp]

    REAL(wp), PARAMETER  :: ai1ml1(0:25) = [ &
        1.99679361896789136501_wp, -0.190663261409686132e-2_wp,  &
       -0.36094622410174481e-3_wp, -0.6841847304599820e-4_wp, -0.1299008228509426e-4_wp,  &
       -0.247152188705765e-5_wp, -0.47147839691972e-6_wp, -0.9020819982592e-7_wp,  &
       -0.1730458637504e-7_wp, -0.332323670159e-8_wp, -0.63736421735e-9_wp,  &
       -0.12180239756e-9_wp, -0.2317346832e-10_wp, -0.439068833e-11_wp, -0.82847110e-12_wp,  &
       -0.15562249e-12_wp, -0.2913112e-13_wp, -0.543965e-14_wp, -0.101177e-14_wp,  &
       -0.18767e-15_wp, -0.3484e-16_wp, -0.643e-17_wp, -0.118e-17_wp, -0.22e-18_wp, -0.4e-19_wp,  &
       -0.1e-19_wp]

    ! START CALCULATION
    x = ABS(xvalue)

    ! Compute the machine-dependent constants.
    test = EPSILON(zero)
    t = test / onehun
    IF (x <= sixten) THEN
      DO  nterm1 = 26, 0, -1
        IF (ABS(arl1(nterm1)) > t) EXIT
      END DO
      xlow1 = SQRT(thirty*test/two)
      xlow2 = SQRT((four+one)*TINY(zero))
    ELSE
      DO  nterm2 = 16, 0, -1
        IF (ABS(arl1as(nterm2)) > t) EXIT
      END DO
      DO  nterm3 = 25, 0, -1
        IF (ABS(ai1ml1(nterm3)) > t) EXIT
      END DO
      xmax = HUGE(zero)
      xhigh2 = thirty / test
      xhigh1 = SQRT(xhigh2)
    END IF

    ! CODE FOR |XVALUE| <= 16
    IF (x <= sixten) THEN
      IF (x <= xlow2) THEN
        fn_val = zero
      ELSE
        xsq = x * x
        IF (x < xlow1) THEN
          fn_val = xsq / pi3by2
        ELSE
          t = (four*x-twent4) / (x+twent4)
          fn_val = xsq * cheval(nterm1,arl1,t) * EXP(x) / pi3by2
        END IF
      END IF
    ELSE
      ! CODE FOR |XVALUE| > 16
      IF (x > xhigh2) THEN
        ch1 = one
      ELSE
        t = (x-thirty) / (two-x)
        ch1 = cheval(nterm2,arl1as,t)
      END IF
      IF (x > xhigh1) THEN
        ch2 = one
      ELSE
        xsq = x * x
        t = (atehun-xsq) / (two88+xsq)
        ch2 = cheval(nterm3,ai1ml1,t)
      END IF
      test = LOG(ch1) - lnr2pi - LOG(x) / two + x
      IF (test > LOG(xmax)) THEN
        CALL errprn(fnname,errmsg)
        fn_val = xmax
      ELSE
        fn_val = EXP(test) - twobpi * ch2
      END IF
    END IF
    RETURN
  END FUNCTION strvl1


  FUNCTION i0ml0(xvalue) RESULT(fn_val)
    !! CALGO 757 Difference between the modified Bessel and Struve functions
    !! \(I_0(x) - \mathbf{L}_0(x)\)
    !
    !   DESCRIPTION:
    !      This program calculates the function I0ML0 defined as
    !
    !                I0ML0(x) = I0(x) - L0(x)
    !
    !      where I0(x) is the modified Bessel function of the first kind of
    !      order 0, and L0(x) is the modified Struve function of order 0.
    !
    !      The code uses Chebyshev expansions with the coefficients
    !      given to an accuracy of 20D.
    !
    !   ERROR RETURNS:
    !      The coefficients are only suitable for XVALUE >= 0.0. If
    !      XVALUE < 0.0, an error message is printed and the function
    !      returns the value 0.0
    !
    !   MACHINE-DEPENDENT PARAMETERS:
    !      NTERM1 - INTEGER - The number of terms required for the array
    !                         AI0L0. The recommended value is such that
    !                              ABS(AI0L0(NTERM1)) < EPS/100
    !      NTERM2 - INTEGER - The number of terms required for the array
    !                         AI0L0A. The recommended value is such that
    !                              ABS(AI0L0A(NTERM2)) < EPS/100
    !      XLOW - REAL(wp) - The value below which I0ML0(x) = 1 to machine
    !                    precision. The recommended value is
    !                               EPSNEG
    !      XHIGH - REAL(wp) - The value above which I0ML0(x) = 2/(pi*x) to
    !                     machine precision. The recommended value is
    !                               SQRT(800/EPS)
    !
    !      For values of EPS, and EPSNEG see the file MACHCON.TXT
    !
    !      The machine-dependent constants are computed internally by
    !      using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !      SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: nterm1, nterm2
    REAL(wp)  :: t, x, xhigh, xlow, xsq
    CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
    CHARACTER (LEN= 6)  :: fnname = 'I0ML0 '

    REAL(wp), PARAMETER  :: zero = 0.0_wp, one = 1.0_wp, six = 6.0_wp,  &
                             sixten = 16.0_wp, forty = 40.0_wp, two88 = 288.0_wp, &
                             onehun = 100.0_wp, atehun = 800.0_wp,  &
                             twobpi = 0.63661977236758134308_wp

    REAL(wp), PARAMETER  :: ai0l0(0:23) = [ &
        0.52468736791485599138_wp, -0.35612460699650586196_wp,  &
        0.20487202864009927687_wp, -0.10418640520402693629_wp,  &
        0.4634211095548429228e-1_wp, -0.1790587192403498630e-1_wp,  &
        0.597968695481143177e-2_wp, -0.171777547693565429e-2_wp,  &
        0.42204654469171422e-3_wp, -0.8796178522094125e-4_wp, 0.1535434234869223e-4_wp,  &
       -0.219780769584743e-5_wp, 0.24820683936666e-6_wp, -0.2032706035607e-7_wp,  &
        0.90984198421e-9_wp, 0.2561793929e-10_wp, -0.710609790e-11_wp, 0.32716960e-12_wp,  &
        0.2300215e-13_wp, -0.292109e-14_wp, -0.3566e-16_wp, 0.1832e-16_wp, -0.10e-18_wp,  &
       -0.11e-18_wp]

    REAL(wp), PARAMETER  :: ai0l0a(0:23) = [ &
        2.00326510241160643125_wp, 0.195206851576492081e-2_wp,  &
        0.38239523569908328e-3_wp, 0.7534280817054436e-4_wp,  &
        0.1495957655897078e-4_wp, 0.299940531210557e-5_wp, 0.60769604822459e-6_wp,  &
        0.12399495544506e-6_wp, 0.2523262552649e-7_wp, 0.504634857332e-8_wp,  &
        0.97913236230e-9_wp, 0.18389115241e-9_wp, 0.3376309278e-10_wp, 0.611179703e-11_wp,  &
        0.108472972e-11_wp, 0.18861271e-12_wp, 0.3280345e-13_wp, 0.565647e-14_wp,  &
        0.93300e-15_wp, 0.15881e-15_wp, 0.2791e-16_wp, 0.389e-17_wp, 0.70e-18_wp, 0.16e-18_wp]

    ! Start computation
    x = xvalue

    ! Error test
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    xsq = EPSILON(zero)
    t = xsq / onehun
    IF (x <= sixten) THEN
      DO  nterm1 = 23, 0, -1
        IF (ABS(ai0l0(nterm1)) > t) EXIT
      END DO
      xlow = xsq
    ELSE
      DO  nterm2 = 23, 0, -1
        IF (ABS(ai0l0a(nterm2)) > t) EXIT
      END DO
      xhigh = SQRT(atehun/xsq)
    END IF

    ! Code for x <= 16
    IF (x <= sixten) THEN
      IF (x < xlow) THEN
        fn_val = one
        RETURN
      ELSE
        t = (six*x-forty) / (x+forty)
        fn_val = cheval(nterm1,ai0l0,t)
        RETURN
      END IF
    ELSE
      ! Code for x > 16
      IF (x > xhigh) THEN
        fn_val = twobpi / x
      ELSE
        xsq = x * x
        t = (atehun-xsq) / (two88+xsq)
        fn_val = cheval(nterm2,ai0l0a,t) * twobpi / x
      END IF
    END IF
    RETURN
  END FUNCTION i0ml0


  FUNCTION i1ml1(xvalue) RESULT(fn_val)
    !! CALGO 757 Difference between the modified Bessel and Struve functions
    !! \(I_1(x) - \mathbf{L}_1(x)\)
    !
    ! DESCRIPTION:
    !
    !    This program calculates the function I1ML1 defined as
    !
    !              I1ML1(x) = I1(x) - L1(x)
    !
    !    where I1(x) is the modified Bessel function of the first kind of
    !    order 1, and L1(x) is the modified Struve function of order 1.
    !
    !    The code uses Chebyshev expansions with the coefficients
    !    given to an accuracy of 20D.
    !
    ! ERROR RETURNS:
    !    The coefficients are only suitable for XVALUE >= 0.0. If
    !    XVALUE < 0.0, an error message is printed and the function
    !    returns the value 0.0
    !
    ! MACHINE-DEPENDENT PARAMETERS:
    !    NTERM1 - INTEGER - The number of terms required for the array
    !                       AI1L1. The recommended value is such that
    !                            ABS(AI1L1(NTERM1)) < EPS/100
    !    NTERM2 - INTEGER - The number of terms required for the array
    !                       AI1L1A. The recommended value is such that
    !                            ABS(AI1L1A(NTERM2)) < EPS/100
    !    XLOW - REAL(wp) - The value below which I1ML1(x) = x/2 to machine
    !                  precision. The recommended value is
    !                             2*EPSNEG
    !    XHIGH - REAL(wp) - The value above which I1ML1(x) = 2/pi to
    !                   machine precision. The recommended value is
    !                             SQRT(800/EPS)
    !
    !    For values of EPS, and EPSNEG see the file MACHCON.TXT
    !
    !    The machine-dependent constants are computed internally by
    !    using the D1MACH subroutine.
    !
    ! INTRINSIC FUNCTIONS USED:
    !    ABS , SQRT
    !
    ! OTHER MISCFUN SUBROUTINES USED:
    !        CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: nterm1, nterm2
    REAL(wp)  :: t, x, xhigh, xlow, xsq
    CHARACTER (LEN=14) :: errmsg = 'ARGUMENT < 0.0'
    CHARACTER (LEN= 6) :: fnname = 'I1ML1 '

    REAL(wp), PARAMETER  :: zero = 0.0_wp, two = 2.0_wp, six = 6.0_wp,  &
                             sixten = 16.0_wp, forty = 40.0_wp, onehun = 100.0_wp, &
                             two88 = 288.0_wp, atehun = 800.0_wp, &
                             twobpi = 0.63661977236758134308_wp

    REAL(wp), PARAMETER  :: ai1l1(0:23) = [ &
        0.67536369062350576137_wp, -0.38134971097266559040_wp,  &
        0.17452170775133943559_wp, -0.7062105887235025061e-1_wp,  &
        0.2517341413558803702e-1_wp, -0.787098561606423321e-2_wp,  &
        0.214814368651922006e-2_wp, -0.50862199717906236e-3_wp, 0.10362608280442330e-3_wp, &
       -0.1795447212057247e-4_wp, 0.259788274515414e-5_wp, -0.30442406324667e-6_wp,  &
        0.2720239894766e-7_wp, -0.158126144190e-8_wp, 0.1816209172e-10_wp,  &
        0.647967659e-11_wp, -0.54113290e-12_wp, -0.308311e-14_wp, 0.305638e-14_wp,  &
       -0.9717e-16_wp, -0.1422e-16_wp, 0.84e-18_wp, 0.7e-19_wp, -0.1e-19_wp]

    REAL(wp), PARAMETER  :: ai1l1a(0:25) = [ &
        1.99679361896789136501_wp, -0.190663261409686132e-2_wp,  &
       -0.36094622410174481e-3_wp, -0.6841847304599820e-4_wp, -0.1299008228509426e-4_wp,  &
       -0.247152188705765e-5_wp, -0.47147839691972e-6_wp, -0.9020819982592e-7_wp,  &
       -0.1730458637504e-7_wp, -0.332323670159e-8_wp, -0.63736421735e-9_wp,  &
       -0.12180239756e-9_wp, -0.2317346832e-10_wp, -0.439068833e-11_wp, -0.82847110e-12_wp,  &
       -0.15562249e-12_wp, -0.2913112e-13_wp, -0.543965e-14_wp, -0.101177e-14_wp,  &
       -0.18767e-15_wp, -0.3484e-16_wp, -0.643e-17_wp, -0.118e-17_wp, -0.22e-18_wp, -0.4e-19_wp,  &
       -0.1e-19_wp]

    ! Start computation
    x = xvalue

    ! Error test
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    xsq = EPSILON(zero)
    t = xsq / onehun
    IF (x <= sixten) THEN
      DO  nterm1 = 23, 0, -1
        IF (ABS(ai1l1(nterm1)) > t) EXIT
      END DO
      xlow = xsq + xsq
    ELSE
      DO  nterm2 = 25, 0, -1
        IF (ABS(ai1l1a(nterm2)) > t) EXIT
      END DO
      xhigh = SQRT(atehun/xsq)
    END IF

    ! Code for x <= 16
    IF (x <= sixten) THEN
      IF (x < xlow) THEN
        fn_val = x / two
        RETURN
      ELSE
        t = (six*x-forty) / (x+forty)
        fn_val = cheval(nterm1,ai1l1,t) * x / two
        RETURN
      END IF
    ELSE
      ! Code for x > 16
      IF (x > xhigh) THEN
        fn_val = twobpi
      ELSE
        xsq = x * x
        t = (atehun-xsq) / (two88+xsq)
        fn_val = cheval(nterm2,ai1l1a,t) * twobpi
      END IF
    END IF
    RETURN
  END FUNCTION i1ml1


  FUNCTION synch1(xvalue) RESULT(fn_val)
    !! CALGO 757 Synchrotron radiation function
    !! \(f_1(x) = x \int_{x}^{\infty} K_{5/3}(t) \, dt\)
    !
    !   DESCRIPTION:
    !      This function calculates the synchrotron radiation function defined as
    !
    !         SYNCH1(x) = x * Integral{x to inf} K(5/3)(t) dt,
    !
    !      where K(5/3) is a modified Bessel function of order 5/3.
    !
    !      The code uses Chebyshev expansions, the coefficients of which
    !      are given to 20 decimal places.
    !
    !   ERROR RETURNS:
    !      The function is undefined if x < 0.0. If XVALUE < 0.0,
    !      an error message is printed and the function returns the value 0.0.
    !
    !   MACHINE-DEPENDENT CONSTANTS:
    !      NTERM1 - INTEGER - The no. of terms needed from the array
    !                         ASYNC1. The recommended value is such that
    !                            ABS(ASYNC1(NTERM1)) < EPS/100.
    !      NTERM2 - INTEGER - The no. of terms needed from the array
    !                         ASYNC2. The recommended value is such that
    !                            ABS(ASYNC2(NTERM2)) < EPS/100.
    !      NTERM3 - INTEGER - The no. of terms needed from the array
    !                         ASYNCA. The recommended value is such that
    !                            ABS(ASYNCA(NTERM3)) < EPS/100.
    !      XLOW - REAL(wp) - The value below which
    !                        SYNCH1(x) = 2.14952.. * (x**(1/3))
    !                    to machine precision. The recommended value
    !                    is     sqrt (8*EPSNEG)
    !      XHIGH1 - REAL(wp) - The value above which
    !                          SYNCH1(x) = 0.0
    !                      to machine precision. The recommended value
    !                      is     -8*LN(XMIN)/7
    !      XHIGH2 - REAL(wp) - The value of LN(XMIN). This is used
    !                      to prevent underflow in calculations
    !                      for large x.
    !     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
    !     The machine-dependent constants are computed internally by
    !     using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !      EXP , LOG , SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: nterm1, nterm2, nterm3
    REAL(wp)  :: cheb1, cheb2, t, x, xhigh1, xhigh2, xlow, xpowth
    CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
    CHARACTER (LEN= 6)  :: fnname = 'SYNCH1'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                             three = 3.0_wp, four = 4.0_wp, eight = 8.0_wp,  &
                             twelve = 12.0_wp, onehun = 100.0_wp,  &
                             conlow = 2.14952824153447863671_wp,  &
                             pibrt3 = 1.81379936423421785059_wp,  &
                             lnrtp2 = 0.22579135264472743236_wp

    REAL(wp), PARAMETER  :: async1(0:13) = [ &
        30.36468298250107627340_wp, 17.07939527740839457449_wp,  &
        4.56013213354507288887_wp, 0.54928124673041997963_wp,  &
        0.3729760750693011724e-1_wp, 0.161362430201041242e-2_wp,  &
        0.4819167721203707e-4_wp, 0.105124252889384e-5_wp,  &
        0.1746385046697e-7_wp, 0.22815486544e-9_wp, 0.240443082e-11_wp,  &
        0.2086588e-13_wp, 0.15167e-15_wp, 0.94e-18_wp]

    REAL(wp), PARAMETER  :: async2(0:11) = [ &
        0.44907216235326608443_wp, 0.8983536779941872179e-1_wp,  &
        0.810445737721512894e-2_wp, 0.42617169910891619e-3_wp,  &
        0.1476096312707460e-4_wp, 0.36286336153998e-6_wp, 0.666348074984e-8_wp,  &
        0.9490771655e-10_wp, 0.107912491e-11_wp, 0.1002201e-13_wp, 0.7745e-16_wp, 0.51e-18_wp]

    REAL(wp), PARAMETER  :: asynca(0:24) = [ &
        2.13293051613550009848_wp, 0.7413528649542002401e-1_wp,  &
        0.869680999099641978e-2_wp, 0.117038262487756921e-2_wp, 0.16451057986191915e-3_wp, &
        0.2402010214206403e-4_wp, 0.358277563893885e-5_wp, 0.54477476269837e-6_wp,  &
        0.8388028561957e-7_wp, 0.1306988268416e-7_wp, 0.205309907144e-8_wp,  &
        0.32518753688e-9_wp, 0.5179140412e-10_wp, 0.830029881e-11_wp, 0.133527277e-11_wp,  &
        0.21591498e-12_wp, 0.3499673e-13_wp, 0.569942e-14_wp, 0.92906e-15_wp, 0.15222e-15_wp,  &
        0.2491e-16_wp, 0.411e-17_wp, 0.67e-18_wp, 0.11e-18_wp, 0.2e-19_wp]

    ! Start calculation
    x = xvalue
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    cheb1 = EPSILON(zero)
    t = cheb1 / onehun
    IF (x <= four) THEN
      DO  nterm1 = 13, 0, -1
        IF (ABS(async1(nterm1)) > t) EXIT
      END DO
      DO  nterm2 = 11, 0, -1
        IF (ABS(async2(nterm2)) > t) EXIT
      END DO
      xlow = SQRT(eight*cheb1)
    ELSE
      DO  nterm3 = 24, 0, -1
        IF (ABS(asynca(nterm3)) > t) EXIT
      END DO
      xhigh2 = LOG(TINY(zero))
      xhigh1 = -eight * xhigh2 / (eight-one)
    END IF

    ! Code for 0 <= x <= 4
    IF (x <= four) THEN
      xpowth = x ** (one/three)
      IF (x < xlow) THEN
        fn_val = conlow * xpowth
      ELSE
        t = (x*x/eight-half) - half
        cheb1 = cheval(nterm1,async1,t)
        cheb2 = cheval(nterm2,async2,t)
        t = xpowth * cheb1 - (xpowth**11) * cheb2
        fn_val = t - pibrt3 * x
      END IF
    ELSE
      IF (x > xhigh1) THEN
        fn_val = zero
      ELSE
        t = (twelve-x) / (x+four)
        cheb1 = cheval(nterm3,asynca,t)
        t = lnrtp2 - x + LOG(SQRT(x)*cheb1)
        IF (t < xhigh2) THEN
          fn_val = zero
        ELSE
          fn_val = EXP(t)
        END IF
      END IF
    END IF
    RETURN
  END FUNCTION synch1


  FUNCTION synch2(xvalue) RESULT(fn_val)
    !! CALGO 757 Synchrotron radiation function
    !! \(f_2(x) = x K_{2/3}(x)\)
    !
    !   DESCRIPTION:
    !      This function calculates the synchrotron radiation function
    !      defined as
    !
    !         SYNCH2(x) = x * K(2/3)(x)
    !
    !      where K(2/3) is a modified Bessel function of order 2/3.
    !
    !      The code uses Chebyshev expansions, the coefficients of which
    !      are given to 20 decimal places.
    !
    !   ERROR RETURNS:
    !      The function is undefined if x < 0.0. If XVALUE < 0.0,
    !      an error message is printed and the function returns the value 0.0.
    !
    !   MACHINE-DEPENDENT CONSTANTS:
    !      NTERM1 - INTEGER - The no. of terms needed from the array
    !                         ASYNC1. The recommended value is such that
    !                            ABS(ASYN21(NTERM1)) < EPS/100.
    !      NTERM2 - INTEGER - The no. of terms needed from the array
    !                         ASYNC2. The recommended value is such that
    !                            ABS(ASYN22(NTERM2)) < EPS/100.
    !      NTERM3 - INTEGER - The no. of terms needed from the array
    !                         ASYNCA. The recommended value is such that
    !                            ABS(ASYN2A(NTERM3)) < EPS/100.
    !      XLOW - REAL(wp) - The value below which
    !                        SYNCH2(x) = 1.074764... * (x**(1/3))
    !                    to machine precision. The recommended value
    !                    is     sqrt (8*EPSNEG)
    !      XHIGH1 - REAL(wp) - The value above which
    !                          SYNCH2(x) = 0.0
    !                      to machine precision. The recommended value
    !                      is     -8*LN(XMIN)/7
    !      XHIGH2 - REAL(wp) - The value of LN(XMIN). This is used
    !                      to prevent underflow in calculations
    !                      for large x.
    !     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
    !     The machine-dependent constants are computed internally by
    !     using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !      EXP , LOG , SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: nterm1, nterm2, nterm3
    REAL(wp)  :: cheb1, cheb2, t, x, xhigh1, xhigh2, xlow, xpowth
    CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
    CHARACTER (LEN= 6)  :: fnname = 'SYNCH2'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                             two = 2.0_wp, three = 3.0_wp, four = 4.0_wp,  &
                             eight = 8.0_wp, ten = 10.0_wp, onehun = 100.0_wp,  &
                             conlow = 1.07476412076723931836_wp,  &
                             lnrtp2 = 0.22579135264472743236_wp

    REAL(wp), PARAMETER  :: asyn21(0:14) = [ &
        38.61783992384308548014_wp, 23.03771559496373459697_wp,  &
        5.38024998683357059676_wp, 0.61567938069957107760_wp,  &
        0.4066880046688955843e-1_wp, 0.172962745526484141e-2_wp,  &
        0.5106125883657699e-4_wp, 0.110459595022012e-5_wp,  &
        0.1823553020649e-7_wp, 0.23707698034e-9_wp, 0.248872963e-11_wp,  &
        0.2152868e-13_wp, 0.15607e-15_wp, 0.96e-18_wp, 0.1e-19_wp]

    REAL(wp), PARAMETER  :: asyn22(0:13) = [ &
        7.90631482706608042875_wp, 3.13534636128534256841_wp,  &
        0.48548794774537145380_wp, 0.3948166758272372337e-1_wp,  &
        0.196616223348088022e-2_wp, 0.6590789322930420e-4_wp,  &
        0.158575613498559e-5_wp, 0.2868653011233e-7_wp, 0.40412023595e-9_wp,  &
        0.455684443e-11_wp, 0.4204590e-13_wp, 0.32326e-15_wp, 0.210e-17_wp, 0.1e-19_wp]

    REAL(wp), PARAMETER  :: asyn2a(0:18) = [ &
        2.02033709417071360032_wp, 0.1095623712180740443e-1_wp,  &
        0.85423847301146755e-3_wp, 0.7234302421328222e-4_wp,  &
        0.631244279626992e-5_wp, 0.56481931411744e-6_wp, 0.5128324801375e-7_wp,  &
        0.471965329145e-8_wp, 0.43807442143e-9_wp, 0.4102681493e-10_wp,  &
        0.386230721e-11_wp, 0.36613228e-12_wp, 0.3480232e-13_wp, 0.333010e-14_wp,  &
        0.31856e-15_wp, 0.3074e-16_wp, 0.295e-17_wp, 0.29e-18_wp, 0.3e-19_wp]

    ! Start calculation
    x = xvalue
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    cheb1 = EPSILON(zero)
    t = cheb1 / onehun
    IF (x <= four) THEN
      DO  nterm1 = 14, 0, -1
        IF (ABS(asyn21(nterm1)) > t) EXIT
      END DO
      DO  nterm2 = 13, 0, -1
        IF (ABS(asyn22(nterm2)) > t) EXIT
      END DO
      xlow = SQRT(eight*cheb1)
    ELSE
      DO  nterm3 = 18, 0, -1
        IF (ABS(asyn2a(nterm3)) > t) EXIT
      END DO
      xhigh2 = LOG(TINY(zero))
      xhigh1 = -eight * xhigh2 / (eight-one)
    END IF

    ! Code for 0 <= x <= 4
    IF (x <= four) THEN
      xpowth = x ** (one/three)
      IF (x < xlow) THEN
        fn_val = conlow * xpowth
      ELSE
        t = (x*x/eight-half) - half
        cheb1 = cheval(nterm1,asyn21,t)
        cheb2 = cheval(nterm2,asyn22,t)
        fn_val = xpowth * cheb1 - (xpowth**5) * cheb2
      END IF
    ELSE
      IF (x > xhigh1) THEN
        fn_val = zero
      ELSE
        t = (ten-x) / (x+two)
        cheb1 = cheval(nterm3,asyn2a,t)
        t = lnrtp2 - x + LOG(SQRT(x)*cheb1)
        IF (t < xhigh2) THEN
          fn_val = zero
        ELSE
          fn_val = EXP(t)
        END IF
      END IF
    END IF
    RETURN
  END FUNCTION synch2


  FUNCTION tran02(xvalue) RESULT(fn_val)
    !! CALGO 757 Transport integral \(J_2(x) =
    !! \int_{0}^{x} \frac{t^2 e^t}{(e^t-1)^2} \, dt\)
    !
    !  DESCRIPTION:
    !    This program calculates the transport integral of order 2, defined as
    !
    !      TRAN02(X) = integral 0 to X { t**2 exp(t)/[exp(t)-1]**2 } dt
    !
    !    The program uses a Chebyshev series, the coefficients of which are
    !    given to an accuracy of 20 decimal places.
    !
    !  ERROR RETURNS:
    !    If XVALUE < 0.0, an error message is printed, and the program
    !    returns the value 0.0.
    !
    !  MACHINE-DEPENDENT CONSTANTS:
    !    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
    !                       The recommended value is such that
    !                             ATRAN(NTERMS) < EPS/100
    !    XLOW1 - REAL(wp) - The value below which TRAN02 = x to
    !                   machine precision. The recommended value is
    !                             sqrt(8*EPSNEG)
    !    XHIGH1 - REAL(wp) - The value above which the exponential series for
    !                    large x contains only one term. The recommended value
    !                    is        - ln(EPS).
    !    XHIGH2 - REAL(wp) - The value above which
    !                       TRAN02 = VALINF  -  x**2 exp(-x)
    !                    The recommended value is 2/EPS
    !    XHIGH3 - REAL(wp) - The value of ln(EPSNEG). Used to prevent overflow
    !                    for large x.
    !    For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
    !    The machine-dependent constants are computed internally by
    !    using the D1MACH subroutine.
    !
    !  INTRINSIC FUNCTIONS USED:
    !     EXP, INT, LOG, SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: k1, k2, nterms, numexp, numjn = 2
    REAL(wp)  :: rk, sumexp, sum2, t, x, xhigh1, xhigh2, xhigh3, xk, xk1, xlow1
    CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
    CHARACTER (LEN= 6)  :: fnname = 'TRAN02'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                             four = 4.0_wp, eight = 8.0_wp, onehun = 100.0_wp,  &
                             rnumjn = 2.0_wp, valinf = 3.2898681336964528729_wp

    REAL(wp), PARAMETER  :: atran(0:19) = [ &
        1.67176044643453850301_wp, -0.14773535994679448986_wp,  &
        0.1482138199469363384e-1_wp, -0.141953303263056126e-2_wp,  &
        0.13065413244157083e-3_wp, -0.1171557958675790e-4_wp,  &
        0.103334984457557e-5_wp, -0.9019113042227e-7_wp, 0.781771698331e-8_wp,  &
       -0.67445656840e-9_wp, 0.5799463945e-10_wp, -0.497476185e-11_wp,  &
        0.42596097e-12_wp, -0.3642189e-13_wp, 0.311086e-14_wp, -0.26547e-15_wp,  &
        0.2264e-16_wp, -0.193e-17_wp, 0.16e-18_wp, -0.1e-19_wp]

    ! Start execution
    x = xvalue

    ! Error test
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    xk = EPSILON(zero)
    t = xk / onehun
    IF (x <= four) THEN
      DO  nterms = 19, 0, -1
        IF (ABS(atran(nterms)) > t) EXIT
      END DO
      xlow1 = SQRT(eight*xk)
    ELSE
      xhigh1 = -LOG(2*EPSILON(zero))
      xhigh2 = one / (half*xk)
      xhigh3 = LOG(xk)
    END IF

    ! Code for x < =  4.0
    IF (x <= four) THEN
      IF (x < xlow1) THEN
        fn_val = (x**(numjn-1)) / (rnumjn-one)
      ELSE
        t = (((x*x)/eight)-half) - half
        fn_val = (x**(numjn-1)) * cheval(nterms,atran,t)
      END IF
    ELSE
      ! Code for x > 4.0
      IF (x > xhigh2) THEN
        sumexp = one
      ELSE
        IF (x <= xhigh1) THEN
          numexp = INT(xhigh1/x) + 1
          t = EXP(-x)
        ELSE
          numexp = 1
          t = one
        END IF
        rk = zero
        DO  k1 = 1, numexp
          rk = rk + one
        END DO
        sumexp = zero
        DO  k1 = 1, numexp
          sum2 = one
          xk = one / (rk*x)
          xk1 = one
          DO  k2 = 1, numjn
            sum2 = sum2 * xk1 * xk + one
            xk1 = xk1 + one
          END DO
          sumexp = sumexp * t + sum2
          rk = rk - one
        END DO
      END IF
      t = rnumjn * LOG(x) - x + LOG(sumexp)
      IF (t < xhigh3) THEN
        fn_val = valinf
      ELSE
        fn_val = valinf - EXP(t)
      END IF
    END IF
    RETURN
  END FUNCTION tran02


  FUNCTION tran03(xvalue) RESULT(fn_val)
    !! CALGO 757 Transport integral \(J_3(x) =
    !! \int_{0}^{x} \frac{t^3 e^t}{(e^t-1)^2} \, dt\)
    !
    !  DESCRIPTION:
    !    This program calculates the transport integral of order 3, defined as
    !
    !      TRAN03(X) = integral 0 to X { t**3 exp(t)/[exp(t)-1]**2 } dt
    !
    !    The program uses a Chebyshev series, the coefficients of which are
    !    given to an accuracy of 20 decimal places.
    !
    !  ERROR RETURNS:
    !    If XVALUE < 0.0, an error message is printed, and the program
    !    returns the value 0.0.
    !
    !  MACHINE-DEPENDENT CONSTANTS:
    !    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
    !                       The recommended value is such that
    !                             ATRAN(NTERMS) < EPS/100
    !    XLOW2 - REAL(wp) - The value below which TRAN03 = 0.0 to machine
    !                    precision. The recommended value is
    !                          square root of (2*XMIN)
    !    XLOW1 - REAL(wp) - The value below which TRAN03 = X**2/2 to
    !                   machine precision. The recommended value is
    !                             sqrt(8*EPSNEG)
    !    XHIGH1 - REAL(wp) - The value above which the exponential series for
    !                    large X contains only one term. The recommended value
    !                    is        - ln(EPS).
    !    XHIGH2 - REAL(wp) - The value above which
    !                       TRAN03 = VALINF  -  X**3 exp(-X)
    !                    The recommended value is 3/EPS
    !    XHIGH3 - REAL(wp) - The value of ln(EPSNEG). Used to prevent overflow
    !                    for large x.
    !    For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
    !    The machine-dependent constants are computed internally by
    !    using the D1MACH subroutine.
    !
    !  INTRINSIC FUNCTIONS USED:
    !     EXP, INT, LOG, SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: k1, k2, nterms, numexp, numjn = 3
    REAL(wp)  :: rk, sumexp, sum2, t, x, xhigh1, xhigh2, xhigh3, xk, xk1,  &
                  xlow1, xlow2
    CHARACTER (LEN=14) :: errmsg = 'ARGUMENT < 0.0'
    CHARACTER (LEN= 6) :: fnname = 'TRAN03'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                             four = 4.0_wp, eight = 8.0_wp, onehun = 100.0_wp,  &
                             rnumjn = 3.0_wp, valinf = 7.2123414189575657124_wp

    REAL(wp), PARAMETER  :: atran(0:19) = [ &
        0.76201254324387200657_wp, -0.10567438770505853250_wp,  &
        0.1197780848196578097e-1_wp, -0.121440152036983073e-2_wp,  &
        0.11550997693928547e-3_wp, -0.1058159921244229e-4_wp, 0.94746633853018e-6_wp,  &
       -0.8362212128581e-7_wp, 0.731090992775e-8_wp, -0.63505947788e-9_wp,  &
        0.5491182819e-10_wp, -0.473213954e-11_wp, 0.40676948e-12_wp, -0.3489706e-13_wp,  &
        0.298923e-14_wp, -0.25574e-15_wp, 0.2186e-16_wp, -0.187e-17_wp, 0.16e-18_wp, -0.1e-19_wp]

    ! Start execution
    x = xvalue

    ! Error test
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    xk = EPSILON(zero)
    t = xk / onehun
    IF (x <= four) THEN
      DO  nterms = 19, 0, -1
        IF (ABS(atran(nterms)) > t) EXIT
      END DO
      xlow1 = SQRT(eight*xk)
      xlow2 = SQRT(TINY(zero)/half)
    ELSE
      xhigh1 = -LOG(2*EPSILON(zero))
      xhigh2 = rnumjn / xk
      xhigh3 = LOG(xk)
    END IF

    ! Code for x < =  4.0
    IF (x <= four) THEN
      IF (x < xlow2) THEN
        fn_val = zero
      ELSE
        IF (x < xlow1) THEN
          fn_val = (x**(numjn-1)) / (rnumjn-one)
        ELSE
          t = (((x*x)/eight)-half) - half
          fn_val = (x**(numjn-1)) * cheval(nterms,atran,t)
        END IF
      END IF
    ELSE
      ! Code for x > 4.0
      IF (x > xhigh2) THEN
        sumexp = one
      ELSE
        IF (x <= xhigh1) THEN
          numexp = INT(xhigh1/x) + 1
          t = EXP(-x)
        ELSE
          numexp = 1
          t = one
        END IF
        rk = zero
        DO  k1 = 1, numexp
          rk = rk + one
        END DO
        sumexp = zero
        DO  k1 = 1, numexp
          sum2 = one
          xk = one / (rk*x)
          xk1 = one
          DO  k2 = 1, numjn
            sum2 = sum2 * xk1 * xk + one
            xk1 = xk1 + one
          END DO
          sumexp = sumexp * t + sum2
          rk = rk - one
        END DO
      END IF
      t = rnumjn * LOG(x) - x + LOG(sumexp)
      IF (t < xhigh3) THEN
        fn_val = valinf
      ELSE
        fn_val = valinf - EXP(t)
      END IF
    END IF
    RETURN
  END FUNCTION tran03


  FUNCTION tran04(xvalue) RESULT(fn_val)
    !! CALGO 757 Transport integral \(J_4(x) =
    !! \int_{0}^{x} \frac{t^4 e^t}{(e^t-1)^2} \, dt\)
    !
    !  DESCRIPTION:
    !    This program calculates the transport integral of order 4, defined as
    !
    !      TRAN04(X) = integral 0 to X { t**4 exp(t)/[exp(t)-1]**2 } dt
    !
    !    The program uses a Chebyshev series, the coefficients of which are
    !    given to an accuracy of 20 decimal places.
    !
    !  ERROR RETURNS:
    !    If XVALUE < 0.0, an error message is printed, and the program
    !    returns the value 0.0.
    !
    !  MACHINE-DEPENDENT CONSTANTS:
    !    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
    !                       The recommended value is such that
    !                             ATRAN(NTERMS) < EPS/100
    !    XLOW2 - REAL(wp) - The value below which TRAN04 = 0.0 to machine
    !                   precision. The recommended value is
    !                          cube root of (3*XMIN)
    !    XLOW1 - REAL(wp) - The value below which TRAN04 = X**3/3 to
    !                   machine precision. The recommended value is
    !                             sqrt(8*EPSNEG)
    !    XHIGH1 - REAL(wp) - The value above which the exponential series for
    !                    large X contains only one term. The recommended value
    !                    is        - ln(EPS).
    !    XHIGH2 - REAL(wp) - The value above which
    !                       TRAN04 = VALINF  -  X**4 exp(-X)
    !                    The recommended value is 4/EPS
    !    XHIGH3 - REAL(wp) - The value of ln(EPSNEG). Used to prevent overflow
    !                    for large x.
    !    For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
    !    The machine-dependent constants are computed internally by
    !    using the D1MACH subroutine.
    !
    !  INTRINSIC FUNCTIONS USED:
    !     EXP, INT, LOG, SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: k1, k2, nterms, numexp, numjn = 4
    REAL(wp)  :: rk, sumexp, sum2, t, x, xhigh1, xhigh2, xhigh3, xk, xk1,  &
                  xlow1, xlow2
    CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
    CHARACTER (LEN= 6)  :: fnname = 'TRAN04'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                             four = 4.0_wp, eight = 8.0_wp, onehun = 100.0_wp,  &
                             rnumjn = 4.0_wp, valinf = 25.975757609067316596_wp

    REAL(wp), PARAMETER  :: atran(0:19) = [ &
        0.48075709946151105786_wp, -0.8175378810321083956e-1_wp,  &
        0.1002700665975162973e-1_wp, -0.105993393598201507e-2_wp,  &
        0.10345062450304053e-3_wp, -0.964427054858991e-5_wp, 0.87455444085147e-6_wp,  &
       -0.7793212079811e-7_wp, 0.686498861410e-8_wp, -0.59995710764e-9_wp,  &
        0.5213662413e-10_wp, -0.451183819e-11_wp, 0.38921592e-12_wp, -0.3349360e-13_wp,  &
        0.287667e-14_wp, -0.24668e-15_wp, 0.2113e-16_wp, -0.181e-17_wp, 0.15e-18_wp, -0.1e-19_wp]

    ! Start execution
    x = xvalue

    ! Error test
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    xk = EPSILON(zero)
    t = xk / onehun
    IF (x <= four) THEN
      DO  nterms = 19, 0, -1
        IF (ABS(atran(nterms)) > t) EXIT
      END DO
      xlow1 = SQRT(eight*xk)
      xk1 = rnumjn - one
      xlow2 = (xk1*TINY(zero)) ** (one/xk1)
    ELSE
      xhigh1 = -LOG(2*EPSILON(zero))
      xhigh2 = rnumjn / xk
      xhigh3 = LOG(xk)
    END IF

    ! Code for x < =  4.0
    IF (x <= four) THEN
      IF (x < xlow2) THEN
        fn_val = zero
      ELSE
        IF (x < xlow1) THEN
          fn_val = (x**(numjn-1)) / (rnumjn-one)
        ELSE
          t = (((x*x)/eight)-half) - half
          fn_val = (x**(numjn-1)) * cheval(nterms,atran,t)
        END IF
      END IF
    ELSE
      ! Code for x > 4.0
      IF (x > xhigh2) THEN
        sumexp = one
      ELSE
        IF (x <= xhigh1) THEN
          numexp = INT(xhigh1/x) + 1
          t = EXP(-x)
        ELSE
          numexp = 1
          t = one
        END IF
        rk = zero
        DO  k1 = 1, numexp
          rk = rk + one
        END DO
        sumexp = zero
        DO  k1 = 1, numexp
          sum2 = one
          xk = one / (rk*x)
          xk1 = one
          DO  k2 = 1, numjn
            sum2 = sum2 * xk1 * xk + one
            xk1 = xk1 + one
          END DO
          sumexp = sumexp * t + sum2
          rk = rk - one
        END DO
      END IF
      t = rnumjn * LOG(x) - x + LOG(sumexp)
      IF (t < xhigh3) THEN
        fn_val = valinf
      ELSE
        fn_val = valinf - EXP(t)
      END IF
    END IF
    RETURN
  END FUNCTION tran04


  FUNCTION tran05(xvalue) RESULT(fn_val)
    !! CALGO 757 Transport integral \(J_5(x) =
    !! \int_{0}^{x} \frac{t^5 e^t}{(e^t-1)^2} \, dt\)
    !
    !  DESCRIPTION:
    !    This program calculates the transport integral of order 5, defined as
    !
    !      TRAN05(X) = integral 0 to X { t**5 exp(t)/[exp(t)-1]**2 } dt
    !
    !    The program uses a Chebyshev series, the coefficients of which are
    !    given to an accuracy of 20 decimal places.
    !
    !  ERROR RETURNS:
    !    If XVALUE < 0.0, an error message is printed, and the program
    !    returns the value 0.0.
    !
    !  MACHINE-DEPENDENT CONSTANTS:
    !    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
    !                       The recommended value is such that
    !                             ATRAN(NTERMS) < EPS/100
    !    XLOW2 - REAL(wp) - The value below which TRAN05 = 0.0 to machine
    !                   precision. The recommended value is
    !                          4th root of (4*XMIN)
    !    XLOW1 - REAL(wp) - The value below which TRAN05 = X**4/4 to
    !                   machine precision. The recommended value is
    !                             sqrt(8*EPSNEG)
    !    XHIGH1 - REAL(wp) - The value above which the exponential series for
    !                    large X contains only one term. The recommended value
    !                    is        - ln(EPS).
    !    XHIGH2 - REAL(wp) - The value above which
    !                       TRAN05 = VALINF  -  X**5 exp(-X)
    !                    The recommended value is 5/EPS
    !    XHIGH3 - REAL(wp) - The value of ln(EPSNEG). Used to prevent overflow
    !                    for large x.
    !     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
    !     The machine-dependent constants are computed internally by
    !     using the D1MACH subroutine.
    !
    !  INTRINSIC FUNCTIONS USED:
    !     EXP, INT, LOG, SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val
    INTEGER    :: k1, k2, nterms, numexp, numjn = 5
    REAL(wp)  :: rk, sumexp, sum2, t, x, xhigh1, xhigh2, xhigh3, xk, xk1,  &
                  xlow1, xlow2
    CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
    CHARACTER (LEN= 6)  :: fnname = 'TRAN05'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                             four = 4.0_wp, eight = 8.0_wp, onehun = 100.0_wp,  &
                             rnumjn = 5.0_wp,  &
                             valinf = 0.12443133061720439116D3
    REAL(wp), PARAMETER  :: atran(0:19) = [ &
        0.34777777713391078928_wp, -0.6645698897605042801e-1_wp,  &
        0.861107265688330882e-2_wp, -0.93966822237555384e-3_wp,  &
        0.9363248060815134e-4_wp, -0.885713193408328e-5_wp, 0.81191498914503e-6_wp,  &
       -0.7295765423277e-7_wp, 0.646971455045e-8_wp, -0.56849028255e-9_wp,  &
        0.4962559787e-10_wp, -0.431093996e-11_wp, 0.37310094e-12_wp, -0.3219769e-13_wp,  &
        0.277220e-14_wp, -0.23824e-15_wp, 0.2044e-16_wp, -0.175e-17_wp, 0.15e-18_wp, -0.1e-19_wp]

    ! Start execution
    x = xvalue

    ! Error test
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    xk = EPSILON(zero)
    t = xk / onehun
    IF (x <= four) THEN
      DO  nterms = 19, 0, -1
        IF (ABS(atran(nterms)) > t) EXIT
      END DO
      xlow1 = SQRT(eight*xk)
      xk1 = rnumjn - one
      xlow2 = (xk1*TINY(zero)) ** (one/xk1)
    ELSE
      xhigh1 = -LOG(2*EPSILON(zero))
      xhigh2 = rnumjn / xk
      xhigh3 = LOG(xk)
    END IF

    ! Code for x < =  4.0
    IF (x <= four) THEN
      IF (x < xlow2) THEN
        fn_val = zero
      ELSE
        IF (x < xlow1) THEN
          fn_val = (x**(numjn-1)) / (rnumjn-one)
        ELSE
          t = (((x*x)/eight)-half) - half
          fn_val = (x**(numjn-1)) * cheval(nterms,atran,t)
        END IF
      END IF
    ELSE
      ! Code for x > 4.0
      IF (x > xhigh2) THEN
        sumexp = one
      ELSE
        IF (x <= xhigh1) THEN
          numexp = INT(xhigh1/x) + 1
          t = EXP(-x)
        ELSE
          numexp = 1
          t = one
        END IF
        rk = zero
        DO  k1 = 1, numexp
          rk = rk + one
        END DO
        sumexp = zero
        DO  k1 = 1, numexp
          sum2 = one
          xk = one / (rk*x)
          xk1 = one
          DO  k2 = 1, numjn
            sum2 = sum2 * xk1 * xk + one
            xk1 = xk1 + one
          END DO
          sumexp = sumexp * t + sum2
          rk = rk - one
        END DO
      END IF
      t = rnumjn * LOG(x) - x + LOG(sumexp)
      IF (t < xhigh3) THEN
        fn_val = valinf
      ELSE
        fn_val = valinf - EXP(t)
      END IF
    END IF
    RETURN
  END FUNCTION tran05


  FUNCTION tran06(xvalue) RESULT(fn_val)
    !! CALGO 757 Transport integral \(J_6(x) =
    !! \int_{0}^{x} \frac{t^6 e^t}{(e^t-1)^2} \, dt\)
    !
    !  DESCRIPTION:
    !    This program calculates the transport integral of order 6, defined as
    !
    !      TRAN06(X) = integral 0 to X { t**6 exp(t)/[exp(t)-1]**2 } dt
    !
    !    The program uses a Chebyshev series, the coefficients of which are
    !    given to an accuracy of 20 decimal places.
    !
    !  ERROR RETURNS:
    !    If XVALUE < 0.0, an error message is printed, and the program
    !    returns the value 0.0.
    !
    !  MACHINE-DEPENDENT CONSTANTS:
    !    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
    !                       The recommended value is such that
    !                             ATRAN(NTERMS) < EPS/100
    !    XLOW2 - REAL(wp) - The value below which TRAN06 = 0.0 to machine
    !                   precision. The recommended value is
    !                          5th root of (5*XMIN)
    !    XLOW1 - REAL(wp) - The value below which TRAN06 = X**5/5 to
    !                   machine precision. The recommended value is
    !                             sqrt(8*EPSNEG)
    !    XHIGH1 - REAL(wp) - The value above which the exponential series for
    !                    large X contains only one term. The recommended value
    !                    is        - ln(EPS).
    !    XHIGH2 - REAL(wp) - The value above which
    !                       TRAN06 = VALINF  -  X**6 exp(-X)
    !                    The recommended value is 6/EPS
    !    XHIGH3 - REAL(wp) - The value of ln(EPSNEG). Used to prevent overflow
    !                    for large x.
    !     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
    !     The machine-dependent constants are computed internally by
    !     using the D1MACH subroutine.
    !
    !  INTRINSIC FUNCTIONS USED:
    !     EXP, INT, LOG, SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val
    INTEGER    :: k1, k2, nterms, numexp, numjn = 6
    REAL(wp)  :: rk, sumexp, sum2, t, x, xhigh1, xhigh2, xhigh3, xk, xk1,  &
                  xlow1, xlow2
    CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
    CHARACTER (LEN= 6)  :: fnname = 'TRAN06'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                             four = 4.0_wp, eight = 8.0_wp, onehun = 100.0_wp,  &
                             rnumjn = 6.0_wp, valinf = 732.48700462880338059_wp

    REAL(wp), PARAMETER  :: atran(0:19) = [ &
        0.27127335397840008227_wp, -0.5588610553191453393e-1_wp,  &
        0.753919513290083056e-2_wp, -0.84351138579211219e-3_wp, 0.8549098079676702e-4_wp,  &
       -0.818715493293098e-5_wp, 0.75754240427986e-6_wp, -0.6857306541831e-7_wp,  &
        0.611700376031e-8_wp, -0.54012707024e-9_wp, 0.4734306435e-10_wp, -0.412701055e-11_wp, &
        0.35825603e-12_wp, -0.3099752e-13_wp, 0.267501e-14_wp, -0.23036e-15_wp,  &
        0.1980e-16_wp, -0.170e-17_wp, 0.15e-18_wp, -0.1e-19_wp]

    ! Start execution
    x = xvalue

    ! Error test
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    xk = EPSILON(zero)
    t = xk / onehun
    IF (x <= four) THEN
      DO  nterms = 19, 0, -1
        IF (ABS(atran(nterms)) > t) EXIT
      END DO
      xlow1 = SQRT(eight*xk)
      xk1 = rnumjn - one
      xlow2 = (xk1*TINY(zero)) ** (one/xk1)
    ELSE
      xhigh1 = -LOG(2*EPSILON(zero))
      xhigh2 = rnumjn / xk
      xhigh3 = LOG(xk)
    END IF

    ! Code for x < =  4 .0
    IF (x <= four) THEN
      IF (x < xlow2) THEN
        fn_val = zero
      ELSE
        IF (x < xlow1) THEN
          fn_val = (x**(numjn-1)) / (rnumjn-one)
        ELSE
          t = (((x*x)/eight)-half) - half
          fn_val = (x**(numjn-1)) * cheval(nterms,atran,t)
        END IF
      END IF
    ELSE
      ! Code for x > 4 .0
      IF (x > xhigh2) THEN
        sumexp = one
      ELSE
        IF (x <= xhigh1) THEN
          numexp = INT(xhigh1/x) + 1
          t = EXP(-x)
        ELSE
          numexp = 1
          t = one
        END IF
        rk = zero
        DO  k1 = 1, numexp
          rk = rk + one
        END DO
        sumexp = zero
        DO  k1 = 1, numexp
          sum2 = one
          xk = one / (rk*x)
          xk1 = one
          DO  k2 = 1, numjn
            sum2 = sum2 * xk1 * xk + one
            xk1 = xk1 + one
          END DO
          sumexp = sumexp * t + sum2
          rk = rk - one
        END DO
      END IF
      t = rnumjn * LOG(x) - x + LOG(sumexp)
      IF (t < xhigh3) THEN
        fn_val = valinf
      ELSE
        fn_val = valinf - EXP(t)
      END IF
    END IF
    RETURN
  END FUNCTION tran06


  FUNCTION tran07(xvalue) RESULT(fn_val)
    !! CALGO 757 Transport integral \(J_7(x) =
    !! \int_{0}^{x} \frac{t^7 e^t}{(e^t-1)^2} \, dt\)
    !
    !  DESCRIPTION:
    !    This program calculates the transport integral of order 7, defined as
    !
    !      TRAN07(X) = integral 0 to X { t**7 exp(t)/[exp(t)-1]**2 } dt
    !
    !    The program uses a Chebyshev series, the coefficients of which are
    !    given to an accuracy of 20 decimal places.
    !
    !  ERROR RETURNS:
    !    If XVALUE < 0.0, an error message is printed, and the program
    !    returns the value 0.0.
    !
    !  MACHINE-DEPENDENT CONSTANTS:
    !    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
    !                       The recommended value is such that
    !                             ATRAN(NTERMS) < EPS/100
    !    XLOW2 - REAL(wp) - The value below which TRAN07 = 0.0 to machine
    !                   precision. The recommended value is
    !                          6th root of (6*XMIN)
    !    XLOW1 - REAL(wp) - The value below which TRAN07 = X**6/6 to
    !                   machine precision. The recommended value is
    !                             sqrt(8*EPSNEG)
    !    XHIGH1 - REAL(wp) - The value above which the exponential series for
    !                    large X contains only one term. The recommended value
    !                    is        - ln(EPS).
    !    XHIGH2 - REAL(wp) - The value above which
    !                       TRAN07 = VALINF  -  X**7 exp(-X)
    !                    The recommended value is 7/EPS
    !    XHIGH3 - REAL(wp) - The value of ln(EPSNEG). Used to prevent overflow
    !                    for large x.
    !     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
    !     The machine-dependent constants are computed internally by
    !     using the D1MACH subroutine.
    !
    !  INTRINSIC FUNCTIONS USED:
    !     EXP, INT, LOG, SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val
    INTEGER    :: k1, k2, nterms, numexp, numjn = 7
    REAL(wp)  :: rk, sumexp, sum2, t, x, xhigh1, xhigh2, xhigh3, xk, xk1,  &
                  xlow1, xlow2
    CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
    CHARACTER (LEN= 6)  :: fnname = 'TRAN07'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                             four = 4.0_wp, eight = 8.0_wp, onehun = 100.0_wp,  &
                             rnumjn = 7.0_wp,  &
                             valinf = 0.50820803580048910473D4
    REAL(wp), PARAMETER  :: atran(0:19) = [ &
        0.22189250734010404423_wp, -0.4816751061177993694e-1_wp,  &
        0.670092448103153629e-2_wp, -0.76495183443082557e-3_wp,  &
        0.7863485592348690e-4_wp, -0.761025180887504e-5_wp,  &
        0.70991696299917e-6_wp, -0.6468025624903e-7_wp, 0.580039233960e-8_wp,  &
        -0.51443370149e-9_wp, 0.4525944183e-10_wp, -0.395800363e-11_wp,  &
        0.34453785e-12_wp, -0.2988292e-13_wp, 0.258434e-14_wp, -0.22297e-15_wp,  &
        0.1920e-16_wp, -0.165e-17_wp, 0.14e-18_wp, -0.1e-19_wp]

    ! Start execution
    x = xvalue

    ! Error test
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    xk = EPSILON(zero)
    t = xk / onehun
    IF (x <= four) THEN
      DO  nterms = 19, 0, -1
        IF (ABS(atran(nterms)) > t) EXIT
      END DO
      xlow1 = SQRT(eight*xk)
      xk1 = rnumjn - one
      xlow2 = (xk1*TINY(zero)) ** (one/xk1)
    ELSE
      xhigh1 = -LOG(2*EPSILON(zero))
      xhigh2 = rnumjn / xk
      xhigh3 = LOG(xk)
    END IF

    ! Code for x <= 4.0
    IF (x <= four) THEN
      IF (x < xlow2) THEN
        fn_val = zero
      ELSE
        IF (x < xlow1) THEN
          fn_val = (x**(numjn-1)) / (rnumjn-one)
        ELSE
          t = (((x*x)/eight)-half) - half
          fn_val = (x**(numjn-1)) * cheval(nterms,atran,t)
        END IF
      END IF
    ELSE
      ! Code for x > 4.0
      IF (x > xhigh2) THEN
        sumexp = one
      ELSE
        IF (x <= xhigh1) THEN
          numexp = INT(xhigh1/x) + 1
          t = EXP(-x)
        ELSE
          numexp = 1
          t = one
        END IF
        rk = zero
        DO  k1 = 1, numexp
          rk = rk + one
        END DO
        sumexp = zero
        DO  k1 = 1, numexp
          sum2 = one
          xk = one / (rk*x)
          xk1 = one
          DO  k2 = 1, numjn
            sum2 = sum2 * xk1 * xk + one
            xk1 = xk1 + one
          END DO
          sumexp = sumexp * t + sum2
          rk = rk - one
        END DO
      END IF
      t = rnumjn * LOG(x) - x + LOG(sumexp)
      IF (t < xhigh3) THEN
        fn_val = valinf
      ELSE
        fn_val = valinf - EXP(t)
      END IF
    END IF
    RETURN
  END FUNCTION tran07


  FUNCTION tran08(xvalue) RESULT(fn_val)
    !! CALGO 757 Transport integral \(J_8(x) =
    !! \int_{0}^{x} \frac{t^8 e^t}{(e^t-1)^2} \, dt\)
    !
    !  DESCRIPTION:
    !    This program calculates the transport integral of order 8, defined as
    !
    !      TRAN08(X) = integral 0 to X { t**8 exp(t)/[exp(t)-1]**2 } dt
    !
    !    The program uses a Chebyshev series, the coefficients of which are
    !    given to an accuracy of 20 decimal places.
    !
    !  ERROR RETURNS:
    !    If XVALUE < 0.0, an error message is printed, and the program
    !    returns the value 0.0.
    !
    !  MACHINE-DEPENDENT CONSTANTS:
    !    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
    !                       The recommended value is such that
    !                             ATRAN(NTERMS) < EPS/100
    !    XLOW2 - REAL(wp) - The value below which TRAN08 = 0.0 to machine
    !                   precision. The recommended value is
    !                          7th root of (7*XMIN)
    !    XLOW1 - REAL(wp) - The value below which TRAN08 = X**7/7 to
    !                   machine precision. The recommended value is
    !                             sqrt(8*EPSNEG)
    !    XHIGH1 - REAL(wp) - The value above which the exponential series for
    !                    large X contains only one term. The recommended value
    !                    is        - ln(EPS).
    !    XHIGH2 - REAL(wp) - The value above which
    !                       TRAN08 = VALINF  -  X**8 exp(-X)
    !                    The recommended value is 8/EPS
    !    XHIGH3 - REAL(wp) - The value of ln(EPSNEG). Used to prevent overflow
    !                    for large x.
    !     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
    !     The machine-dependent constants are computed internally by
    !     using the D1MACH subroutine.
    !
    ! INTRINSIC FUNCTIONS USED:
    !     EXP, INT, LOG, SQRT
    !
    ! OTHER MISCFUN SUBROUTINES USED:
    !     CHEVAL, ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: k1, k2, nterms, numexp, numjn = 8
    REAL(wp)  :: rk, sumexp, sum2, t, x, xhigh1, xhigh2, xhigh3, xk, xk1,  &
                  xlow1, xlow2
    CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
    CHARACTER (LEN= 6)  :: fnname = 'TRAN08'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                             four = 4.0_wp, eight = 8.0_wp, onehun = 100.0_wp,  &
                             rnumjn = 8.0_wp, valinf = 0.40484399001901115764D5
    REAL(wp), PARAMETER  :: atran(0:19) = [ &
        0.18750695774043719233_wp, -0.4229527646093673337e-1_wp,  &
        0.602814856929065592e-2_wp, -0.69961054811814776e-3_wp,  &
        0.7278482421298789e-4_wp, -0.710846250050067e-5_wp, 0.66786706890115e-6_wp,  &
       -0.6120157501844e-7_wp, 0.551465264474e-8_wp, -0.49105307052e-9_wp,  &
        0.4335000869e-10_wp, -0.380218700e-11_wp, 0.33182369e-12_wp, -0.2884512e-13_wp,  &
        0.249958e-14_wp, -0.21605e-15_wp, 0.1863e-16_wp, -0.160e-17_wp, 0.14e-18_wp, -0.1e-19_wp]

    ! Start execution
    x = xvalue

    ! Error test
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    xk = EPSILON(zero)
    t = xk / onehun
    IF (x <= four) THEN
      DO  nterms = 19, 0, -1
        IF (ABS(atran(nterms)) > t) EXIT
      END DO
      xlow1 = SQRT(eight*xk)
      xk1 = rnumjn - one
      xlow2 = (xk1*TINY(zero)) ** (one/xk1)
    ELSE
      xhigh1 = -LOG(2*EPSILON(zero))
      xhigh2 = rnumjn / xk
      xhigh3 = LOG(xk)
    END IF

    ! Code for x < =  4.0
    IF (x <= four) THEN
      IF (x < xlow2) THEN
        fn_val = zero
      ELSE
        IF (x < xlow1) THEN
          fn_val = (x**(numjn-1)) / (rnumjn-one)
        ELSE
          t = (((x*x)/eight)-half) - half
          fn_val = (x**(numjn-1)) * cheval(nterms,atran,t)
        END IF
      END IF
    ELSE
      ! Code for x > 4.0
      IF (x > xhigh2) THEN
        sumexp = one
      ELSE
        IF (x <= xhigh1) THEN
          numexp = INT(xhigh1/x) + 1
          t = EXP(-x)
        ELSE
          numexp = 1
          t = one
        END IF
        rk = zero
        DO  k1 = 1, numexp
          rk = rk + one
        END DO
        sumexp = zero
        DO  k1 = 1, numexp
          sum2 = one
          xk = one / (rk*x)
          xk1 = one
          DO  k2 = 1, numjn
            sum2 = sum2 * xk1 * xk + one
            xk1 = xk1 + one
          END DO
          sumexp = sumexp * t + sum2
          rk = rk - one
        END DO
      END IF
      t = rnumjn * LOG(x) - x + LOG(sumexp)
      IF (t < xhigh3) THEN
        fn_val = valinf
      ELSE
        fn_val = valinf - EXP(t)
      END IF
    END IF
    RETURN
  END FUNCTION tran08


  FUNCTION tran09(xvalue) RESULT(fn_val)
    !! CALGO 757 Transport integral \(J_9(x) =
    !! \int_{0}^{x} \frac{t^9 e^t}{(e^t-1)^2} \, dt\)
    !
    ! DESCRIPTION:
    !   This program calculates the transport integral of order 9, defined as
    !
    !      TRAN09(X) = integral 0 to X { t**9 exp(t)/[exp(t)-1]**2 } dt
    !
    !   The program uses a Chebyshev series, the coefficients of which are
    !   given to an accuracy of 20 decimal places.
    !
    ! ERROR RETURNS:
    !   If XVALUE < 0.0, an error message is printed, and the program
    !   returns the value 0.0.
    !
    ! MACHINE-DEPENDENT CONSTANTS:
    !   NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
    !                      The recommended value is such that
    !                           ATRAN(NTERMS) < EPS/100
    !   XLOW2 - REAL(dp) - The value below which TRAN09 = 0.0 to machine
    !                      precision. The recommended value is
    !                           8th root of (8*XMIN)
    !   XLOW1 - REAL(dp) - The value below which TRAN09 = X**8/8 to
    !                      machine precision. The recommended value is
    !                           sqrt(8*EPSNEG)
    !   XHIGH1 - REAL(dp) - The value above which the exponential series for
    !                       large X contains only one term. The recommended value
    !                       is - ln(EPS).
    !   XHIGH2 - REAL(dp) - The value above which
    !                           TRAN09 = VALINF - X**9 exp(-X)
    !                       The recommended value is 9/EPS
    !   XHIGH3 - REAL(dp) - The value of ln(EPSNEG). Used to prevent overflow
    !                       for large x.
    !   For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
    !   The machine-dependent constants are computed internally by
    !   using the D1MACH subroutine.
    !
    ! INTRINSIC FUNCTIONS USED:
    !     EXP, INT, LOG, SQRT
    !
    ! OTHER MISCFUN SUBROUTINES USED:
    !     CHEVAL, ERRPRN, D1MACH

    REAL(wp), INTENT(IN) :: xvalue
    REAL(wp)             :: fn_val
    INTEGER  :: k1, k2, nterms, numexp, numjn = 9
    REAL(wp) :: rk, sumexp, sum2, t, x, xhigh1, xhigh2, xhigh3, xk, xk1, xlow1, xlow2
    CHARACTER (LEN=14) :: errmsg = 'ARGUMENT < 0.0'
    CHARACTER (LEN= 6) :: fnname = 'TRAN09'

    REAL(wp), PARAMETER :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                           four = 4.0_wp, eight = 8.0_wp, onehun = 100.0_wp, &
                           rnumjn = 9.0_wp, valinf = 0.36360880558872871397D6

    REAL(wp), PARAMETER  :: atran(0:19) = [ &
        0.16224049991949846835_wp, -0.3768351452195937773e-1_wp,  &
        0.547669715917719770e-2_wp, -0.64443945009449521e-3_wp, 0.6773645285280983e-4_wp,  &
       -0.666813497582042e-5_wp, 0.63047560019047e-6_wp, -0.5807478663611e-7_wp,  &
        0.525551305123e-8_wp, -0.46968861761e-9_wp, 0.4159395065e-10_wp, -0.365808491e-11_wp, &
        0.32000794e-12_wp, -0.2787651e-13_wp, 0.242017e-14_wp, -0.20953e-15_wp, 0.1810e-16_wp,  &
       -0.156e-17_wp, 0.13e-18_wp, -0.1e-19_wp]

    ! Start execution
    x = xvalue

    ! Error test
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF
! Compute the machine-dependent constants.
    xk = EPSILON(zero)
    t = xk / onehun
    IF (x <= four) THEN
      DO  nterms = 19, 0, -1
        IF (ABS(atran(nterms)) > t) EXIT
      END DO
      xlow1 = SQRT(eight*xk)
      xk1 = rnumjn - one
      xlow2 = (xk1*TINY(zero)) ** (one/xk1)
    ELSE
      xhigh1 = -LOG(2*EPSILON(zero))
      xhigh2 = rnumjn / xk
      xhigh3 = LOG(xk)
    END IF

    ! Code for x < =  4.0
    IF (x <= four) THEN
      IF (x < xlow2) THEN
        fn_val = zero
      ELSE
        IF (x < xlow1) THEN
          fn_val = (x**(numjn-1)) / (rnumjn-one)
        ELSE
          t = (((x*x)/eight)-half) - half
          fn_val = (x**(numjn-1)) * cheval(nterms,atran,t)
        END IF
      END IF
    ELSE
      ! Code for x > 4.0
      IF (x > xhigh2) THEN
        sumexp = one
      ELSE
        IF (x <= xhigh1) THEN
          numexp = INT(xhigh1/x) + 1
          t = EXP(-x)
        ELSE
          numexp = 1
          t = one
        END IF
        rk = zero
        DO  k1 = 1, numexp
          rk = rk + one
        END DO
        sumexp = zero
        DO  k1 = 1, numexp
          sum2 = one
          xk = one / (rk*x)
          xk1 = one
          DO  k2 = 1, numjn
            sum2 = sum2 * xk1 * xk + one
            xk1 = xk1 + one
          END DO
          sumexp = sumexp * t + sum2
          rk = rk - one
        END DO
      END IF
      t = rnumjn * LOG(x) - x + LOG(sumexp)
      IF (t < xhigh3) THEN
        fn_val = valinf
      ELSE
        fn_val = valinf - EXP(t)
      END IF
    END IF
    RETURN
  END FUNCTION tran09


  FUNCTION atnint(xvalue) RESULT(fn_val)
    !! CALGO 757 Inverse-tangent integral \(\int_{0}^{x} \frac{\arctan t}{t} \, dt\)
    !
    ! DESCRIPTION:
    !   The function ATNINT calculates the value of the
    !   inverse-tangent integral defined by
    !
    !       ATNINT(x) = integral 0 to x ( (arctan t)/t ) dt
    !
    !   The approximation uses Chebyshev series with the coefficients
    !   given to an accuracy of 20D.
    !
    ! ERROR RETURNS:
    !   There are no error returns from this program.
    !
    ! MACHINE-DEPENDENT CONSTANTS:
    !   NTERMS - INTEGER - The no. of terms of the array ATNINTT.
    !                      The recommended value is such that
    !                          ATNINA(NTERMS) < EPS/100
    !   XLOW - REAL(wp) - A bound below which ATNINT(x) = x to machine
    !                 precision. The recommended value is
    !                     sqrt(EPSNEG/2).
    !   XUPPER - REAL(wp) - A bound on x, above which, to machine precision
    !                   ATNINT(x) = (pi/2)ln x
    !                   The recommended value is 1/EPS.
    !
    !     For values of EPSNEG and EPS for various machine/compiler
    !     combinations refer to the text file MACHCON.TXT
    !
    !     The machine-dependent constants are computed internally by
    !     using the D1MACH subroutine.
    !
    ! INTRINSIC FUNCTIONS USED:
    !    ABS , LOG
    ! 
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL ,  D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: ind, nterms
    REAL(wp)  :: t, x, xlow, xupper

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                             onehun = 100.0_wp, twobpi = 0.63661977236758134308_wp

    REAL(wp), PARAMETER  :: atnina(0:22) = [ &
        1.91040361296235937512_wp, -0.4176351437656746940e-1_wp,  &
        0.275392550786367434e-2_wp, -0.25051809526248881e-3_wp, 0.2666981285121171e-4_wp, &
       -0.311890514107001e-5_wp, 0.38833853132249e-6_wp, -0.5057274584964e-7_wp,  &
        0.681225282949e-8_wp, -0.94212561654e-9_wp, 0.13307878816e-9_wp,  &
       -0.1912678075e-10_wp, 0.278912620e-11_wp, -0.41174820e-12_wp, 0.6142987e-13_wp,  &
       -0.924929e-14_wp, 0.140387e-14_wp, -0.21460e-15_wp, 0.3301e-16_wp, -0.511e-17_wp,  &
        0.79e-18_wp, -0.12e-18_wp, 0.2e-19_wp]

    ! Compute the machine-dependent constants.
    t = 2*EPSILON(zero) / onehun
    DO  nterms = 22, 0, -1
      IF (ABS(atnina(nterms)) > t) EXIT
    END DO
    t = EPSILON(zero)
    xlow = SQRT(t/(one+one))
    xupper = one / t

    ! Start calculation
    ind = 1
    x = xvalue
    IF (x < zero) THEN
      x = -x
      ind = -1
    END IF

    ! Code for X < =  1.0
    IF (x <= one) THEN
      IF (x < xlow) THEN
        fn_val = x
      ELSE
        t = x * x
        t = (t-half) + (t-half)
        fn_val = x * cheval(nterms,atnina,t)
      END IF
    ELSE
      ! Code for X > 1.0
      IF (x > xupper) THEN
        fn_val = LOG(x) / twobpi
      ELSE
        t = one / (x*x)
        t = (t-half) + (t-half)
        fn_val = LOG(x) / twobpi + cheval(nterms,atnina,t) / x
      END IF
    END IF
    IF (ind < 0) fn_val = -fn_val
    RETURN
  END FUNCTION atnint


  FUNCTION clausn(xvalue) RESULT(fn_val)
    !! CALGO 757 Clausen's integral
    !! \(-\int_{0}^{x} \ln{\lvert 2 \sin\frac{t}{2} \rvert} \, dt\)
    !
    ! DESCRIPTION:
    !   This program calculates Clausen's integral defined by
    !
    !          CLAUSN(x) = integral 0 to x of (-ln(2*sin(t/2))) dt
    !
    !   The code uses Chebyshev expansions with the coefficients
    !   given to 20 decimal places.
    !
    ! ERROR RETURNS:
    !   If |x| is too large it is impossible to reduce the argument
    !   to the range [0,2*pi] with any precision. An error message
    !   is printed and the program returns the value 0.0
    !
    ! MACHINE-DEPENDENT CONSTANTS:
    !   NTERMS - INTEGER - the no. of terms of the array ACLAUS
    !                      to be used. The recommended value is
    !                      such that ABS(ACLAUS(NTERMS)) < EPS/100
    !                      subject to 1 <= NTERMS <= 15
    !   XSMALL - REAL(wp) - the value below which Cl(x) can be
    !                   approximated by x (1-ln x). The recommended
    !                   value is pi*sqrt(EPSNEG/2).
    !   XHIGH - REAL(wp) - The value of |x| above which we cannot
    !                  reliably reduce the argument to [0,2*pi].
    !                  The recommended value is   1/EPS.
    !
    !     For values of EPS and EPSNEG refer to the file MACHCON.TXT
    !
    !     The machine-dependent constants are computed internally by
    !     using the D1MACH subroutine.
    !
    ! INTRINSIC FUNCTIONS USED:
    !   AINT , LOG , SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: indx, nterms
    REAL(wp)  :: t, x, xhigh, xsmall
    CHARACTER (LEN=26)  :: errmsg = 'ARGUMENT TOO LARGE IN SIZE'
    CHARACTER (LEN= 6)  :: fnname = 'CLAUSN'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                             onehun = 100.0_wp, pi = 3.1415926535897932385_wp,  &
                             pisq = 9.8696044010893586188_wp,  &
                             twopi = 6.2831853071795864769_wp,  &
                             twopia = 6.28125_wp, twopib = 0.19353071795864769253e-2_wp
    REAL(wp), PARAMETER  :: aclaus(0:15) = [ &
        2.14269436376668844709_wp, 0.7233242812212579245e-1_wp,  &
        0.101642475021151164e-2_wp, 0.3245250328531645e-4_wp,  &
        0.133315187571472e-5_wp, 0.6213240591653e-7_wp, 0.313004135337e-8_wp,  &
        0.16635723056e-9_wp, 0.919659293e-11_wp, 0.52400462e-12_wp,  &
        0.3058040e-13_wp, 0.181969e-14_wp, 0.11004e-15_wp, 0.675e-17_wp, 0.42e-18_wp, 0.3e-19_wp]

    ! Start execution
    x = xvalue

    ! Compute the machine-dependent constants.
    t = EPSILON(zero)
    xhigh = one / t

    ! Error test
    IF (ABS(x) > xhigh) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Continue with machine-dependent constants
    xsmall = pi * SQRT(half*t)
    t = t / onehun
    DO  nterms = 15, 0, -1
      IF (ABS(aclaus(nterms)) > t) GO TO 20
    END DO

    ! Continue with computation
  20 indx = 1
    IF (x < zero) THEN
      x = -x
      indx = -1
    END IF

    ! Argument reduced using simulated extra precision
    IF (x > twopi) THEN
      t = AINT(x/twopi)
      x = (x-t*twopia) - t * twopib
    END IF
    IF (x > pi) THEN
      x = (twopia-x) + twopib
      indx = -indx
    END IF

    ! Set result to zero if X multiple of PI
    IF (x == zero) THEN
      fn_val = zero
      RETURN
    END IF

    ! Code for X < XSMALL
    IF (x < xsmall) THEN
      fn_val = x * (one-LOG(x))
    ELSE
      ! Code for XSMALL < =  X < =  PI
      t = (x*x) / pisq - half
      t = t + t
      IF (t > one) t = one
      fn_val = x * cheval(nterms,aclaus,t) - x * LOG(x)
    END IF
    IF (indx < 0) fn_val = -fn_val
    RETURN
  END FUNCTION clausn


  FUNCTION exp3(xvalue) RESULT(fn_val)
    !! CALGO 757 \(\int_{0}^{x} e^{-t^3} \, dt\)
    !
    !   DESCRIPTION
    !      This function calculates
    !
    !           EXP3(X) = integral 0 to X  (exp(-t*t*t)) dt
    !
    !      The code uses Chebyshev expansions, whose coefficients are
    !      given to 20 decimal places.
    !
    !   ERROR RETURNS
    !      If XVALUE < 0, an error message is printed and the function
    !      returns the value 0.
    !
    !   MACHINE-DEPENDENT CONSTANTS
    !      NTERM1 - INTEGER - The no. of terms of the array AEXP3,
    !                         The recommended value is such that
    !                               AEXP3(NTERM1) < EPS/100.
    !      NTERM2 - INTEGER - The no. of terms of the array AEXP3A.
    !                         The recommended value is such that
    !                               AEXP3A(NTERM2) < EPS/100.
    !      XLOW - REAL(wp) - The value below which EXP3(X) = X to machine
    !                    precision. The recommended value is
    !                          cube root(4*EPSNEG)
    !      XUPPER - REAL(wp) - The value above which EXP3(X) = 0.89297...
    !                      to machine precision. The recommended value is
    !                           cube root(-ln(EPSNEG))
    !
    !      For values of EPS and EPSNEG for various machine/compiler
    !      combinations refer to the file MACHCON.TXT.
    !
    !      The machine-dependent constants are computed internally by
    !      using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED
    !      EXP, LOG
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: nterm1, nterm2
    REAL(wp)  :: t, x, xlow, xupper
    CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
    CHARACTER (LEN= 6)  :: fnname = 'EXP3  '

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                             two = 2.0_wp, three = 3.0_wp, four = 4.0_wp, &
                             sixten = 16.0_wp, onehun = 100.0_wp,  &
                             funinf = 0.89297951156924921122_wp
    REAL(wp), PARAMETER  :: aexp3(0:24) = [ &
        1.26919841422112601434_wp, -0.24884644638414098226_wp,  &
        0.8052622071723104125e-1_wp, -0.2577273325196832934e-1_wp,  &
        0.759987887307377429e-2_wp, -0.203069558194040510e-2_wp,  &
        0.49083458669932917e-3_wp, -0.10768223914202077e-3_wp, 0.2155172626428984e-4_wp, &
       -0.395670513738429e-5_wp, 0.66992409338956e-6_wp, -0.10513218080703e-6_wp, &
        0.1536258019825e-7_wp, -0.209909603636e-8_wp, 0.26921095381e-9_wp,  &
       -0.3251952422e-10_wp, 0.371148157e-11_wp, -0.40136518e-12_wp, 0.4123346e-13_wp, &
       -0.403375e-14_wp, 0.37658e-15_wp, -0.3362e-16_wp, 0.288e-17_wp, -0.24e-18_wp, 0.2e-19_wp]

    REAL(wp), PARAMETER  :: aexp3a(0:24) = [ &
        1.92704649550682737293_wp, -0.3492935652048138054e-1_wp, &
        0.145033837189830093e-2_wp, -0.8925336718327903e-4_wp, 0.705423921911838e-5_wp, &
       -0.66717274547611e-6_wp, 0.7242675899824e-7_wp, -0.878258256056e-8_wp, &
        0.116722344278e-8_wp, -0.16766312812e-9_wp, 0.2575501577e-10_wp, -0.419578881e-11_wp, &
        0.72010412e-12_wp, -0.12949055e-12_wp, 0.2428703e-13_wp, -0.473311e-14_wp,  &
        0.95531e-15_wp, -0.19914e-15_wp, 0.4277e-16_wp, -0.944e-17_wp, 0.214e-17_wp, -0.50e-18_wp, &
        0.12e-18_wp, -0.3e-19_wp, 0.1e-19_wp]

    ! Start calculation
    x = xvalue

    ! Error test
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    t = EPSILON(zero)
    xlow = (four*t) ** (one/three)
    xupper = (-LOG(t)) ** (one/three)
    t = t / onehun
    IF (x <= two) THEN
      DO  nterm1 = 24, 0, -1
        IF (ABS(aexp3(nterm1)) > t) EXIT
      END DO
    ELSE
      DO  nterm2 = 24, 0, -1
        IF (ABS(aexp3a(nterm2)) > t) EXIT
      END DO
    END IF

    ! Code for XVALUE < =  2
    IF (x <= two) THEN
      IF (x < xlow) THEN
        fn_val = x
      ELSE
        t = ((x*x*x/four)-half) - half
        fn_val = x * cheval(nterm1,aexp3,t)
      END IF
    ELSE
      ! Code for XVALUE > 2
      IF (x > xupper) THEN
        fn_val = funinf
      ELSE
        t = ((sixten/(x*x*x))-half) - half
        t = cheval(nterm2,aexp3a,t)
        t = t * EXP(-x*x*x) / (three*x*x)
        fn_val = funinf - t
      END IF
    END IF
    RETURN
  END FUNCTION exp3


  FUNCTION goodst(xvalue) RESULT(fn_val)
    !! CALGO 757 \(\int_{0}^{\infty} \frac{e^{-t^2}}{t + x} \, dt\)
    !
    !   DESCRIPTION:
    !      This function calculates the function defined as
    !
    !        GOODST(x) = {integral 0 to inf} ( exp(-u*u)/(u+x) ) du
    !
    !      The code uses Chebyshev expansions whose coefficients are
    !      given to 20 decimal places.
    !
    !   ERROR RETURNS:
    !      If XVALUE <= 0.0, an error message is printed, and the
    !      code returns the value 0.0.
    !
    !   MACHINE-DEPENDENT CONSTANTS:
    !      NTERM1 - The no. of terms to be used in the array AGOST.
    !                The recommended value is such that
    !                    AGOST(NTERM1) < EPS/100,
    !      NTERM2 - The no. of terms to be used in the array AGOSTA.
    !                The recommended value is such that
    !                    AGOSTA(NTERM2) < EPS/100,
    !      XLOW - The value below which f(x) = -(gamma/2) - ln(x)
    !             to machine precision. The recommended value is
    !                EPSNEG
    !      XHIGH - The value above which f(x) = sqrt(pi)/(2x) to
    !              machine precision. The recommended value is
    !                 2 / EPSNEG
    !
    !      For values of EPS and EPSNEG refer to the file MACHCON.TXT
    !
    !      The machine-dependent constants are computed internally by
    !      using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !       EXP , LOG
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: nterm1, nterm2
    REAL(wp)  :: fval, t, x, xhigh, xlow
    CHARACTER (LEN=15) :: errmsg = 'ARGUMENT <= 0.0'
    CHARACTER (LEN= 6) :: fnname = 'GOODST'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp,  &
                             two = 2.0_wp, six = 6.0_wp, onehun = 100.0_wp,  &
                             gamby2 = 0.28860783245076643030_wp,  &
                             rtpib2 = 0.88622692545275801365_wp

    REAL(wp), PARAMETER  :: agost(0:28) = [ &
        0.63106560560398446247_wp, 0.25051737793216708827_wp,  &
       -0.28466205979018940757_wp, 0.8761587523948623552e-1_wp,  &
        0.682602267221252724e-2_wp, -0.1081129544192254677e-1_wp,  &
        0.169101244117152176e-2_wp, 0.50272984622615186e-3_wp, -0.18576687204100084e-3_wp, &
       -0.428703674168474e-5_wp, 0.1009598903202905e-4_wp, -0.86529913517382e-6_wp,  &
       -0.34983874320734e-6_wp, 0.6483278683494e-7_wp, 0.757592498583e-8_wp,  &
       -0.277935424362e-8_wp, -0.4830235135e-10_wp, 0.8663221283e-10_wp, -0.394339687e-11_wp, &
       -0.209529625e-11_wp, 0.21501759e-12_wp, 0.3959015e-13_wp, -0.692279e-14_wp,  &
       -0.54829e-15_wp, 0.17108e-15_wp, 0.376e-17_wp, -0.349e-17_wp, 0.7e-19_wp, 0.6e-19_wp]

    REAL(wp), PARAMETER  :: agosta(0:23) = [ &
        1.81775467984718758767_wp, -0.9921146570744097467e-1_wp,  &
       -0.894058645254819243e-2_wp, -0.94955331277726785e-3_wp, -0.10971379966759665e-3_wp, &
       -0.1346694539578590e-4_wp, -0.172749274308265e-5_wp, -0.22931380199498e-6_wp,  &
       -0.3127844178918e-7_wp, -0.436197973671e-8_wp, -0.61958464743e-9_wp,  &
       -0.8937991276e-10_wp, -0.1306511094e-10_wp, -0.193166876e-11_wp, -0.28844270e-12_wp,  &
       -0.4344796e-13_wp, -0.659518e-14_wp, -0.100801e-14_wp, -0.15502e-15_wp, -0.2397e-16_wp,  &
       -0.373e-17_wp, -0.58e-18_wp, -0.9e-19_wp, -0.1e-19_wp]

    ! Start computation
    x = xvalue

    ! Error test
    IF (x <= zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    fval = EPSILON(zero)
    t = fval / onehun
    IF (x <= two) THEN
      DO  nterm1 = 28, 0, -1
        IF (ABS(agost(nterm1)) > t) EXIT
      END DO
      xlow = fval
    ELSE
      DO  nterm2 = 23, 0, -1
        IF (ABS(agosta(nterm2)) > t) EXIT
      END DO
      xhigh = two / fval
    END IF

    ! Computation for 0 < x <= 2
    IF (x <= two) THEN
      IF (x < xlow) THEN
        fn_val = -gamby2 - LOG(x)
      ELSE
        t = (x-half) - half
        fn_val = cheval(nterm1,agost,t) - EXP(-x*x) * LOG(x)
      END IF
    ELSE
      ! Computation for x > 2
      fval = rtpib2 / x
      IF (x > xhigh) THEN
        fn_val = fval
      ELSE
        t = (six-x) / (two+x)
        fn_val = fval * cheval(nterm2,agosta,t)
      END IF
    END IF
    RETURN
  END FUNCTION goodst


  FUNCTION lobach(xvalue) RESULT(fn_val)
    !! CALGO 757 Lobachevski's integral
    !! \(-\int_{0}^{\infty} \log \lvert \cos t \rvert \, dt\)
    !
    !   DESCRIPTION:
    !      This function calculates the Lobachewsky function L(x), defined as
    !
    !         LOBACH(x) = {integral 0 to x} ( -ln ( | cos t | ) dt
    !
    !      The code uses Chebyshev expansions whose coefficients are given
    !      to 20 decimal places.
    !
    !   ERROR RETURNS:
    !      If |x| too large, it is impossible to accurately reduce the
    !      argument to the range [0,pi]. An error message is printed
    !      and the program returns the value 0.0
    !
    !   MACHINE-DEPENDENT CONSTANTS:
    !      NTERM1 - INTEGER - The no. of terms to be used of the array ARLOB1.
    !                          The recommended value is such that
    !                          ABS(ARLOB1(NTERM1)) < EPS/100
    !      NTERM2 - INTEGER - The no. of terms to be used of the array ARLOB2.
    !                          The recommended value is such that
    !                          ABS(ARLOB2(NTERM2)) < EPS/100
    !      XLOW1 - DOUBLE PRECISION - The value below which L(x) = 0.0 to
    !                          machine-precision.
    !                     The recommended value is
    !                              cube-root ( 6*XMIN )
    !      XLOW2 - REAL(wp) - The value below which L(x) = x**3/6 to
    !                     machine-precision. The recommended value is
    !                              sqrt ( 10*EPS )
    !      XLOW3 - REAL(wp) - The value below which
    !                         L(pi/2) - L(pi/2-x) = x ( 1 - log(x) )
    !                     to machine-precision. The recommended value is
    !                               sqrt ( 18*EPS )
    !      XHIGH - REAL(wp) - The value of |x| above which it is impossible
    !                     to accurately reduce the argument. The
    !                     recommended value is   1 / EPS.
    !
    !      For values of EPS, and XMIN, refer to the file MACHCON.TXT
    !
    !      The machine-dependent constants are computed internally by
    !      using the D1MACH subroutine.
    !
    !   INTRINSIC FUNCTIONS USED:
    !      INT , LOG , SQRT
    !
    !   OTHER MISCFUN SUBROUTINES USED:
    !          CHEVAL , ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: indpi2, indsgn, npi, nterm1, nterm2
    REAL(wp)  :: fval, fval1, lbpb21, lobpi1, pi, piby2, piby21, piby4,  &
                  pi1, t, x, xcub, xhigh, xlow1, xlow2, xlow3, xr
    CHARACTER (LEN=26)  :: errmsg = 'ARGUMENT TOO LARGE IN SIZE'
    CHARACTER (LEN= 6)  :: fnname = 'LOBACH'

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                             two = 2.0_wp, six = 6.0_wp, ten = 10.0_wp,  &
                             onehun = 100.0_wp, lobpia = 1115.0_wp,  &
                             lobpib = 512.0_wp,  &
                             lobpi2 = -1.48284696397869499311e-4_wp, &
                             lbpb22 = -7.41423481989347496556e-5_wp, &
                             pi11 = 201.0_wp, pi12 = 64.0_wp,  &
                             pi2 = 9.67653589793238462643e-4_wp,  &
                             piby22 = 4.83826794896619231322e-4_wp,  &
                             tcon = 3.24227787655480868620_wp

    REAL(wp), PARAMETER  :: arlob1(0:15) = [ &
        0.34464884953481300507_wp, 0.584198357190277669e-2_wp,  &
        0.19175029694600330e-3_wp, 0.787251606456769e-5_wp,  &
        0.36507477415804e-6_wp, 0.1830287272680e-7_wp, 0.96890333005e-9_wp,  &
        0.5339055444e-10_wp, 0.303408025e-11_wp, 0.17667875e-12_wp,  &
        0.1049393e-13_wp, 0.63359e-15_wp, 0.3878e-16_wp, 0.240e-17_wp, 0.15e-18_wp, 0.1e-19_wp]
    REAL(wp), PARAMETER  :: arlob2(0:10) = [ &
        2.03459418036132851087_wp, 0.1735185882027407681e-1_wp,  &
        0.5516280426090521e-4_wp, 0.39781646276598e-6_wp, 0.369018028918e-8_wp,  &
        0.3880409214e-10_wp, 0.44069698e-12_wp, 0.527674e-14_wp, 0.6568e-16_wp, 0.84e-18_wp,  &
        0.1e-19_wp]

    ! Start computation
    x = ABS(xvalue)
    indsgn = 1
    IF (xvalue < zero) THEN
      indsgn = -1
    END IF

    ! Compute the machine-dependent constants.
    xr = EPSILON(zero)
    xhigh = one / xr

    ! Error test
    IF (x > xhigh) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! continue with constants
    t = xr / onehun
    DO  nterm1 = 15, 0, -1
      IF (ABS(arlob1(nterm1)) > t) EXIT
    END DO
    DO  nterm2 = 10, 0, -1
      IF (ABS(arlob2(nterm2)) > t) EXIT
    END DO
    xlow1 = (six*TINY(zero)) ** (two/six)
    xlow2 = SQRT(ten*xr)
    t = two * ten - two
    xlow3 = SQRT(t*xr)

    ! Reduce argument to [0,pi]
    pi1 = pi11 / pi12
    pi = pi1 + pi2
    piby2 = pi / two
    piby21 = pi1 / two
    piby4 = piby2 / two
    npi = INT(x/pi)
    xr = (x-npi*pi1) - npi * pi2

    ! Reduce argument to [0,pi/2]
    indpi2 = 0
    IF (xr > piby2) THEN
      indpi2 = 1
      xr = (pi1-xr) + pi2
    END IF

    ! Code for argument in [0,pi/4]
    IF (xr <= piby4) THEN
      IF (xr < xlow1) THEN
        fval = zero
      ELSE
        xcub = xr * xr * xr
        IF (xr < xlow2) THEN
          fval = xcub / six
        ELSE
          t = (tcon*xr*xr-half) - half
          fval = xcub * cheval(nterm1,arlob1,t)
        END IF
      END IF
    ELSE
      ! Code for argument in [pi/4,pi/2]
      xr = (piby21-xr) + piby22
      IF (xr == zero) THEN
        fval1 = zero
      ELSE
        IF (xr < xlow3) THEN
          fval1 = xr * (one-LOG(xr))
        ELSE
          t = (tcon*xr*xr-half) - half
          fval1 = xr * (cheval(nterm2,arlob2,t)-LOG(xr))
        END IF
      END IF
      lbpb21 = lobpia / (lobpib+lobpib)
      fval = (lbpb21-fval1) + lbpb22
    END IF
    lobpi1 = lobpia / lobpib

    ! Compute value for argument in [pi/2,pi]
    IF (indpi2 == 1) THEN
      fval = (lobpi1-fval) + lobpi2
    END IF
    fn_val = fval

    ! Scale up for arguments > pi
    IF (npi > 0) THEN
      fn_val = (fval+npi*lobpi2) + npi * lobpi1
    END IF
    IF (indsgn == -1) THEN
      fn_val = -fn_val
    END IF
    RETURN
  END FUNCTION lobach


  FUNCTION strom(xvalue) RESULT(fn_val)
    !! CALGO 757 Stromgren's integral
    !! \(\frac{15}{4 \pi^4} \int_{0}^{x} \frac{t^7 e^{2 t}}{(e^t - 1)^3} \, dt\)
    !
    !  DESCRIPTION:
    !    This program calculates Stromgren's integral, defined as
    !
    !      STROM(X) = integral 0 to X { t**7 exp(2t)/[exp(t)-1]**3 } dt
    !
    !    The code uses a Chebyshev series, the coefficients of which are
    !    given to an accuracy of 20 decimal places.
    !
    !  ERROR RETURNS:
    !    If XVALUE < 0.0, an error message is printed, and the program
    !    returns the value 0.0.
    !
    !  MACHINE-DEPENDENT CONSTANTS:
    !    NTERMS - INTEGER - The number of terms of the array ASTROM to be used.
    !                       The recommended value is such that
    !                             ASTROM(NTERMS) < EPS/100
    !    XLOW0 - REAL(wp) - The value below which STROM = 0.0 to machine
    !                    precision. The recommended value is
    !                          5th root of (130*XMIN)
    !    XLOW1 - REAL(wp) - The value below which STROM = 3*(X**5)/(4*(pi**4))
    !                   to machine precision. The recommended value is
    !                             2*EPSNEG
    !    EPSLN - REAL(wp) - The value of ln(EPS). Used to determine the no.
    !                   of exponential terms for large X.
    !    EPNGLN - REAL(wp) - The value of ln(EPSNEG). Used to prevent
    !                    overflow for large X.
    !    XHIGH - REAL(wp) - The value above which
    !                           STROM = 196.52 - 15*(x**7)*exp(-x)/(4pi**4)
    !                   to machine precision. The recommended value is
    !                             7 / EPS
    !
    !     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
    !
    !     The machine-dependent constants are computed internally by
    !     using the D1MACH subroutine.
    !
    ! INTRINSIC FUNCTIONS USED:
    !     EXP, INT, LOG
    !
    ! OTHER MISCFUN SUBROUTINES USED:
    !     CHEVAL, ERRPRN, D1MACH

    REAL(wp), INTENT(IN)  :: xvalue
    REAL(wp)              :: fn_val

    INTEGER    :: k1, k2, nterms, numexp
    REAL(wp)  :: epngln, epsln, rk, sumexp, sum2, t, x, xhigh, xk, xk1, xlow0, &
                  xlow1
    CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
    CHARACTER (LEN= 6)  :: fnname = 'STROM '

    REAL(wp), PARAMETER  :: zero = 0.0_wp, half = 0.5_wp, one = 1.0_wp,  &
                             two = 2.0_wp, four = 4.0_wp, seven = 7.0_wp,  &
                             onehun = 100.0_wp, one30 = 130.0_wp, one5ln = 0.4055_wp, &
                             f15bp4 = 0.38497433455066256959e-1_wp,  &
                             pi4b3 = 1.29878788045336582982D2,    &
                             valinf = 196.51956920868988261257_wp

    REAL(wp), PARAMETER  :: astrom(0:26) = [ &
        0.56556120872539155290_wp, 0.4555731969101785525e-1_wp,  &
       -0.4039535875936869170e-1_wp, -0.133390572021486815e-2_wp,   &
        0.185862506250538030e-2_wp, -0.4685555868053659e-4_wp, -0.6343475643422949e-4_wp,  &
        0.572548708143200e-5_wp, 0.159352812216822e-5_wp, -0.28884328431036e-6_wp,  &
       -0.2446633604801e-7_wp, 0.1007250382374e-7_wp, -0.12482986104e-9_wp,  &
       -0.26300625283e-9_wp, 0.2490407578e-10_wp, 0.485454902e-11_wp, -0.105378913e-11_wp,  &
       -0.3604417e-13_wp, 0.2992078e-13_wp, -0.163971e-14_wp, -0.61061e-15_wp, 0.9335e-16_wp,  &
        0.709e-17_wp, -0.291e-17_wp, 0.8e-19_wp, 0.6e-19_wp, -0.1e-19_wp]

    ! Start execution
    x = xvalue

    ! Error test
    IF (x < zero) THEN
      CALL errprn(fnname,errmsg)
      fn_val = zero
      RETURN
    END IF

    ! Compute the machine-dependent constants.
    xk = EPSILON(zero)
    t = xk / onehun
    IF (x <= four) THEN
      DO  nterms = 26, 0, -1
        IF (ABS(astrom(nterms)) > t) EXIT
      END DO
      xlow0 = (one30*TINY(zero)) ** (one/(seven-two))
      xlow1 = two * xk
    ELSE
      epsln = LOG(2*EPSILON(zero))
      epngln = LOG(xk)
      xhigh = seven / xk
    END IF

    ! Code for x < =  4.0
    IF (x <= four) THEN
      IF (x < xlow0) THEN
        fn_val = zero
      ELSE
        IF (x < xlow1) THEN
          fn_val = (x**5) / pi4b3
        ELSE
          t = ((x/two)-half) - half
          fn_val = (x**5) * cheval(nterms,astrom,t) * f15bp4
        END IF
      END IF
    ELSE
      ! Code for x > 4.0
      IF (x > xhigh) THEN
        sumexp = one
      ELSE
        numexp = INT(epsln/(one5ln-x)) + 1
        IF (numexp > 1) THEN
          t = EXP(-x)
        ELSE
          t = one
        END IF
        rk = zero
        DO  k1 = 1, numexp
          rk = rk + one
        END DO
        sumexp = zero
        DO  k1 = 1, numexp
          sum2 = one
          xk = one / (rk*x)
          xk1 = one
          DO  k2 = 1, 7
            sum2 = sum2 * xk1 * xk + one
            xk1 = xk1 + one
          END DO
          sum2 = sum2 * (rk+one) / two
          sumexp = sumexp * t + sum2
          rk = rk - one
        END DO
      END IF
      t = seven * LOG(x) - x + LOG(sumexp)
      IF (t < epngln) THEN
        fn_val = valinf
      ELSE
        fn_val = valinf - EXP(t) * f15bp4
      END IF
    END IF
    RETURN
  END FUNCTION strom


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
    !
    ! MACHINE-DEPENDENT CONSTANTS: NONE
    !
    ! INTRINSIC FUNCTIONS USED:
    !    ABS

    INTEGER, INTENT(IN)  :: n
    REAL(dp), INTENT(IN) :: a(0:n)
    REAL(dp), INTENT(IN) :: t
    REAL(dp)             :: fn_val
    INTEGER  :: i
    REAL(dp) :: d1, d2, tt, u0, u1, u2
    REAL(dp), PARAMETER :: zero = 0.0_wp, half = 0.5_wp, test = 0.6_wp, two = 2.0_wp

    u1 = zero

    !   If ABS ( T )  < 0.6 use the standard Clenshaw method
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
      ! If ABS ( T )  > =  0.6 use the Reinsch modification
      d1 = zero
      
      ! T > =  0.6 code
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
        ! T < =  -0.6 code
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


  SUBROUTINE errprn(fnname, errmsg)
    ! DESCRIPTION:
    !    This subroutine prints out an error message if
    !    an error has occurred in one of the MISCFUN functions.
    !
    ! INPUT PARAMETERS:
    !    FNNAME - CHARACTER - The name of the function with the error.
    !    ERRMSG - CHARACTER - The message to be printed out.
    !
    ! MACHINE-DEPENDENT PARAMETER:
    !    OUTSTR - INTEGER - The numerical value of the output stream to be used
    !                       for printing the error message.
    !                       The subroutine has the default value   OUTSTR = 6.

    CHARACTER(LEN = 6), INTENT(IN) :: fnname
    CHARACTER(LEN = *), INTENT(IN) :: errmsg

    INTEGER :: outstr = 6

    WRITE(outstr, 5000) fnname
    WRITE(outstr, 5100) errmsg
    RETURN

  5000 FORMAT(/t6, 'ERROR IN MISCFUN FUNCTION  ', a6)
  5100 FORMAT(/t6, a50)

  END SUBROUTINE errprn

end module calgo_757
