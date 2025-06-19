! Licensed under the ACM Software License Agreement
! Copyright © 1970–2012 Association for Computing Machinery (ACM)
! See ColSpecF LICENSE file for details.

MODULE calgo_683
!* # CALGO 683
! Algorithm 683. 
!
! Procedures:
!
! - `cexint`: Exponential integral \(\mathrm{E}_n(z)\)
!
! Untested procedures:  
!
! - `g8`: Gauss-Legendre quadrature with 8 points
! - `gaus8`: Adaptive Gauss-Legendre quadrature with 8 points
! - `psixn`: Digamma function \(\psi(n)\) for a positive integer \(n\)
!
! ## Author
! Donald E. Amos
!
! ## History
! - 1987-05-15 - Donald E. Amos
!     - Original code.
! - 1999-12-28 - Donald E. Amos
!     - F77 code distributed by ACM: <https://calgo.acm.org/683.zip>
! - YYYY-mm-dd - Alan Miller
!     - F90 code adaptation by Alan Miller:
!       <https://jblevins.org/mirror/amiller/toms683.f90>
! - 2025-06-05 - Rodrigo Castro (GitHub: rodpcastro)
!     - Removed procedures `cexqad` and `fqcex`, which were originally used for
!       testing.
! - 2025-06-06 - Rodrigo Castro (GitHub: rodpcastro)
!     - Created abstract interface for single-variable function, which is used by `g8`
!       and `gaus8`.
!     - Replaced `dp` (double precision) by `wp` (working precision)
!     - Fixed typo at `gaus8`:
!         - `anib = LOG10(DBLE(RADIX(0.0_wp))) * k / 0.30102000_wp`
!         - to `anib = LOG10(DBLE(RADIX(0.0_wp))) * k / 0.30103000_wp`
!
! ## References
! 1. Donald E. Amos. 1990. Algorithms 683: a portable FORTRAN subroutine for
!    exponential integrals of a complex argument. ACM Trans. Math. Softw. 16,
!*   2 (June 1990), 178–182. <https://doi.org/10.1145/78928.78934>

  use csf_kinds, only: wp

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cexint, g8, gaus8, psixn

  REAL(wp), SAVE :: fnm, gln, x, y
  INTEGER, SAVE  :: iflag

  ! Signature of single-variable function.
  abstract interface
    function funx(x) result(y)
      import :: wp
      real(wp), intent(in) :: x
      real(wp) :: y
    end function funx
  end interface

CONTAINS

  FUNCTION g8(fun, x, h) RESULT(fn_val)
    !! CALGO 683 Gauss-Legendre quadrature with 8 points.
    !
    ! Reference: https://en.wikipedia.org/wiki/Gaussian_quadrature#Change_of_interval
    !
    !* To evaluate \(\int_{a}^{b} f(x) dx\), the user must provide:
    !
    ! - `fun`: Single-variable function \(f(x)\)
    ! - `x = (a + b) / 2.0_wp`
    ! - `h = (b - a) / 2.0_wp`
    !*

    procedure(funx)      :: fun
    REAL(wp), INTENT(IN) :: x, h
    REAL(wp)             :: fn_val

    REAL(wp), PARAMETER :: x1 = 1.83434642495649805e-1_wp, &
                           x2 = 5.25532409916328986e-1_wp, &
                           x3 = 7.96666477413626740e-1_wp, &
                           x4 = 9.60289856497536232e-1_wp
    REAL(wp), PARAMETER :: w1 = 3.62683783378361983e-1_wp, &
                           w2 = 3.13706645877887287e-1_wp, &
                           w3 = 2.22381034453374471e-1_wp, &
                           w4 = 1.01228536290376259e-1_wp

    fn_val = h * ( &
        w1*(fun(x-x1*h) + fun(x+x1*h)) &
      + w2*(fun(x-x2*h) + fun(x+x2*h)) &
      + w3*(fun(x-x3*h) + fun(x+x3*h)) &
      + w4*(fun(x-x4*h) + fun(x+x4*h)) &
    )
    RETURN
  END FUNCTION g8

  SUBROUTINE gaus8(fun, a, b, ERR, ans, ierr)
    !! CALGO 683 Adaptive Gauss-Legendre quadrature with 8 points.
    !
    !! The description of all parameters and outputs can be found in the source code
    !! docstring.
    !
    ! WRITTEN BY R.E. JONES
    !
    ! ABSTRACT
    !    GAUS8 INTEGRATES REAL FUNCTIONS OF ONE VARIABLE OVER FINITE INTERVALS
    !    USING AN ADAPTIVE 8-POINT LEGENDRE-GAUSS ALGORITHM.
    !    GAUS8 IS INTENDED PRIMARILY FOR HIGH ACCURACY INTEGRATION OR
    !    INTEGRATION OF SMOOTH FUNCTIONS.
    !
    ! DESCRIPTION OF ARGUMENTS
    !
    !    INPUT--
    !    FUN - NAME OF EXTERNAL FUNCTION TO BE INTEGRATED.  THIS NAME MUST BE
    !          IN AN EXTERNAL STATEMENT IN THE CALLING PROGRAM.
    !          FUN MUST BE A FUNCTION OF ONE REAL ARGUMENT.  THE VALUE OF THE
    !          ARGUMENT TO FUN IS THE VARIABLE OF INTEGRATION WHICH RANGES
    !          FROM A TO B.
    !    A   - LOWER LIMIT OF INTEGRAL
    !    B   - UPPER LIMIT OF INTEGRAL (MAY BE LESS THAN A)
    !    ERR - IS A REQUESTED PSEUDORELATIVE ERROR TOLERANCE.  NORMALLY PICK A
    !          VALUE OF ABS(ERR) SO THAT STOL < ABS(ERR) <= 1.0E-3 WHERE STOL
    !          IS THE DOUBLE PRECISION UNIT ROUNDOFF = EPSILON(0.0_wp).
    !          ANS WILL NORMALLY HAVE NO MORE ERROR THAN ABS(ERR) TIMES THE
    !          INTEGRAL OF THE ABSOLUTE VALUE OF FUN(X).
    !          USUALLY, SMALLER VALUES FOR ERR YIELD MORE ACCURACY AND
    !          REQUIRE MORE FUNCTION EVALUATIONS.
    !   
    !          A NEGATIVE VALUE FOR ERR CAUSES AN ESTIMATE OF THE
    !          ABSOLUTE ERROR IN ANS TO BE RETURNED IN ERR.  NOTE THAT
    !          ERR MUST BE A VARIABLE (NOT A CONSTANT) IN THIS CASE.
    !          NOTE ALSO THAT THE USER MUST RESET THE VALUE OF ERR
    !          BEFORE MAKING ANY MORE CALLS THAT USE THE VARIABLE ERR.
    !
    !    OUTPUT--
    !    ERR - WILL BE AN ESTIMATE OF THE ABSOLUTE ERROR IN ANS IF THE
    !          INPUT VALUE OF ERR WAS NEGATIVE.  (ERR IS UNCHANGED IF
    !          THE INPUT VALUE OF ERR WAS NONNEGATIVE.)  THE ESTIMATED
    !          ERROR IS SOLELY FOR INFORMATION TO THE USER AND SHOULD
    !          NOT BE USED AS A CORRECTION TO THE COMPUTED INTEGRAL.
    !    ANS - COMPUTED VALUE OF INTEGRAL
    !    IERR- A STATUS CODE
    !        --NORMAL CODES
    !           1 ANS MOST LIKELY MEETS REQUESTED ERROR TOLERANCE, OR A=B.
    !          -1 A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMAL INTEGRATION.
    !             ANS IS SET TO ZERO.
    !        --ABNORMAL CODE
    !           2 ANS PROBABLY DOES NOT MEET REQUESTED ERROR TOLERANCE.

    procedure(funx)          :: fun
    REAL(wp), INTENT(IN)     :: a
    REAL(wp), INTENT(IN)     :: b
    REAL(wp), INTENT(IN OUT) :: ERR
    REAL(wp), INTENT(OUT)    :: ans
    INTEGER, INTENT(OUT)     :: ierr

    INTEGER  :: k, l, lmn, lmx, lr(30), mxl, nbits, nib, nlmx
    REAL(wp) :: aa(30), ae, anib, area, c, ce, ee, ef, eps, &
                est, gl, glr, gr(30), hh(30), tol, vl(30), vr
    INTEGER, PARAMETER  :: nlmn = 1, kmx = 5000, kml = 6
    INTEGER, SAVE       :: icall = 0
    real(wp), parameter :: sq2 = 1.41421356_wp

    ! INITIALIZE

    IF (icall /= 0) THEN
      WRITE(*, *) 'GAUS8- GAUS8 CALLED RECURSIVELY; NOT ALLOWED HERE'
      RETURN
    END IF

    icall = 1
    k = DIGITS(0.0_wp)
    anib = LOG10(DBLE(RADIX(0.0_wp))) * k / 0.30103000_wp
    nbits = INT(anib)
    nlmx = (nbits*5) / 8
    ans = 0.0_wp
    ierr = 1
    ce = 0.0_wp
    IF (a /= b) THEN
      lmx = nlmx
      lmn = nlmn
      IF (b /= 0.0_wp) THEN
        IF (SIGN(1.0_wp,b)*a > 0.0_wp) THEN
          c = ABS(1.0_wp-a/b)
          IF (c <= 0.1_wp) THEN
            IF (c <= 0.0_wp) GO TO 100
            anib = 0.5_wp - LOG(c) / 0.69314718_wp
            nib = anib
            lmx = MIN(nlmx, nbits-nib-7)
            IF (lmx < 1) GO TO 90
            lmn = MIN(lmn,lmx)
          END IF
        END IF
      END IF
      tol = MAX(ABS(ERR), 2.0_wp**(5-nbits)) / 2.0_wp
      IF (ERR == 0.0_wp) tol = SQRT(EPSILON(0.0_wp))
      eps = tol
      hh(1) = (b-a) / 4.0_wp
      aa(1) = a
      lr(1) = 1
      l = 1
      est = g8(fun, aa(l)+2.0_wp*hh(l), 2.0_wp*hh(l))
      k = 8
      area = ABS(est)
      ef = 0.5_wp
      mxl = 0
      
    ! COMPUTE REFINED ESTIMATES, ESTIMATE THE ERROR, ETC.
      
  10  gl = g8(fun, aa(l)+hh(l), hh(l))
      gr(l) = g8(fun, aa(l)+3.0_wp*hh(l), hh(l))
      k = k + 16
      area = area + (ABS(gl) + ABS(gr(l)) - ABS(est))
    ! IF (L < LMN) GO TO 11
      glr = gl + gr(l)
      ee = ABS(est-glr) * ef
      ae = MAX(eps*area, tol*ABS(glr))
      IF (ee-ae > 0.0) THEN
        GO TO 40
      ELSE
        GO TO 30
      END IF
  20  mxl = 1
  30  ce = ce + (est-glr)
      IF (lr(l) > 0) THEN
        GO TO 70
      ELSE
        GO TO 50
      END IF
      
    ! CONSIDER THE LEFT HALF OF THIS LEVEL
      
  40  IF (k > kmx) lmx = kml
      IF (l >= lmx) GO TO 20
      l = l + 1
      eps = eps * 0.5_wp
      ef = ef / sq2
      hh(l) = hh(l-1) * 0.5_wp
      lr(l) = -1
      aa(l) = aa(l-1)
      est = gl
      GO TO 10
      
    ! PROCEED TO RIGHT HALF AT THIS LEVEL
      
  50  vl(l) = glr
  60  est = gr(l-1)
      lr(l) = 1
      aa(l) = aa(l) + 4.0_wp * hh(l)
      GO TO 10
      
    ! RETURN ONE LEVEL
      
  70    vr = glr
  80    IF (l > 1) THEN
          l = l - 1
          eps = eps * 2.0_wp
          ef = ef * sq2
        IF (lr(l) <= 0) THEN
          vl(l) = vl(l+1) + vr
          GO TO 60
        END IF
        vr = vl(l+1) + vr
        GO TO 80
      END IF
      
    ! EXIT
      
      ans = vr
      IF (mxl == 0 .OR. ABS(ce) <= 2.0_wp*tol*area) GO TO 100
      ierr = 2
      WRITE(*, *) 'GAUS8- ANS IS PROBABLY INSUFFICIENTLY ACCURATE.'
      GO TO 100

  90  ierr = -1
      WRITE(*, *) 'GAUS8- THE FOLLOWING TEMPORARY DIAGNOSTIC WILL APPEAR ONLY ONCE.'
      WRITE(*, *) 'A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMAL INTEGRATION.'
      WRITE(*, *) 'ANS IS SET TO ZERO, AND IERR=-1.'
    END IF
  100 icall = 0
    IF (ERR < 0.0_wp) ERR = ce
    RETURN
  END SUBROUTINE gaus8

  SUBROUTINE cexint(z, n, kode, tol, m, cy, ierr)
    !! CALGO 683 Exponential integral \(\mathrm{E}_n(z)\).
    !
    !! \(n \geq 1,\thinspace
    !! \lbrace z \in \mathbb{C} \mid -\pi \lt \arg(z) \leq \pi \rbrace \)
    !
    !! The description of all parameters and outputs can be found in the source code
    !! docstring.
    !
    ! ON KODE=1, CEXINT COMPUTES AN M MEMBER SEQUENCE OF COMPLEX(wp)
    ! EXPONENTIAL INTEGRALS CY(J)=E(N+J-1,Z), J=1,...,M, FOR
    ! POSITIVE ORDERS N,...,N+M-1 AND COMPLEX(wp) Z IN THE CUT PLANE
    ! -PI < ARG(Z) <= PI (N=1 AND Z=CMPLX(0.0,0.0) CANNOT HOLD AT
    ! THE SAME TIME).  ON KODE=2, CEXINT COMPUTES SCALED FUNCTIONS
    !
    !                  CY(J)=E(N+J-1,Z)*CEXP(Z),      J=1,...,M,
    !
    ! WHICH REMOVES THE EXPONENTIAL BEHAVIOR IN BOTH THE LEFT AND
    ! RIGHT HALF PLANES.  DEFINITIONS AND NOTATION ARE FOUND IN THE
    ! NBS HANDBOOK OF MATHEMATICAL FUNCTIONS (REF. 1).
    !
    ! INPUT
    !   Z      - Z=CMPLX(X,Y), -PI < ARG(Z) <= PI
    !   N      - INTEGER ORDER OF INITIAL E FUNCTION, N=1,2,...
    !            (N=1 AND Z=CMPLX(0.0,0.0) IS AN ERROR)
    !   KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
    !            KODE= 1  RETURNS
    !                     CY(J)=E(N+J-1,Z),          J=1,...,M
    !                = 2  RETURNS
    !                     CY(J)=E(N+J-1,Z)*CEXP(Z),  J=1,...,M
    !   TOL    - PRECISION (ACCURACY) DESIRED FOR THE SEQUENCE,
    !            URND <= TOL <= 1.0E-3, WHERE URND IS LIMITED BY
    !            URND = MAX(UNIT ROUNDOFF,1.0E-18) AND UNIT
    !            ROUNDOFF = EPSILON(0.0_wp)
    !   M      - NUMBER OF E FUNCTIONS IN THE SEQUENCE, M >= 1
    !
    ! OUTPUT
    !   CY     - A COMPLEX(wp) VECTOR WHOSE FIRST M COMPONENTS CONTAIN
    !            VALUES FOR THE SEQUENCE
    !            CY(J)=E(N+J-1,Z)  OR
    !            CY(J)=E(N+J-1,Z)*CEXP(Z), J=1,...,M
    !            DEPENDING ON KODE.
    !   IERR   - ERROR FLAG
    !            IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
    !            IERR=1, INPUT ERROR   - NO COMPUTATION
    !            IERR=2, UNDERFLOW     - FIRST M COMPONENTS OF CY
    !                    SET TO ZERO, CY(J)=CMPLX(0.0,0.0), J=1,M,
    !                    REAL(Z) > 0.0 TOO LARGE ON KODE=1
    !            IERR=3, OVERFLOW      - NO COMPUTATION,
    !                    REAL(Z) < 0.0 TOO SMALL ON KODE=1
    !            IERR=4, CABS(Z) OR N+M-1 LARGE - COMPUTATION DONE
    !                    BUT LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION
    !                    PRODUCE LESS THAN HALF OF MACHINE ACCURACY
    !            IERR=5, CABS(Z) OR N+M-1 TOO LARGE - NO COMPUTATION BECAUSE
    !                    OF COMPLETE LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION
    !            IERR=6, ERROR         - NO COMPUTATION,
    !                    ALGORITHM TERMINATION CONDITION NOT MET.
    !                    SEE LONG DESCRIPTION ABOUT PARAMETER ICMAX.
    !            IERR=7, ERROR         - NO COMPUTATION,
    !                    DISCRIMINATION ERROR.  THIS CONDITION SHOULD NEVER OCCUR.
    !
    ! LONG DESCRIPTION
    !
    ! CEXINT USES A COMBINATION OF POWER SERIES AND BACKWARD RECURRENCE
    ! DESCRIBED IN REF. 2 FOR THE COMPLEX(wp) Z PLANE EXCEPT FOR A STRIP
    ! 2*YB WIDE ABOUT THE NEGATIVE REAL AXIS, WHERE ANALYTIC CONTINUATION IS
    ! CARRIED OUT BY LIMITED USE OF TAYLOR SERIES.
    ! THE SWITCH FROM BACKWARD RECURRENCE TO TAYLOR SERIES IS NECESSARY BECAUSE
    ! BACKWARD RECURRENCE IS SLOWLY CONVERGENT NEAR THE NEGATIVE REAL AXIS.
    ! THE BOUNDARIES Y=-YB AND Y=YB WERE DETERMINED SO THAT BACKWARD RECURRENCE
    ! WOULD CONVERGE EASILY WITH N AS LARGE AS 100 AND TOL AS SMALL AS 1.0e-18.
    ! SUBROUTINE CEXENZ DOES THE BACKWARD RECURRENCE AND SUBROUTINE CACEXI DOES
    ! THE ANALYTIC CONTINUATION.  TO START THE CONTINUATION, CACEXI CALLS CEXENZ
    ! WITH ZB=CMPLX(X,YB).
    ! IF CEXENZ RETURNS IERR=6, THEN YB IS INCREASED BY 0.5 UNTIL CEXENZ RETURNS
    ! IERR=0 OR 10 TRIES, WHICHEVER COMES FIRST.
    ! WHEN IERR=0, THEN THE ANALYTIC CONTINUATION PROCEEDS VERTICALLY DOWN FROM
    ! ZB=CMPLX(X,YB) TO Z=CMPLX(X,Y), 0 <= Y < YB.
    ! CONJUGATION IS USED FOR Y < 0.  YB INCREASES AS TOL DECREASES TO KEEP
    ! CONVERGENCE RATES UP AND RECURRENCE DOWN.
    !
    ! PARAMETER ICDIM=250 ALLOCATES STORAGE FOR THE COEFFICIENTS OF THE BACKWARD
    ! RECURRENCE ALGORITHM.  IF THE ALGORITHM TERMINATION CONDITION IS NOT MET
    ! IN ICDIM STEPS, THEN RECURRENCE PROCEEDS WITH NO ADDITIONAL STORAGE UNTIL
    ! THE TERMINATION CONDITION IS MET OR THE LIMIT ICMAX=2000 IS EXCEEDED.
    ! THE PURPOSE OF STORAGE IS TO MAKE THE ALGORITHM MORE EFFICIENT.
    ! THE TERMINATION CONDITION IS MET IN LESS THAN 250 STEPS OVER MOST OF
    ! THE COMPLEX(wp) PLANE EXCLUDING THE STRIP ABS(Y) < YB, X < 0.
    ! EXCEPTIONS TO THIS RULE ARE GENERATED NEAR STRIP BOUNDARIES WHEN N+M-1
    ! AND ABS(Z) ARE LARGE AND NEARLY EQUAL.  IN THESE CASES, THE CONVERGENCE
    ! IS VERY SLOW AND ADDITIONAL RECURRENCE (UP TO ICMAX) MUST BE USED.
    ! ON THE OTHERHAND, THESE REGIONS OF SLOW CONVERGENCE ARE KEPT SMALL BY
    ! ADJUSTING YB AS A FUNCTION OF TOL.  THESE REGIONS COULD BE ELIMINATED
    ! ENTIRELY BY MAKING YB SUFFICIENTLY LARGE, BUT THE EXPENSE AND INSTABILITY
    ! OF CONTINUATION BY TAYLOR SERIES NOT ONLY GOES UP, BUT THE COMPUTATIONAL
    ! EXPENSE BECOMES EXCESSIVELY LARGE IN OTHER PARTS OF THE LEFT HALF PLANE
    ! (Y < YB) WHERE THE BACKWARD RECURRENCE ALGORITHM WOULD CONVERGE RAPIDLY.
    !
    ! DERIVATIVES FOR SUCCESSIVE POWER SERIES ARE NOT COMPUTED BY EVALUATING
    ! DERIVATIVES OF A PREVIOUS POWER SERIES.  BECAUSE OF THE RELATION
    !
    !   (1)           DE(N,Z)/DZ = - E(N-1,Z),
    !
    ! SUCCESSIVE DERIVATIVES AT Z ARE GIVEN BY LOWER ORDER FUNCTIONS AND CAN BE
    ! COMPUTED IN A STABLE FASHION BY BACKWARD RECURRENCE USING (2) PROVIDED
    ! THAT THE BEGINNING ORDER NUB IS SMALLER THAN THE ARGUMENT.
    ! TO ACHIEVE THIS FOR ALL INTERMEDIATE VALUES ZZ BETWEEN ZB AND Z, WE TAKE
    ! NUB=MINO(N+M-1,INT(CABS(Z)+0.5)).
    ! TO START, E(NUB,ZB) IS EVALUATED BY THE BACKWARD RECURRENCE ALGORITHM OF
    ! REF. 3.  TO CONTINUE THE FUNCTION FROM ZB TO Z VIA INTERMEDIATE VALUES ZZ,
    ! DERIVATIVES OF E(NUB,ZB) ARE COMPUTED BY BACKWARD RECURRENCE ON (2).
    ! THIS ALLOWS A STEP (NO LARGER THAN 0.5) TO ZZ FOR E(NUB,ZZ) USING THE
    ! TAYLOR SERIES ABOUT ZB.  NOW, WE APPLY (2) AGAIN STARTING AT E(NUB,ZZ) FOR
    ! THE DERIVATIVES AT ZZ AND TAKE ANOTHER STEP, ETC.  NOTICE THAT THE
    ! STABILITY CONDITION FOR BACKWARD RECURRENCE, NUB <= ABS(Z) <= ABS(ZZ)
    ! <= ABS(ZB), IS SATISFIED FOR ALL INTERMEDIATE VALUES ZZ.  THE FINAL
    ! SEQUENCE FOR ORDERS N,...,N+M-1 IS GENERATED FROM (2) BY BACKWARD
    ! RECURRENCE, FORWARD RECURRENCE OR BOTH ONCE E(NUB,Z) HAS BEEN COMPUTED.
    !
    ! RECURRENCE WITH THE RELATION
    !
    !    (2)     N*E(N+1,Z) + Z*E(N,Z) = CEXP(-Z)
    !
    ! IN A DIRECTION AWAY FROM THE INTEGER CLOSEST TO ABS(Z) IS STABLE.
    ! FOR NEGATIVE ORDERS, THE RECURRENCE
    !
    !       E( 0,Z) = CEXP(-Z)/Z
    !       E(-N,Z) = ( CEXP(-Z)+N*E(-N+1,Z) )/Z   ,N=1,2,...
    !
    ! IS NUMERICALLY STABLE FOR ALL Z.
    !
    ! THE (CAPITAL) SINE AND COSINE INTEGRALS CAN BE COMPUTED FROM
    !
    !         SI(Z) =  (E(1,I*Z)-E(1,-I*Z))/(2*I) + PI/2
    !         CI(Z) = -(E(1,I*Z)+E(1,-I*Z))/2
    !
    ! IN -PI/2 < ARG(Z) <= PI/2, (I**2=-1), WHILE THE PRINCIPAL
    ! VALUED EXPONENTIAL INTEGRAL EI(X) CAN BE COMPUTED FROM
    !
    !     EI( X) = -(E(1,-X+I*0)+E(1,-X-I*0))/2 = -REAL(E(1,-X))
    !     EI(-X) = -REAL(E(1,X))
    !
    ! FOR X > 0.0 TO AN ACCURACY TOL.  IF Z = X > 0 THEN THE REAL SINE AND
    ! COSINE INTEGRALS ARE GIVEN BY
    !
    !         SI(X) = AIMAG(E(1,I*X)) + PI/2
    !         CI(X) = -REAL(E(1,I*X)) .
    !
    ! THE ANALYTIC CONTINUATION TO OTHER SHEETS CAN BE DONE BY THE RELATIONS
    !
    ! E(N,Z*CEXP(2*PI*M*I)) = E(N,Z) - 2*PI*M*I*(-Z)**(N-1)/(N-1)!
    !
    ! WHERE M=+1 OR M=-1 AND I**2=-1.
    !
    ! REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND
    !         I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF COMMERCE, 1955.
    !
    !       COMPUTATION OF EXPONENTIAL INTEGRALS OF COMPLEX(wp) ARGUMENT
    !         BY D. E. AMOS, ACM TRANS. MATH. SOFTWARE
    !
    !       COMPUTATION OF EXPONENTIAL INTEGRALS BY D. E. AMOS, ACM TRANS.
    !         MATH. SOFTWARE, VOL 6, NO. 3 SEPTEMBER 1980, PP. 365-377;
    !         ALGORITHM 556, EXPONENTIAL INTEGRALS, PP. 420-428.
    !
    !       REMARK ON ALGORITHM 556
    !         BY D. E. AMOS, ACM TRANS. MATH. SOFTWARE, VOL 9, NO. 4
    !         DECEMBER 1983, P. 525.
    !
    !       UNIFORM ASYMPTOTIC EXPANSIONS FOR EXPONENTIAL INTEGRALS E(N,X)
    !         AND BICKLEY FUNCTIONS KI(N,X) BY D. E. AMOS, ACM TRANS. MATH.
    !         SOFTWARE, VOL 9, NO. 4, DECEMBER. 1983, PP. 467-479;
    !         ALGORITHM 609, A PORTABLE FORTRAN SUBROUTINE FOR BICKLEY
    !         FUNCTIONS KI(N,X), PP. 480-493.
    !
    ! ROUTINES CALLED CACEXI, CEXENZ

    COMPLEX(wp), INTENT(IN)  :: z
    INTEGER, INTENT(IN)      :: n
    INTEGER, INTENT(IN)      :: kode
    REAL(wp), INTENT(IN)     :: tol
    INTEGER, INTENT(IN)      :: m
    COMPLEX(wp), INTENT(OUT) :: cy(m)
    INTEGER, INTENT(OUT)     :: ierr

    INTEGER     :: i, k, k1, k2
    REAL(wp)    :: aa, alim, az, bb, ca(250), d, elim, &
                   fn, rbry, r1m5, urnd, x, y, yb, htol
    COMPLEX(wp) :: cb(250)
    !-----------------------------------------------------------------------
    ! DIMENSION CA(ICDIM),CB(ICDIM)
    !-----------------------------------------------------------------------
    INTEGER, PARAMETER :: icdim = 250

    ! FIRST EXECUTABLE STATEMENT  CEXINT
    ierr = 0
    x = REAL(z, KIND=wp)
    y = AIMAG(z)
    IF (x == 0.0_wp .AND. y == 0.0_wp .AND. n == 1) ierr = 1
    IF (n < 1) ierr = 1
    IF (kode < 1 .OR. kode > 2) ierr = 1
    IF (m < 1) ierr = 1
    urnd = MAX(EPSILON(0.0_wp), 1.0E-18_wp)
    IF (tol < urnd .OR. tol > 1.0E-3_wp) ierr = 1
    IF (ierr /= 0) RETURN
    IF (x /= 0.0_wp .OR. y /= 0.0_wp .OR. n <= 1) THEN
    !-----------------------------------------------------------------------
    ! SET PARAMETERS RELATED TO MACHINE CONSTANTS.
    ! URND IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
    ! ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
    ! EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/URND    AND
    ! EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*URND       ARE INTERVALS NEAR
    ! UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
    !-----------------------------------------------------------------------
      k1 = MINEXPONENT(0.0_wp)
      k2 = MAXEXPONENT(0.0_wp)
      r1m5 = LOG10(DBLE(RADIX(0.0_wp)))
      k = MIN(ABS(k1),ABS(k2))
      elim = 2.303_wp * (k*r1m5 - 3.0_wp)
      k1 = DIGITS(0.0_wp) - 1
      aa = r1m5 * k1
      aa = aa * 2.303
      alim = elim + MAX(-aa, -41.45_wp)
      rbry = 2.0
      IF (urnd > 1.0E-8_wp) rbry = 1.0_wp
    !-----------------------------------------------------------------------
    ! TEST VARIABLES FOR RANGE. ABS(Z) CANNOT BE LARGER THAN THE ARGUMENT
    ! OF THE INT( ) FUNCTION.
    !-----------------------------------------------------------------------
      az = ABS(z)
      fn = n+m-1
      aa = 0.5_wp / urnd
      bb = HUGE(0.0_wp) * 0.5_wp
      aa = MIN(aa,bb)
      IF (az > aa) GO TO 20
      IF (fn > aa) GO TO 20
      aa = SQRT(aa)
      IF (az > aa) ierr = 4
      IF (fn > aa) ierr = 4
      IF (x >= 0.0_wp) THEN
    !-----------------------------------------------------------------------
    ! BACKWARD RECURRENCE FOR THE RIGHT HALF PLANE, X >= 0.0E0
    !-----------------------------------------------------------------------
        CALL cexenz(z, n, kode, m, cy, ierr, rbry, tol, elim, alim, icdim, ca, cb)
        RETURN
      END IF
      IF (az <= rbry) THEN
    !-----------------------------------------------------------------------
    ! POWER SERIES FOR ABS(Z) <= RBRY AND X < 0.0E0
    !-----------------------------------------------------------------------
        CALL cexenz(z, n, kode, m, cy, ierr, rbry, tol, elim, alim, icdim, ca, cb)
        RETURN
      END IF
      d = -0.4342945_wp * LOG(tol)
      yb = 10.5_wp - 0.538460_wp * (18.0_wp-d)
      IF (ABS(y) >= yb) THEN
    !-----------------------------------------------------------------------
    ! BACKWARD RECURRENCE EXTERIOR TO THE STRIP ABS(Y) < YB, X < 0.0
    !-----------------------------------------------------------------------
        htol = 0.125_wp * tol
        CALL cexenz(z, n, kode, m, cy, ierr, rbry, htol, elim, alim, icdim, ca, cb)
        RETURN
      END IF
    !-----------------------------------------------------------------------
    ! TAYLOR SERIES IN CACEXI FOR ANALYTIC CONTINUATION
    !-----------------------------------------------------------------------
      CALL cacexi(z, n, kode, m, cy, ierr, yb, rbry, tol, elim, alim, icdim, ca, cb )
      RETURN
    END IF
    DO i = 1, m
      cy(i) = 1.0_wp / (n+i-2)
    END DO
    RETURN

  20 ierr = 5
    RETURN
  END SUBROUTINE cexint

  SUBROUTINE cexenz(z, n, kode, m, cy, ierr, rbry, tol, elim, alim, icdim, ca, cb)
    ! REFER TO CEXINT
    !
    ! CEXENZ COMPUTES THE EXPONENTIAL INTEGRAL BY MEANS OF POWER SERIES AND A
    ! BACKWARD RECURRENCE ALGORITHM FOR THE CONFLUENT HYPERGEOMETRIC
    ! REPRESENTATION.
    !
    ! ROUTINES CALLED PSIXN

    COMPLEX(wp), INTENT(IN)  :: z
    INTEGER, INTENT(IN)      :: n
    INTEGER, INTENT(IN)      :: kode
    INTEGER, INTENT(IN)      :: m
    COMPLEX(wp), INTENT(OUT) :: cy(m)
    INTEGER, INTENT(OUT)     :: ierr
    REAL(wp), INTENT(IN)     :: rbry
    REAL(wp), INTENT(IN)     :: tol
    REAL(wp), INTENT(IN)     :: elim
    REAL(wp), INTENT(IN)     :: alim
    INTEGER, INTENT(IN)      :: icdim
    REAL(wp), INTENT(OUT)    :: ca(icdim)
    COMPLEX(wp), INTENT(OUT) :: cb(icdim)

    INTEGER     :: i, ic, icase, ict, ik, ind, iz, jset, &
                   k, kflag, kn, ks, ml, mu, nd, nm
    REAL(wp)    :: aa, aam, aams, aem, ah, ak, ap1, at, az, bk, bt, &
                   dk, ERR, fc, fnm, rtola, tola, x, xtol, y, ck
    COMPLEX(wp) :: cp1, cp2, cpt, cat, cbt, cy1, cy2, cyy(2), cnorm, &
                   cs, cak, emz, caa, tz, fz, ct, scle, rscle

    REAL(wp), PARAMETER    :: euler = -5.77215664901532861e-1_wp
    COMPLEX(wp), PARAMETER :: czero = (0.0_wp, 0.0_wp), cone = (1.0_wp, 0.0_wp)
    INTEGER, PARAMETER     :: icmax = 2000

    ierr = 0
    scle = cone
    rscle = cone
    x = REAL(z, KIND=wp)
    y = AIMAG(z)
    az = ABS(z)
    IF (az <= rbry) THEN
    !-----------------------------------------------------------------------
    ! SERIES FOR E(N,Z) FOR ABS(Z) <= RBRY
    !-----------------------------------------------------------------------
      iz = az + 0.5_wp
    !-----------------------------------------------------------------------
    ! ICASE=1 MEANS INTEGER CLOSEST TO ABS(Z) IS 2 AND N=1
    ! ICASE=2 MEANS INTEGER CLOSEST TO ABS(Z) IS 0,1, OR 2 AND N >= 2
    !-----------------------------------------------------------------------
      icase = 2
      IF (iz > n) icase = 1
      nm = n - icase + 1
      nd = nm + 1
      ind = 3 - icase
      mu = m - ind
      ml = 1
      ks = nd
      fnm = nm
      cs = czero
      xtol = 0.3333_wp * tol
      aam = 1.0_wp
      IF (nd /= 1) THEN
        aam = 1.0_wp / fnm
        cs = aam
      END IF
      caa = cone
      aa = 1.0_wp
      ak = 1.0_wp
    !-----------------------------------------------------------------------
    ! LIMIT INDEX I TO IK TO PREVENT UNDERFLOW ON SMALL VALUES OF Z
    !-----------------------------------------------------------------------
      ik = 35
      IF (az < xtol*aam) ik = 1
      DO i = 1, ik
        at = 1.0_wp / ak
        caa = -caa * z * at
        aa = aa * az * at
        IF (i /= nm) THEN
          cs = cs - caa / (ak-fnm)
        ELSE
          cs = cs + caa * (-LOG(z) + psixn(nd))
        END IF
        IF (aa <= xtol*ABS(cs)) EXIT
        ak = ak + 1.0_wp
      END DO

      IF (nd == 1) cs = cs + (-LOG(z) + euler)
      IF (kode == 2) cs = cs * EXP(z)
      cy(1) = cs
      ct = cs
      emz = cone
      IF (m /= 1) THEN
        cy(ind) = cs
        ak = ks
        IF (kode == 1) emz = EXP(-z)
        SELECT CASE ( icase )
          CASE (    1)
            GO TO 140
          CASE (    2)
            GO TO 160
        END SELECT
      END IF
      IF (icase == 2) RETURN
      IF (kode == 1) emz = EXP(-z)
      cy(1) = (emz-cs) / z
      RETURN
    END IF
    !-----------------------------------------------------------------------
    ! BACKWARD RECURSIVE MILLER ALGORITHM FOR
    !   E(N,Z)=EXP(-Z)*(Z**(N-1))*U(N,N,Z)
    ! WITH RECURSION AWAY FROM N=INTEGER CLOSEST TO ABS(Z)
    ! U(A,B,Z) IS THE SECOND CONFLUENT HYPERGEOMETRIC FUNCTION
    !-----------------------------------------------------------------------
    emz = cone
    IF (kode /= 2) THEN
    !-----------------------------------------------------------------------
    ! SCALE NEAR EXPONENT EXTREMES ON KODE=1
    !-----------------------------------------------------------------------
      IF (x >= 0.0_wp) THEN
        at = x + (n+m-1)
        ct = CMPLX(at, y, KIND=wp)
        aa = ABS(ct)
        at = x + LOG(aa)
        IF (at > elim) GO TO 30
        kflag = 1
      ELSE
        at = x
        IF (at < (-elim)) GO TO 180
        kflag = 2
      END IF
      IF (ABS(at) >= alim) THEN
        tola = EXP(alim-elim)
        rtola = 1.0_wp / tola
        IF (kflag /= 2) THEN
          scle  = rtola
          rscle = tola
        ELSE
          scle  = tola
          rscle = rtola
        END IF
        emz = scle
      END IF
      emz = emz * EXP(-z)
    END IF
    iz = az + 0.5_wp
    kn = n + m - 1
    IF (kn <= iz) GO TO 50
    IF (n < iz .AND. iz < kn) GO TO 80
    IF (n >= iz) GO TO 70
    ierr = 7
    RETURN

  30 ierr = 2
    cy(1:m) = czero
    RETURN

  50 icase = 1
    ks = kn
    ml = m - 1
    mu = -1
    ind = m
    IF (kn > 1) GO TO 90

  60 ks = 2
    icase = 3
    GO TO 90

  70 icase = 2
    ind = 1
    ks = n
    mu = m - 1
    IF (n > 1) GO TO 90
    IF (kn == 1) GO TO 60
    iz = 2

  80 icase = 1
    ks = iz
    ml = iz - n
    ind = ml + 1
    mu = kn - iz

  90 ik = ks / 2
    ah = ik
    jset = 1 + ks - 2 * ik
    !-----------------------------------------------------------------------
    ! START COMPUTATION FOR
    !   CYY(1) = C*U( A , A ,Z)    JSET=1
    !   CYY(1) = C*U(A+1,A+1,Z)    JSET=2
    ! FOR AN EVEN INTEGER A.
    !-----------------------------------------------------------------------
    ic = 0
    aa = ah + ah
    caa = aa
    aam = aa - 1.0_wp
    aams = aam * aam
    tz = z + z
    fz = tz + tz
    ak = ah
    xtol = tol
    ct = aams + fz * ah
    cak = z + caa
    aem = (ak+1.0_wp) / xtol
    aem = aem / ABS(cak)
    bk = aa
    ck = ah * ah
    !-----------------------------------------------------------------------
    ! FORWARD RECURSION FOR P(IC),P(IC+1) AND INDEX IC FOR BACKWARD
    ! RECURSION
    !-----------------------------------------------------------------------
    cp1 = czero
    cp2 = cone

  100 ic = ic + 1
    IF (ic <= icdim) THEN
      ak = ak + 1.0_wp
      ck = ck + 1.0_wp
      at = bk / (bk+ak+ck)
      bk = bk + ak + ak
      cat = at
      ca(ic) = at
      bt = 1.0_wp / (ak + 1.0_wp)
      cbt = (ak + ak + z) * bt
      cb(ic) = cbt
      cpt = cp2
      cp2 = cbt * cp2 - cat * cp1
      cp1 = cpt
      ct = ct + fz
      aem = aem * at
      bt = ABS(ct)
      ERR = aem / SQRT(bt)
      ap1 = ABS(cp1)
      IF (ERR*(ak+1.0_wp)/ap1 > ap1) GO TO 100
    ELSE
    !-----------------------------------------------------------------------
    ! CONTINUE FORWARD RECURRENCE UNINDEXED WHEN IC EXCEEDS ICDIM
    !-----------------------------------------------------------------------
      ic = ic - 1

  110 ic = ic + 1
      IF (ic > icmax) GO TO 190
      ak = ak + 1.0_wp
      ck = ck + 1.0_wp
      at = bk / (bk+ak+ck)
      bk = bk + ak + ak
      cat = at
      bt = 1.0_wp / (ak+1.0_wp)
      cbt = (ak+ak + z) * bt
      cpt = cp2
      cp2 = cbt * cp2 - cat * cp1
      cp1 = cpt
      ct = ct + fz
      aem = aem * at
      bt = ABS(ct)
      ERR = aem / SQRT(bt)
      ap1 = ABS(cp1)
      IF (ERR*(ak+1.0_wp)/ap1 > ap1) GO TO 110
    END IF
    fc = ic
    at = ((fc+1.0_wp)/(ak+1.0_wp)) * ((ak+ah)/(ak+1.0_wp))
    cat = at * SQRT(ct/(ct+fz))
    cy2 = cat * (cp1/cp2)
    cy1 = cone
    !-----------------------------------------------------------------------
    ! BACKWARD RECURRENCE FOR
    !   CY1=C*U( A ,A,Z)
    !   CY2=C*(A/(1+A/2))*U(A+1,A,Z)
    !-----------------------------------------------------------------------
    IF (ic > icdim) THEN
    !-----------------------------------------------------------------------
    ! BACKWARD RECUR UNINDEXED WHEN IC EXCEEDS ICDIM
    !-----------------------------------------------------------------------
      bt = aa + x

  120 ak = ah + fc
      bk = ak + 1.0_wp
      ck = aam + fc
      dk = bt + fc + fc
      at = (ak/fc) * (bk/ck)
      cbt = CMPLX(dk/bk, y/bk, KIND=wp)
      cpt = cy1
      cy1 = (cbt*cy1 - cy2) * at
      cy2 = cpt
      fc = fc - 1.0_wp
      ic = ic - 1
      IF (ic > icdim) GO TO 120
    END IF
    ict = ic
    DO k = 1, ict
      at = 1.0_wp / ca(ic)
      cpt = cy1
      cy1 = (cb(ic)*cy1 - cy2) * at
      cy2 = cpt
      ic = ic - 1
    END DO
    !-----------------------------------------------------------------------
    ! THE CONTIGUOUS RELATION
    !   Z*U(B,C+1,Z)=(C-B)*U(B,C,Z)+U(B-1,C,Z)
    ! WITH  B=A+1 , C=A IS USED FOR
    !   CYY(2) = C * U(A+1,A+1,Z)
    ! Z IS INCORPORATED INTO THE NORMALIZING RELATION FOR CNORM.
    !-----------------------------------------------------------------------
    cpt = cy2 / cy1
    cnorm = cone - cpt * (ah + 1.0_wp) / aa
    cyy(1) = cone / (cnorm*caa+z)
    cyy(2) = cnorm * cyy(1)
    IF (icase /= 3) THEN
      ct = emz * cyy(jset)
      cy(ind) = ct * rscle
      cs = ct
      IF (m == 1) RETURN
      ak = ks
      SELECT CASE ( icase )
        CASE (    1)
          GO TO 140
        CASE (    2)
          GO TO 160
      END SELECT
    END IF
    !-----------------------------------------------------------------------
    ! RECURSION SECTION  N*E(N+1,Z) + Z*E(N,Z)=EMZ
    !-----------------------------------------------------------------------
    ct = emz * (cone - cyy(1)) / z
    cy(1) = ct * rscle
    RETURN

  140 caa = ak
    tz = cone / z
    k = ind - 1
    DO i = 1, ml
      caa = caa - cone
      ct = (emz-caa*ct) * tz
      cy(k) = ct * rscle
      k = k - 1
    END DO
    IF (mu <= 0) RETURN
    ak = ks

  160 k = ind
    DO i = 1, mu
      cs = (emz - z*cs) / ak
      cy(k+1) = cs * rscle
      ak = ak + 1.0_wp
      k = k + 1
    END DO
    RETURN

  180 ierr = 3
    RETURN

  190 ierr = 6
    RETURN
  END SUBROUTINE cexenz

  SUBROUTINE cacexi(z, nu, kode, n, y, ierr, yb, rbry, tol, elim, alim, icdim, ca, cb)
    ! REFER TO  CEXINT
    !
    ! CACEXI COMPUTES THE ANALYTIC CONTINUATION OF THE EXPONENTIAL INTEGRAL
    ! FOR X < 0 AND -YB < Y < YB BY TAYLOR SERIES IN INCREMENTS OF HALF A UNIT.
    ! THE CONTINUATION PROCEEDS VERTICALLY DOWN FROM ZB=CMPLX(X,YB) TO
    ! Z=CMPLX(X,Y) FOR 0.0 <= Y < YB.  CONJUGATION IS USED FOR Y < 0.0E0.
    !
    ! ROUTINES CALLED  CEXENZ

    COMPLEX(wp), INTENT(IN)     :: z
    INTEGER, INTENT(IN)         :: nu
    INTEGER, INTENT(IN)         :: kode
    INTEGER, INTENT(IN)         :: n
    COMPLEX(wp), INTENT(OUT)    :: y(n)
    INTEGER, INTENT(OUT)        :: ierr
    REAL(wp), INTENT(IN OUT)    :: yb
    REAL(wp), INTENT(IN OUT)    :: rbry
    REAL(wp), INTENT(IN)        :: tol
    REAL(wp), INTENT(IN)        :: elim
    REAL(wp), INTENT(IN)        :: alim
    INTEGER, INTENT(IN)         :: icdim
    REAL(wp), INTENT(IN OUT)    :: ca(icdim)
    COMPLEX(wp), INTENT(IN OUT) :: cb(icdim)

    INTEGER     :: i, iaz, il, is, iy, k, kb, kl, kmax, ks, kyb, nb, nflg, nub
    REAL(wp)    :: az, del, fk, rtola, tola, yt, zi, zid, zr, &
                   xtol, rzi, atrm, fj, rw, rq(64), asum, htol
    COMPLEX(wp) :: yy(1), cex, sum, trm, zz, zp(64), zt, scle, rscle, trms, zw, cezt
    COMPLEX(wp), PARAMETER :: cone = (1.0_wp, 0.0_wp)

    scle = cone
    rscle = cone
    zr = REAL(z, KIND=wp)
    zi = AIMAG(z)
    zid = zi
    kyb = 0
    IF (zi < 0.0_wp) zid = -zid
    az = ABS(z)
    iaz = az + 0.5_wp
    !-----------------------------------------------------------------------
    !     SET ORDER NUB=MIN(N+M-1,INT(ABS(Z)+0.5)), GENERATE REMAINDER
    !     OF THE SEQUENCE AFTER THE COMPUTATION OF E(NUB,Z)
    !-----------------------------------------------------------------------
    IF (nu < iaz) THEN
      IF (nu+n-1 <= iaz) GO TO 10
      nub = iaz
      nb = 0
      nflg = 3
      kb = nub - nu + 1
      ks = kb + 1
      kl = n
      is = 1
      il = kb - 1
      GO TO 30
    END IF
    nub = MAX(iaz, 1)
    nb = nu - nub
    nflg = 1
    kb = 1
    ks = 2
    kl = n
    GO TO 30

  10 nub = nu + n - 1
    nb = 0
    nflg = 2
    is = 2
    il = n
    kb = n
    GO TO 30
    !-----------------------------------------------------------------------
    !     SET PARAMETERS FOR ANALYTIC CONTINUATION FROM Y=YB INTO THE REGION
    !     0 <= ZID < YB.
    !-----------------------------------------------------------------------
  20 yb = yb + 0.5_wp
    kyb = kyb + 1
    IF (kyb > 10) RETURN
    ierr = 0

  30 del = yb - zid
    !-----------------------------------------------------------------------
    !     MAKE DEL LARGE ENOUGH TO AVOID UNDERFLOW IN GENERATION OF POWERS
    !-----------------------------------------------------------------------
    IF (ABS(del) <= 1.0E-4_wp) THEN
      yb = yb + 1.0E-4_wp
      del = yb - zid
    END IF
    htol = 0.125_wp * tol
    zz = CMPLX(zr, yb, KIND=wp)
    CALL cexenz(zz, nub, 2, 1, yy, ierr, rbry, htol, elim, alim, icdim, ca, cb)
    IF (ierr == 6) GO TO 20
    !-----------------------------------------------------------------------
    !     ANALYTIC CONTINUATION VIA TAYLOR SERIES FOR ORDER NUB
    !-----------------------------------------------------------------------
    iy = del + del
    yt = del / (iy+1)
    sum = yy(1)
    htol = 0.25_wp * tol
    trm = sum
    cezt = CMPLX(COS(yt), -SIN(yt), KIND=wp)
    zw = cone
    zp(1) = cone
    fk = 1.0_wp
    fj = nub - 1
    zt = cone / zz
    !-----------------------------------------------------------------------
    !     TERMS OF THE SERIES TRM=E(NUB-K,ZZ)*(YT**K)/K!,  K=0,1,... ARE
    !     GENERATED BY BACKWARD RECURRENCE.  E IS SCALED BY CEXP(ZZ).
    !-----------------------------------------------------------------------
    DO k = 2, 64
      rw = yt / fk
      zw = zw * CMPLX(0.0_wp, rw, KIND=wp)
      rw = fj * rw
      trm = (zw - CMPLX(0.0_wp, rw, KIND=wp)*trm) * zt
      sum = sum + trm
      zp(k) = zw
      rq(k) = rw
      asum = ABS(sum)
      atrm = ABS(trm)
      IF (atrm < htol*asum) GO TO 50
      fk = fk + 1.0_wp
      fj = fj - 1.0_wp
    END DO
    k = 64

  50 kmax = k
    sum = sum * cezt
    IF (iy /= 0) THEN
      DO i = 1, iy
        rzi = (iy-i+1) * yt + zid
        zz = CMPLX(zr, rzi, KIND=wp)
        zt = cone / zz
        trm = sum
        DO k = 2, kmax
          trm = (zp(k) - CMPLX(0.0_wp, rq(k), KIND=wp)*trm) * zt
          sum = sum + trm
        END DO
        atrm = ABS(trm)
        asum = ABS(sum)
        xtol = htol * asum
        IF (atrm >= xtol) THEN
          IF (kmax < 64) THEN
            kmax = kmax + 1
            DO k = kmax, 64
              rw = yt / fk
              zw = zw * CMPLX(0.0_wp, rw, KIND=wp)
              rw = fj * rw
              trm = (zw - CMPLX(0.0_wp, rw, KIND=wp)*trm) * zt
              sum = sum + trm
              zp(k) = zw
              rq(k) = rw
              atrm = ABS(trm)
              IF (atrm < xtol) GO TO 80
              fk = fk + 1.0_wp
              fj = fj - 1.0_wp
            END DO
            k = 64
  80        kmax = k
          END IF
        END IF
        sum = sum * cezt
      END DO
    END IF
    !-----------------------------------------------------------------------
    !     FORM SEQUENCE UP OR DOWN FROM ORDER NUB
    !-----------------------------------------------------------------------
    IF (zi < 0.0_wp) sum = CONJG(sum)
    cex = cone
    !-----------------------------------------------------------------------
    !     SCALE NEAR OVERFLOW LIMIT ON KODE=1
    !-----------------------------------------------------------------------
    IF (kode /= 2) THEN
      IF (ABS(zr) >= alim) THEN
        IF (ABS(zr) > elim) GO TO 130
        tola = EXP(alim-elim)
        rtola = 1.0_wp / tola
        scle = tola
        rscle = rtola
      END IF
      cex = scle * EXP(-z)
    END IF
    trm = sum * cex
    y(kb) = trm * rscle
    trms = trm
    IF (n == 1 .AND. nflg /= 1) RETURN
    IF (nflg /= 2) THEN
      fk = nub
      IF (nflg == 1 .AND. nb /= 0) THEN
        DO k = 1, nb
          trm = (cex - z*trm) / fk
          fk = fk + 1.0_wp
        END DO
      END IF
      y(kb) = trm * rscle
      trms = trm
      IF (n == 1) RETURN
      DO k = ks, kl
        trm = (cex - z*trm) / fk
        y(k) = trm * rscle
        fk = fk + 1.0_wp
      END DO
      IF (nflg == 1) RETURN
    END IF
    k = kb - 1
    fk = nub - 1
    zt = cone / z
    DO i = is, il
      trms = (cex - fk*trms) * zt
      y(k) = trms * rscle
      k = k - 1
      fk = fk - 1.0_wp
    END DO
    RETURN

  130 ierr = 3
    RETURN
  END SUBROUTINE cacexi

  FUNCTION psixn(n) RESULT(fn_val)
    !! CALGO 683 Digamma function \(\psi_0(n)\).
    !
    !! \(\lbrace n \in \mathbb{Z} \mid n \geq 1 \rbrace\)
    !
    ! REFER TO CEXENZ
    !
    ! THIS SUBROUTINE RETURNS VALUES OF PSI(X)=DERIVATIVE OF LOG GAMMA(X),
    ! X > 0.0 AT INTEGER ARGUMENTS. A TABLE LOOK-UP IS PERFORMED FOR N <= 100,
    ! AND THE ASYMPTOTIC EXPANSION IS EVALUATED FOR N > 100.

    INTEGER, INTENT(IN) :: n
    REAL(wp)            :: fn_val

    INTEGER  :: k
    REAL(wp) :: ax, fn, rfn2, trm, s, wdtol
    !-----------------------------------------------------------------------
    ! PSIXN(N), N = 1,100
    !-----------------------------------------------------------------------
    REAL(wp), PARAMETER  :: c(100) = (/ -5.77215664901532861e-1_wp,  &
      4.22784335098467139e-1_wp, 9.22784335098467139e-1_wp, 1.25611766843180047_wp, &
      1.50611766843180047_wp, 1.70611766843180047_wp, 1.87278433509846714_wp,  &
      2.01564147795561000_wp, 2.14064147795561000_wp, 2.25175258906672111_wp,  &
      2.35175258906672111_wp, 2.44266167997581202_wp, 2.52599501330914535_wp,  &
      2.60291809023222227_wp, 2.67434666166079370_wp, 2.74101332832746037_wp,  &
      2.80351332832746037_wp, 2.86233685773922507_wp, 2.91789241329478063_wp,  &
      2.97052399224214905_wp, 3.02052399224214905_wp, 3.06814303986119667_wp,  &
      3.11359758531574212_wp, 3.15707584618530734_wp, 3.19874251285197401_wp,  &
      3.23874251285197401_wp, 3.27720405131351247_wp, 3.31424108835054951_wp,  &
      3.34995537406483522_wp, 3.38443813268552488_wp, 3.41777146601885821_wp,  &
      3.45002953053498724_wp, 3.48127953053498724_wp, 3.51158256083801755_wp,  &
      3.54099432554389990_wp, 3.56956575411532847_wp, 3.59734353189310625_wp,  &
      3.62437055892013327_wp, 3.65068634839381748_wp, 3.67632737403484313_wp,  &
      3.70132737403484313_wp, 3.72571761793728215_wp, 3.74952714174680596_wp,  &
      3.77278295570029433_wp, 3.79551022842756706_wp, 3.81773245064978928_wp,  &
      3.83947158108457189_wp, 3.86074817682925274_wp, 3.88158151016258607_wp,  &
      3.90198967342789220_wp, 3.92198967342789220_wp, 3.94159751656514710_wp,  &
      3.96082828579591633_wp, 3.97969621032421822_wp, 3.99821472884273674_wp,  &
      4.01639654702455492_wp, 4.03425368988169777_wp, 4.05179754953082058_wp,  &
      4.06903892884116541_wp, 4.08598808138353829_wp, 4.10265474805020496_wp,  &
      4.11904819067315578_wp, 4.13517722293122029_wp, 4.15105023880423617_wp,  &
      4.16667523880423617_wp, 4.18205985418885155_wp, 4.19721136934036670_wp,  &
      4.21213674247469506_wp, 4.22684262482763624_wp, 4.24133537845082464_wp,  &
      4.25562109273653893_wp, 4.26970559977879245_wp, 4.28359448866768134_wp,  &
      4.29729311880466764_wp, 4.31080663231818115_wp, 4.32413996565151449_wp,  &
      4.33729786038835659_wp, 4.35028487337536958_wp, 4.36310538619588240_wp,  &
      4.37576361404398366_wp, 4.38826361404398366_wp, 4.40060929305632934_wp,  &
      4.41280441500754886_wp, 4.42485260777863319_wp, 4.43675736968339510_wp,  &
      4.44852207556574804_wp, 4.46014998254249223_wp, 4.47164423541605544_wp,  &
      4.48300787177969181_wp, 4.49424382683587158_wp, 4.50535493794698269_wp,  &
      4.51634394893599368_wp, 4.52721351415338499_wp, 4.53796620232542800_wp,  &
      4.54860450019776842_wp, 4.55913081598724211_wp, 4.56954748265390877_wp,  &
      4.57985676100442424_wp, 4.59006084263707730_wp, 4.60016185273808740_wp /)
    !-----------------------------------------------------------------------
    ! COEFFICIENTS OF ASYMPTOTIC EXPANSION
    !-----------------------------------------------------------------------
    REAL(wp), PARAMETER :: b(6) = (/ &
      8.33333333333333333e-2_wp, -8.33333333333333333e-3_wp, &
      3.96825396825396825e-3_wp, -4.16666666666666666e-3_wp, &
      7.57575757575757576e-3_wp, -2.10927960927960928e-2_wp &
    /)

    IF (n <= 100) THEN
      fn_val = c(n)
      RETURN
    END IF
    wdtol = MAX(EPSILON(0.0_wp), 1.0e-18_wp)
    fn = n
    ax = 1.0_wp
    s = -0.5_wp / fn
    IF (ABS(s) > wdtol) THEN
      rfn2 = 1.0_wp / (fn*fn)
      DO k = 1, 6
        ax = ax * rfn2
        trm = -b(k) * ax
        IF (ABS(trm) < wdtol) EXIT
        s = s + trm
      END DO
    END IF

    fn_val = s + LOG(fn)
    RETURN
  END FUNCTION psixn

END MODULE calgo_683
