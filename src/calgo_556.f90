SUBROUTINE expint(x, n, kode, m, tol, en, ierr)                   EXP   10

! Code converted using TO_F90 by Alan Miller
! Date: 2025-05-27  Time: 01:35:20

!     WRITTEN BY D.E. AMOS, SANDIA LABORATORIES, ALBUQUERQUE, NM, 87185

!     REFERENCE
!         COMPUTATION OF EXPONENTIAL INTEGRALS BY D.E. AMOS, ACM
!         TRANS. MATH SOFTWARE, 1980

!     ABSTRACT
!         EXPINT COMPUTES M MEMBER SEQUENCES OF EXPONENTIAL INTEGRALS
!         E(N+K,X), K=0,1,...,M-1 FOR N.GE.1 AND X.GE.0.  THE POWER
!         SERIES IS IMPLEMENTED FOR X.LE.XCUT AND THE CONFLUENT
!         HYPERGEOMETRIC REPRESENTATION

!                     E(A,X) = EXP(-X)*(X**(A-1))*U(A,A,X)

!         IS COMPUTED FOR X.GT.XCUT. SINCE SEQUENCES ARE COMPUTED IN A
!         STABLE FASHION BY RECURRING AWAY FROM X, A IS SELECTED AS THE
!         INTEGER CLOSEST TO X WITHIN THE CONSTRAINT N.LE.A.LE.N+M-1.
!         FOR THE U COMPUTATION  A IS FURTHER MODIFIED TO BE THE
!         NEAREST EVEN INTEGER. INDICES ARE CARRIED FORWARD OR
!         BACKWARD BY THE TWO TERM RECURSION RELATION

!                     K*E(K+1,X) + X*E(K,X) = EXP(-X)

!         ONCE E(A,X) IS COMPUTED. THE U FUNCTION IS COMPUTED BY MEANS
!         OF THE BACKWARD RECURSIVE MILLER ALGORITHM APPLIED TO THE
!         THREE TERM CONTIGUOUS RELATION FOR U(A+K,A,X), K=0,1,...
!         THIS PRODUCES ACCURATE RATIOS AND DETERMINES U(A+K,A,X),AND
!         HENCE E(A,X), TO WITHIN A MULTIPLICATIVE CONSTANT C.
!         ANOTHER CONTIGUOUS RELATION APPLIED TO C*U(A,A,X) AND
!         C*U(A+1,A,X) GETS C*U(A+1,A+1,X), A QUANTITY PROPORTIONAL TO
!         E(A+1,X). THE NORMALIZING CONSTANT C IS OBTAINED FROM THE
!         TWO TERM RECURSION RELATION ABOVE WITH K=A.

!         MACHINE DEPENDENT PARAMETERS - XCUT, XLIM, ETOL, EULER, DIGAM

!         EXPINT WRITES ERROR DIAGNOSTICS TO LOGICAL UNIT 3

!     DESCRIPTION OF ARGUMENTS

!         INPUT
!           X       X.GT.0.0 FOR N=1 AND  X.GE.0.0 FOR N.GE.2
!           N       ORDER OF THE FIRST MEMBER OF THE SEQUENCE, N.GE.1
!           KODE    A SELECTION PARAMETER FOR SCALED VALUES
!                   KODE=1   RETURNS        E(N+K,X), K=0,1,...,M-1.
!                       =2   RETURNS EXP(X)*E(N+K,X), K=0,1,...,M-1.
!           M       NUMBER OF EXPONENTIAL INTEGRALS IN THE SEQUENCE,
!                   M.GE.1
!           TOL     RELATIVE ACCURACY WANTED, ETOL.LE.TOL.LE.0.1
!                   ETOL=1.E-12

!         OUTPUT
!           EN      A VECTOR OF DIMENSION AT LEAST M CONTAINING VALUES
!                   EN(K) = E(N+K-1,X) OR EXP(X)*E(N+K-1,X), K=1,M
!                   DEPENDING ON KODE
!           IERR    UNDERFLOW INDICATOR
!                   IERR=0   A NORMAL RETURN
!                       =1   X EXCEEDS XLIM AND AN UNDERFLOW OCCURS.
!                            EN(K)=0.0 , K=1,M RETURNED ON KODE=1
!                            XLIM=667.

!     ERROR CONDITIONS
!         AN IMPROPER INPUT PARAMETER IS A FATAL ERROR
!         UNDERFLOW IS A NON FATAL ERROR. ZERO ANSWERS ARE RETURNED.


REAL, INTENT(IN OUT)                     :: x
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: kode
INTEGER, INTENT(IN)                      :: m
REAL, INTENT(IN)                         :: tol
REAL, INTENT(OUT)                        :: en(1)
INTEGER, INTENT(OUT)                     :: ierr
DIMENSION  a(99), b(99), y(2)

DATA xcut, xlim, etol /2.0E0,667.0E0,1.0E-12/
DATA euler /-5.77215664901533E-01/
DATA lun /3/

IF (n < 1) GO TO 260
IF (kode < 1 .OR. kode > 2) GO TO 270
IF (m < 1) GO TO 280
IF (tol < etol .OR. tol > 0.1E0) GO TO 290

ierr = 0
IF (x > xcut) GO TO 100
IF (x < 0.0E0) GO TO 300
IF (x == 0.0E0 .AND. n == 1) GO TO 310
IF (x == 0.0E0 .AND. n > 1) GO TO 80

!     SERIES FOR E(N,X) FOR X.LE.XCUT

ix = INT(x+0.5E0)
!     ICASE=1 MEANS INTEGER CLOSEST TO X IS 2 AND N=1
!     ICASE=2 MEANS INTEGER CLOSEST TO X IS 0,1, OR 2 AND N.GE.2
icase = 2
IF (ix > n) icase = 1
nm = n - icase + 1
nd = nm + 1
ind = 3 - icase
mu = m - ind
ml = 1
ks = nd
fnm = FLOAT(nm)
s = 0.0E0
xtol = 3.0E0*tol
IF (nd == 1) GO TO 10
xtol = 0.3333E0*tol
s = 1.0E0/fnm
10 CONTINUE
aa = 1.0E0
ak = 1.0E0
DO  i=1,35
  aa = -aa*x/ak
  IF (i == nm) GO TO 30
  s = s - aa/(ak-fnm)
  IF (ABS(aa) <= xtol*ABS(s)) GO TO 20
  ak = ak + 1.0E0
  CYCLE
  20   CONTINUE
  IF (i < 2) GO TO 40
  IF (nd-2 > i .OR. i > nd-1) GO TO 60
  ak = ak + 1.0E0
  CYCLE
  30   s = s + aa*(-ALOG(x)+digam(nd))
  xtol = 3.0E0*tol
  40   ak = ak + 1.0E0
END DO
GO TO 320
60 IF (nd == 1) s = s + (-ALOG(x)+euler)
IF (kode == 2) s = s*EXP(x)
en(1) = s
emx = 1.0E0
IF (m == 1) GO TO 70
en(ind) = s
aa = FLOAT(ks)
IF (kode == 1) emx = EXP(-x)
SELECT CASE ( icase )
  CASE (    1)
    GO TO 220
  CASE (    2)
    GO TO  240
END SELECT
70 IF (icase == 2) RETURN
IF (kode == 1) emx = EXP(-x)
en(1) = (emx-s)/x
RETURN
80 CONTINUE
DO  i=1,m
  en(i) = 1.0E0/FLOAT(n+i-2)
END DO
RETURN

!     BACKWARD RECURSIVE MILLER ALGORITHM FOR
!              E(N,X)=EXP(-X)*(X**(N-1))*U(N,N,X)
!     WITH RECURSION AWAY FROM N=INTEGER CLOSEST TO X.
!     U(A,B,X) IS THE SECOND CONFLUENT HYPERGEOMETRIC FUNCTION

100 CONTINUE
emx = 1.0E0
IF (kode == 2) GO TO 130
IF (x <= xlim) GO TO 120
ierr = 1
DO  i=1,m
  en(i) = 0.0E0
END DO
RETURN
120 emx = EXP(-x)
130 CONTINUE
ix = INT(x+0.5E0)
kn = n + m - 1
IF (kn <= ix) GO TO 140
IF (n < ix .AND. ix < kn) GO TO 170
IF (n >= ix) GO TO 160
GO TO 340
140 icase = 1
ks = kn
ml = m - 1
mu = -1
ind = m
IF (kn > 1) GO TO 180
150 ks = 2
icase = 3
GO TO 180
160 icase = 2
ind = 1
ks = n
mu = m - 1
IF (n > 1) GO TO 180
IF (kn == 1) GO TO 150
ix = 2
170 icase = 1
ks = ix
ml = ix - n
ind = ml + 1
mu = kn - ix
180 CONTINUE
ik = ks/2
ah = FLOAT(ik)
jset = 1 + ks - (ik+ik)
!     START COMPUTATION FOR
!              EN(IND) = C*U( A , A ,X)    JSET=1
!              EN(IND) = C*U(A+1,A+1,X)    JSET=2
!     FOR AN EVEN INTEGER A.
ic = 0
aa = ah + ah
aams = aa - 1.0E0
aams = aams*aams
tx = x + x
fx = tx + tx
ak = ah
xtol = tol
IF (tol <= 1.0E-3) xtol = 20.0E0*tol
ct = aams + fx*ah
em = (ah+1.0E0)/((x+aa)*xtol*SQRT(ct))
bk = aa
cc = ah*ah
!     FORWARD RECURSION FOR P(IC),P(IC+1) AND INDEX IC FOR BACKWARD
!     RECURSION
p1 = 0.0E0
p2 = 1.0E0
190 CONTINUE
IF (ic == 99) GO TO 330
ic = ic + 1
ak = ak + 1.0E0
at = bk/(bk+ak+cc+FLOAT(ic))
bk = bk + ak + ak
a(ic) = at
bt = (ak+ak+x)/(ak+1.0E0)
b(ic) = bt
pt = p2
p2 = bt*p2 - at*p1
p1 = pt
ct = ct + fx
em = em*at*(1.0E0-tx/ct)
IF (em*(ak+1.0E0) > p1*p1) GO TO 190
ict = ic
kk = ic + 1
bt = tx/(ct+fx)
y2 = (bk/(bk+cc+FLOAT(kk)))*(p1/p2)*(1.0E0-bt+0.375E0*bt*bt)
y1 = 1.0E0
!     BACKWARD RECURRENCE FOR
!              Y1=             C*U( A ,A,X)
!              Y2= C*(A/(1+A/2))*U(A+1,A,X)
DO  k=1,ict
  kk = kk - 1
  yt = y1
  y1 = (b(kk)*y1-y2)/a(kk)
  y2 = yt
END DO
!     THE CONTIGUOUS RELATION
!              X*U(B,C+1,X)=(C-B)*U(B,C,X)+U(B-1,C,X)
!     WITH  B=A+1 , C=A IS USED FOR
!              Y(2) = C * U(A+1,A+1,X)
!     X IS INCORPORATED INTO THE NORMALIZING RELATION FOR CNORM.
pt=y2/y1
cnorm=1.0E0-pt*(ah+1.0E0)/aa
y(1)=1.0E0/(cnorm*aa+x)
y(2)=cnorm*y(1)
IF (icase == 3) GO TO 210
en(ind) =   emx*y(jset)
IF (m == 1) RETURN
aa = FLOAT(ks)
SELECT CASE ( icase )
  CASE (    1)
    GO TO 220
  CASE (    2)
    GO TO  240
END SELECT

!     RECURSION SECTION  N*E(N+1,X) + X*E(N,X)=EMX

210 en(1) = emx*(1.0E0-y(1))/x
RETURN
220 k = ind - 1
DO  i=1,ml
  aa = aa - 1.0E0
  en(k) = (emx-aa*en(k+1))/x
  k = k - 1
END DO
IF (mu <= 0) RETURN
aa = FLOAT(ks)
240 k = ind
DO  i=1,mu
  en(k+1) = (emx-x*en(k))/aa
  aa = aa + 1.0E0
  k = k + 1
END DO
RETURN


260 WRITE (lun,99999)
RETURN
270 WRITE (lun,99998)
RETURN
280 WRITE (lun,99997)
RETURN
290 WRITE (lun,99996)
RETURN
300 WRITE (lun,99995)
RETURN
310 WRITE (lun,99994)
RETURN
320 WRITE (lun,99993)
RETURN
330 WRITE (lun,99992)
RETURN
340 WRITE (lun,99991)
RETURN
99999 FORMAT (32H in expint, n NOT greater than 0)
99998 FORMAT (27H in expint, kode NOT 1 OR 2)
99997 FORMAT (32H in expint, m NOT greater than 0)
99996 FORMAT (33H in expint, tol NOT within limits)
99995 FORMAT (37H in expint, x is NOT zero OR positive)
99994 FORMAT (46H in expint, the exponential integral is NOT de,  &
    21HFINED for x=0 AND n=1)
99993 FORMAT (46H in expint, relative error test for series ter,  &
    28HMINATION NOT met in 36 terms)
99992 FORMAT (46H in expint, termination test for miller algori,  &
    23HTHM NOT met in 99 steps)
99991 FORMAT (46H in expint, an error in placing INT(x+0.5) wit,  &
    47HH respect TO n AND n+m-1 occurred for x > xcut)
END SUBROUTINE expint

FUNCTION digam(n)                                                 dig   10

!     THIS SUBROUTINE RETURNS VALUES OF PSI(X)=DERIVATIVE OF LOG
!     GAMMA(X), X.GT.0.0 AT INTEGER ARGUMENTS. A TABLE LOOK-UP IS
!     PERFORMED FOR N.LE.100, AND THE ASYMPTOTIC EXPANSION IS
!     EVALUATED FOR N.GT.100.


INTEGER, INTENT(IN)                      :: n
DIMENSION b(4), c(100), c1(32), c2(27), c3(22), c4(19)
EQUIVALENCE (c(1),c1(1))
EQUIVALENCE (c(33),c2(1))
EQUIVALENCE (c(60),c3(1))
EQUIVALENCE (c(82),c4(1))

DATA c1 /-5.7721566490153E-01,4.22784335098467E-01,  &
    9.22784335098467E-01,1.25611766843180E+00,1.50611766843180E+00,  &
    1.70611766843180E+00,1.87278433509847E+00,2.01564147795561E+00,  &
    2.14064147795561E+00,2.25175258906672E+00,2.35175258906672E+00,  &
    2.44266167997581E+00,2.52599501330915E+00,2.60291809023222E+00,  &
    2.67434666166079E+00,2.74101332832746E+00,2.80351332832746E+00,  &
    2.86233685773923E+00,2.91789241329478E+00,2.97052399224215E+00,  &
    3.02052399224215E+00,3.06814303986120E+00,3.11359758531574E+00,  &
    3.15707584618531E+00,3.19874251285197E+00,3.23874251285197E+00,  &
    3.27720405131351E+00,3.31424108835055E+00,3.34995537406484E+00,  &
    3.38443813268552E+00,3.41777146601886E+00,3.45002953053499E+00/
DATA c2 /3.48127953053499E+00,3.51158256083802E+00,  &
    3.54099432554390E+00,3.56956575411533E+00,3.59734353189311E+00,  &
    3.62437055892013E+00,3.65068634839382E+00,3.67632737403484E+00,  &
    3.70132737403484E+00,3.72571761793728E+00,3.74952714174681E+00,  &
    3.77278295570029E+00,3.79551022842757E+00,3.81773245064979E+00,  &
    3.83947158108457E+00,3.86074817682925E+00,3.88158151016259E+00,  &
    3.90198967342789E+00,3.92198967342789E+00,3.94159751656515E+00,  &
    3.96082828579592E+00,3.97969621032422E+00,3.99821472884274E+00,  &
    4.01639654702455E+00,4.03425368988170E+00,4.05179754953082E+00,  &
    4.06903892884117E+00/
DATA c3 /4.08598808138354E+00,4.10265474805020E+00,  &
    4.11904819067316E+00,4.13517722293122E+00,4.15105023880424E+00,  &
    4.16667523880424E+00,4.18205985418885E+00,4.19721136934037E+00,  &
    4.21213674247470E+00,4.22684262482764E+00,4.24133537845082E+00,  &
    4.25562109273654E+00,4.26970559977879E+00,4.28359448866768E+00,  &
    4.29729311880467E+00,4.31080663231818E+00,4.32413996565151E+00,  &
    4.33729786038836E+00,4.35028487337537E+00,4.36310538619588E+00,  &
    4.37576361404398E+00,4.38826361404398E+00/
DATA c4 /4.40060929305633E+00,4.41280441500755E+00,  &
    4.42485260777863E+00,4.43675736968340E+00,4.44852207556575E+00,  &
    4.46014998254249E+00,4.47164423541606E+00,4.48300787177969E+00,  &
    4.49424382683587E+00,4.50535493794698E+00,4.51634394893599E+00,  &
    4.52721351415338E+00,4.53796620232543E+00,4.54860450019777E+00,  &
    4.55913081598724E+00,4.56954748265391E+00,4.57985676100442E+00,  &
    4.59006084263708E+00,4.60016185273809E+00/

DATA b /1.66666666666667E-01,-3.33333333333333E-02,  &
    2.38095238095238E-02,-3.33333333333333E-02/

IF (n > 100) GO TO 10
digam = c(n)
RETURN
10 fn = n
ax = 1.0E0
ak = 2.0E0
s = -0.5E0/fn
IF (fn > 1.e+8) GO TO 30
fn2 = fn*fn
DO  k=1,3
  ax = ax*fn2
  s = s - b(k)/(ax*ak)
  ak = ak + 2.0E0
END DO
30 CONTINUE
digam = s + ALOG(fn)
RETURN
END FUNCTION digam
!     PROGRAM TSTEXP(INPUT,OUTPUT,TAPE3=OUTPUT)                         00000010
!                                                                       00000020
!     PROGRAM TO TEST SUBROUTINE EXPINT AGAINST AN ADAPTIVE QUADRATURE. 00000030
!     PARAMETER VALUES ARE PRINTED AND, IN THE EVENT THAT THE RELATIVE  00000040
!     ERROR TEST IS NOT SATISFIED, X, ERROR, N, AND KODE ARE ALSO       00000050
!     PRINTED. AN OUTPUT WITH ONLY PARAMETER VALUES INDICATES THAT ALL  00000060
!     TESTS WERE PASSED. GAUS8 COMPUTES THE QUADRATURES.                00000070
!                                                                       00000080
DIMENSION xtol(10), en(50), ev(50)                                00000090
iout = 3                                                          00000100
nl = 1                                                            00000110
nu = 16                                                           00000120
ninc = 5                                                          00000130
ml = 1                                                            00000140
mu = 25                                                           00000150
minc = 8                                                          00000160
km = 5                                                            00000170
jl = 1                                                            00000180
ju = 40                                                           00000190
jinc = 3                                                          00000200
xtol(1) = 1.0E-2                                                  00000210
DO  i=2,3                                                       00000220
  xtol(i) = xtol(i-1)*1.0E-3                                      00000230
END DO
DO  it=1,3                                                      00000250
  tol = xtol(it)                                                  00000260
  tola = AMAX1(1.0E-12,tol/10.0E0)                                00000270
  btol = tol                                                      00000280
  WRITE (iout,99999) tol                                          00000290
  DO  m=ml,mu,minc                                              00000300
    WRITE (iout,99998) m                                          00000310
    DO  n=nl,nu,ninc                                            00000320
      WRITE (iout,99997) n                                        00000330
      DO  j=jl,ju,jinc                                          00000340
        x = FLOAT(j-1)/5.0E0                                      00000350
        ex = EXP(-x)                                              00000360
        IF (x == 0. .AND. n == 1) CYCLE
        CALL expint(x, n, 1, m, tol, en, ierr)                    00000380
        CALL expint(x, n, 2, m, tol, ev, ierr)                    00000390
        DO  k=1,m,km                                            00000400
          IF (x > 0.) GO TO 20                                   00000410
          IF (n+k == 2) CYCLE
          y = 1.0E0/FLOAT(n+k-2)                                  00000430
          yy = y                                                  00000440
          GO TO 30                                                00000450
          20           CONTINUE                                                00000460
          nn = n + k - 1                                          00000470
          yy = eint(nn,x,tola,2)                                  00000480
          y = yy*ex                                               00000490
          30           CONTINUE                                                00000500
          er = ABS((y-en(k))/y)                                   00000510
          kode = 1                                                00000520
          IF (er > btol) WRITE (iout,99996) x, er, nn, kode      00000530
          kode = 2                                                00000540
          ERR = ABS((yy-ev(k))/yy)                                00000550
          IF (ERR > btol) WRITE (iout,99996) x, ERR, nn, kode    00000560
        END DO
      END DO
    END DO
  END DO
END DO
STOP                                                              00000620
99999 FORMAT (1H0, 5H tol=, e15.4/)                                     00000630
99998 FORMAT (1H0, 2HM=, i5/)                                           00000640
99997 FORMAT (3X, 2HN=, i5)                                             00000650
99996 FORMAT (2E15.6, 2I5)                                              00000660
END                                                               00000670

FUNCTION eint(n, x, tol, kode)                                    00000680

INTEGER, INTENT(IN)                      :: n
REAL, INTENT(IN)                         :: x
REAL, INTENT(IN)                         :: tol
INTEGER, INTENT(IN OUT)                  :: kode
COMMON /geint/ xx, fn
EXTERNAL feint

xx = x
fn = n
sig = 1.0E0
s = 0.0E0
tola = tol
b = x
10 CONTINUE
a = b
rel = tol
b = b + sig
CALL gaus8(feint, a, b, rel, ans, ierr)
s = s + ans
IF (ABS(ans) < s*tola) GO TO 20
GO TO 10
20 eint = s*EXP((fn-1.0E0)*ALOG(x)-FLOAT(2-kode)*x)
RETURN
END FUNCTION eint

FUNCTION feint(t)                                                 00000880


REAL, INTENT(IN)                         :: t
COMMON /geint/ xx, fn
feint = EXP(-t+xx-fn*ALOG(t))
RETURN
END FUNCTION feint

SUBROUTINE gaus8  (fun,a,b,ERR,ans,ierr)                          00000930

!     BY RONDALL E JONES, SANDIA LABORATORIES
!     SALIENT FEATURES -- INTERVAL BISECTION, COMBINED RELATIVE/ABSOLUTE
!     ERROR CONTROL, COMPUTED MAXIMUM REFINEMENT LEVEL WHEN A IS
!     CLOSE TO B.

!     ABSTRACT
!        GAUS8 INTEGRATES REAL FUNCTIONS OF ONE VARIABLE OVER FINITE
!        INTERVALS, USING AN ADAPTIVE 8-POINT LEGENDRE-GAUSS ALGORITHM.
!        GAUS8 IS INTENDED PRIMARILY FOR HIGH ACCURACY INTEGRATION
!        OR INTEGRATION OF SMOOTH FUNCTIONS.  FOR LOWER ACCURACY
!        INTEGRATION OF FUNCTIONS WHICH ARE NOT VERY SMOOTH,
!        EITHER QNC3 OR QNC7 MAY BE MORE EFFICIENT.

!     DESCRIPTION OF ARGUMENTS

!        INPUT--
!        FUN - NAME OF EXTERNAL FUNCTION TO BE INTEGRATED.  THIS NAME
!              MUST BE IN AN EXTERNAL STATEMENT IN THE CALLING PROGRAM.
!              FUN MUST BE A FUNCTION OF ONE REAL ARGUMENT.  THE VALUE
!              OF THE ARGUMENT TO FUN IS THE VARIABLE OF INTEGRATION
!              WHICH RANGES FROM A TO B.
!        A   - LOWER LIMIT OF INTEGRAL
!        B   - UPPER LIMIT OF INTEGRAL (MAY BE LESS THAN A)
!        ERR - IS A REQUESTED ERROR TOLERANCE.  NORMALLY PICK A VALUE OF
!              ABS(ERR).LT.1.E-3.  ANS WILL NORMALLY HAVE NO MORE ERROR
!              THAN ABS(ERR) TIMES THE INTEGRAL OF THE ABSOLUTE VALUE
!              OF FUN(X).  USUALLY, SMALLER VALUES FOR ERR YIELD
!              MORE ACCURACY AND REQUIRE MORE FUNCTION EVALUATIONS.
!              A NEGATIVE VALUE FOR ERR CAUSES AN ESTIMATE OF THE
!              ABSOLUTE ERROR IN ANS TO BE RETURNED IN ERR.

!        OUTPUT--
!        ERR - WILL BE AN ESTIMATE OF THE ERROR IN ANS IF THE INPUT
!              VALUE OF ERR WAS NEGATIVE.  THE ESTIMATED ERROR IS SOLELY
!              FOR INFORMATION TO THE USER AND SHOULD NOT BE USED AS
!              A CORRECTION TO THE COMPUTED INTEGRAL.
!        ANS - COMPUTED VALUE OF INTEGRAL
!        IERR- A STATUS CODE
!            --NORMAL CODES
!               1 ANS MOST LIKELY MEETS REQUESTED ERROR TOLERANCE,
!                 OR A=B.
!              -1 A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMAL
!                 INTEGRATION.  ANS IS SET TO ZERO.
!            --ABNORMAL CODE
!               2 ANS PROBABLY DOES NOT MEET REQUESTED ERROR TOLERANCE.



!     GAUS8  USES SUBROUTINES ERRCHK, ERRGET, ERRPRT, ERXSET, ERSTGT
!     COMPILE DECKS GAUS8, ERRCHK


REAL, INTENT(IN)                         :: fun
REAL, INTENT(IN)                         :: a
REAL, INTENT(IN OUT)                     :: b
REAL, INTENT(OUT)                        :: ERR
REAL, INTENT(OUT)                        :: ans
INTEGER, INTENT(OUT)                     :: ierr
DIMENSION aa(30),hh(30),lr(30),vl(30),gr(30)
DATA x1,x2,x3,x4/0.183434642495650 , 0.525532409916329 ,  &
    0.796666477413627 , 0.960289856497536 /
DATA w1,w2,w3,w4/0.362683783378362 , 0.313706645877887 ,  &
    0.222381034453374 , 0.101228536290376 /
DATA sq2/1.41421356/,icall/0/
DATA nlmn/1/,nlmx/30/,kmx/5000/,kml/6/,nbits/48/

g8(x,h) = h*( (w1*(fun(x-x1*h)+fun(x+x1*h)) +w2*(fun(x-x2*h)+fun(x+x2*h)))  &
    +(w3*(fun(x-x3*h)+fun(x+x3*h)) +w4*(fun(x-x4*h)+fun(x+x4*h))) )

!     INITIALIZE

IF(icall /= 0)CALL errchk(-71,71H*****gaus8 called recursively.  r  &
    ecursive calls are illegal in fortran. )
icall = 1
ans = 0.0
ierr = 1
ce = 0.0
IF (a == b) GO TO 35
lmx = nlmx
lmn = nlmn
IF (b == 0.0) GO TO 4
IF (SIGN(1.0,b)*a <= 0.0) GO TO 4
c = ABS(1.0-a/b)
IF (c > 0.1) GO TO 4
IF (c <= 0.0) GO TO 35
nib = 0.5-ALOG(c)/ALOG(2.0)
lmx = MIN0(nlmx , nbits-nib-7)
IF (lmx < 1) GO TO 32
lmn = MIN0(lmn,lmx)
4 tol = AMAX1(ABS(ERR),2.0**(5-nbits))/2.0
IF (ERR == 0.0) tol = 0.5E-6
eps = tol
hh(1) = (b-a)/4.0
aa(1) = a
lr(1) = 1
l = 1
est = g8(aa(l)+2.0*hh(l),2.0*hh(l))
k = 8
area = ABS(est)
ef = 0.5
mxl = 0

!     COMPUTE REFINED ESTIMATES, ESTIMATE THE ERROR, ETC.

5 gl = g8(aa(l)+hh(l),hh(l))
gr(l) = g8(aa(l)+3.0*hh(l),hh(l))
k = k+16
area = area+(ABS(gl)+ABS(gr(l))-ABS(est))
!     IF (L.LT.LMN) GO TO 11
glr = gl+gr(l)
ee = ABS(est-glr)*ef
ae = AMAX1(eps*area,tol*ABS(glr))
IF (ee-ae > 0.0) THEN
  GO TO    10
ELSE
  GO TO     8
END IF
7 mxl = 1
8 ce = ce + (est-glr)
IF (lr(l) > 0) THEN
  GO TO    20
ELSE
  GO TO    15
END IF

!     CONSIDER THE LEFT HALF OF THIS LEVEL

10 IF (k > kmx) lmx = kml
IF (l >= lmx) GO TO 7
11 l = l+1
eps = eps*0.5
ef = ef/sq2
hh(l) = hh(l-1)*0.5
lr(l) = -1
aa(l) = aa(l-1)
est = gl
GO TO 5

!     PROCEED TO RIGHT HALF AT THIS LEVEL

15 vl(l) = glr
16 est = gr(l-1)
lr(l) = 1
aa(l) = aa(l)+4.0*hh(l)
GO TO 5

!     RETURN ONE LEVEL

20 vr = glr
22 IF (l <= 1) GO TO 30
l = l-1
eps = eps*2.0
ef = ef*sq2
IF (lr(l) > 0) THEN
  GO TO    26
END IF
24 vl(l) = vl(l+1)+vr
GO TO 16
26 vr = vl(l+1)+vr
GO TO 22

!      EXIT

30 ans = vr
IF ((mxl == 0).OR.(ABS(ce) <= 2.0*tol*area)) GO TO 35
ierr = 2
CALL errchk(51,51HIN gaus8 , ans is probably insufficiently accura te.)
GO TO 35
32 ierr =-1
CALL onechk(-70,70HTHE following temporary informative diagnostic  &
    will appear only once. )
callonechk(-102,102HIN gaus8 , a AND b are too nearly equal TO all  &
    ow normal integration.  ans is set TO zero, AND ierr=-1.)
35 icall = 0
IF (ERR < 0.0) ERR = ce
RETURN
END SUBROUTINE gaus8

SUBROUTINE errchk(nchars,narray)                                  00002570

!     SANDIA MATHEMATICAL PROGRAM LIBRARY
!     APPLIED MATHEMATICS DIVISION 2642
!     SANDIA LABORATORIES
!     ALBUQUERQUE, NEW MEXICO 87115

!     SIMPLIFIED VERSION FOR STAND-ALONE USE.     APRIL 1977

!     ABSTRACT
!         THE ROUTINES ERRCHK, ERXSET, AND ERRGET TOGETHER PROVIDE
!         A UNIFORM METHOD WITH SEVERAL OPTIONS FOR THE PROCESSING
!         OF DIAGNOSTICS AND WARNING MESSAGES WHICH ORIGINATE
!         IN THE MATHEMATICAL PROGRAM LIBRARY ROUTINES.
!         ERRCHK IS THE CENTRAL ROUTINE, WHICH ACTUALLY PROCESSES
!         MESSAGES.

!     DESCRIPTION OF ARGUMENTS
!         NCHARS - NUMBER OF CHARACTERS IN HOLLERITH MESSAGE.
!                  IF NCHARS IS NEGATED, ERRCHK WILL UNCONDITIONALLY
!                  PRINT THE MESSAGE AND STOP EXECUTION.  OTHERWISE,
!                  THE BEHAVIOR OF ERRCHK MAY BE CONTROLLED BY
!                  AN APPROPRIATE CALL TO ERXSET.
!         NARRAY - NAME OF ARRAY OR VARIABLE CONTAINING THE MESSAGE,
!                  OR ELSE A LITERAL HOLLERITH CONSTANT CONTAINING
!                  THE MESSAGE.  BY CONVENTION, ALL MESSAGES SHOULD
!                  BEGIN WITH *IN SUBNAM, ...*, WHERE SUBNAM IS THE
!                  NAME OF THE ROUTINE CALLING ERRCHK.

!     EXAMPLES
!         1. TO ALLOW CONTROL BY CALLING ERXSET, USE
!            CALL ERRCHK(30,30HIN QUAD, INVALID VALUE OF ERR.)
!         2. TO UNCONDITIONALLY PRINT A MESSAGE AND STOP EXECUTION, USE
!            CALL ERRCHK(-30,30HIN QUAD, INVALID VALUE OF ERR.)



!     ERRCHK USES SUBROUTINES ERRGET, ERRPRT, ERXSET, ERSTGT
!     COMPILE DECKS ERRCHK



INTEGER, INTENT(IN OUT)                  :: nchars
INTEGER, INTENT(IN OUT)                  :: narray(14)


iout=6
CALL errget(nf,nt)
!     IF ERRCHK WAS CALLED WITH NEGATIVE CHARACTER COUNT, SET FATAL FLAG
IF (nchars < 0) nf = -1
!     IF MESSAGES ARE TO BE SUPPRESSED, RETURN
IF (nf == 0) RETURN
!     IF CHARACTER COUNT IS INVALID, STOP
!     IF (NCHARS.EQ.0) PRINT 5
IF (nchars == 0) WRITE (iout,5)
5 FORMAT(/31H errchk was called incorrectly.)
IF (nchars == 0) STOP
!     PRINT MESSAGE
CALL errprt(IABS(nchars),narray)
!     IF LAST MESSAGE, SAY SO
!     IF (NF.EQ.1) PRINT 10
IF (nf == 1) WRITE (iout,10)
10 FORMAT (30H errchk message limit reached.)
!     PRINT TRACE-BACK IF ASKED TO
!     IF ((NT.GT.0).OR.(NF.LT.0)) CALL SYSTEM ROUTINE FOR TRACEBACK
!     DECREMENT MESSAGE COUNT
IF (nf > 0) nf = nf-1
CALL erxset(nf,nt)
!     IF ALL IS WELL, RETURN
IF (nf >= 0) RETURN
!     IF THIS MESSAGE IS SUPPRESSABLE BY AN ERXSET CALL,
!     THEN EXPLAIN ERXSET USAGE.
!     IF (NCHARS.GT.0) PRINT 15
IF (nchars > 0) WRITE (iout,15)
15 FORMAT (/13H *** note ***  &
    /53H TO make the error message printed above be nonfatal,  &
    /39H OR TO suppress the message completely,  &
    /37H insert an appropriate CALL TO erxset  &
    ,30H at the start of your PROGRAM.  &
    /62H for example, TO PRINT up TO 10 nonfatal warning messages, use  &
    /27H          CALL erxset(10,0)    )
!     PRINT 20
WRITE (iout,20)
20 FORMAT (/28H PROGRAM abort due TO error.)
STOP
END SUBROUTINE errchk

SUBROUTINE onechk(nchars,narray)                                  00003390

!     ABSTRACT
!         ONECHK IS A COMPANION ROUTINE OF ERRCHK.  IT IS CALLED
!         JUST LIKE ERRCHK, AND MESSAGES FROM IT MAY BE SUPPRESSED
!         BY AN APPROPRIATE CALL TO ERXSET.  IT DIFFERS FROM ERRCHK
!         IN THAT EACH CALL TO ONECHK WILL PRODUCE NO MORE THAN ONE
!         PRINTED MESSAGE, REGARDLESS OF HOW MANY TIMES THAT CALL IS
!         EXECUTED, AND ONECHK NEVER TERMINATES EXECUTION.
!         ITS PURPOSE IS TO PROVIDE ONE-TIME-ONLY INFORMATIVE
!         DIAGNOSTICS.

!     DESCRIPTION OF ARGUMENTS
!         NCHARS - NUMBER OF CHARACTERS IN THE MESSAGE.
!                  IF NEGATED, THE MESSAGE WILL BE PRINTED (ONCE) EVEN
!                  IF NFATAL HAS BEEN SET TO 0 (SEE ERXSET).
!         NARRAY - SAME AS IN ERRCHK



!     ONECHK USES SUBROUTINES ERRGET, ERRPRT, ERXSET, ERSTGT
!     COMPILE DECKS ERRCHK


INTEGER, INTENT(IN OUT)                  :: nchars
INTEGER, INTENT(OUT)                     :: narray(14)

DATA nflag/4H.$,*/

IF (narray(1) == nflag) RETURN
CALL errget(nf,nt)
IF ((nf == 0).AND.(nchars > 0)) RETURN
CALL errprt (59,59HTHE following informative diagnostic will appea  &
    r only once.)
CALL errprt(IABS(nchars),narray)
IF (nf > 0) nf = nf-1
CALL erxset(nf,nt)
narray(1) = nflag
RETURN
END SUBROUTINE onechk

SUBROUTINE errprt(nchars,narray)                                  00003750

!     UTILITY ROUTINE TO SIMPLY PRINT THE HOLLERITH MESSAGE IN NARRAY,
!     WHOSE LENGTH IS NCHARS CHARACTERS.



INTEGER, INTENT(IN OUT)                  :: nchars
INTEGER, INTENT(IN OUT)                  :: narray(14)


!     NOTE - NCH MUST BE THE NUMBER OF HOLLERITH CHARACTERS STORED
!     PER WORD.  IF NCH IS CHANGED, FORMAT 1 MUST ALSO BE
!     CHANGED CORRESPONDINGLY.

iout=6
nch = 10
!     FOR LINE PRINTERS, USE
1 FORMAT (1X,13A10)
!     FOR DATA TERMINALS, USE
!   1 FORMAT (1X,7A10)
nwords = (nchars+nch-1)/nch
!     PRINT 1,(NARRAY(I),I=1,NWORDS)
WRITE (iout,1) (narray(i),i=1,nwords)
RETURN
END SUBROUTINE errprt

SUBROUTINE erxset(nfatal,ntrace)                                  00003970

!     ABSTRACT
!         ERXSET IS A COMPANION ROUTINE TO SUBROUTINE ERRCHK.
!         ERXSET ASSIGNS THE VALUES OF NFATAL AND NTRACE RESPECTIVELY
!         TO NF AND NT IN COMMON BLOCK MLBLK0 THEREBY SPECIFYING THE
!         STATE OF THE OPTIONS WHICH CONTROL THE EXECUTION OF ERRCHK.

!     DESCRIPTION OF ARGUMENTS
!         BOTH ARGUMENTS ARE INPUT ARGUMENTS OF DATA TYPE INTEGER.
!         NFATAL - IS A FATAL-ERROR / MESSAGE-LIMIT FLAG. A NEGATIVE
!                  VALUE DENOTES THAT DETECTED DIFFICULTIES ARE TO BE
!                  TREATED AS FATAL ERRORS.  NONNEGATIVE MEANS NONFATAL.
!                  A NONNEGATIVE VALUE IS THE MAXIMUM NUMBER OF NONFATAL
!                  WARNING MESSAGES WHICH WILL BE PRINTED BY ERRCHK,
!                  AFTER WHICH NONFATAL MESSAGES WILL NOT BE PRINTED.
!                  (DEFAULT VALUE IS -1.)
!         NTRACE - .GE.1 WILL CAUSE A TRACE-BACK TO BE GIVEN,
!                        IF THIS FEATURE IS IMPLEMENTED ON THIS SYSTEM.
!                  .LE.0 WILL SUPPRESS ANY TRACE-BACK, EXCEPT FOR
!                        CASES WHEN EXECUTION IS TERMINATED.
!                  (DEFAULT VALUE IS 0.)

!         *NOTE* -- SOME CALLS TO ERRCHK WILL CAUSE UNCONDITIONAL
!         TERMINATION OF EXECUTION.  ERXSET HAS NO EFFECT ON SUCH CALLS.

!     EXAMPLES
!         1. TO PRINT UP TO 100 MESSAGES AS NONFATAL WARNINGS USE
!            CALL ERXSET(100,0)
!         2. TO SUPPRESS ALL MATHLIB WARNING MESSAGES USE
!            CALL ERXSET(0,0)



!     ERXSET USES SUBROUTINES ERSTGT
!     COMPILE DECKS ERRCHK

CALL erstgt(0,nfatal,ntrace)
RETURN
END SUBROUTINE erxset

SUBROUTINE errget(nfatal,ntrace)                                  00004370

!     ABSTRACT
!         ERRGET IS A COMPANION ROUTINE TO SUBROUTINE ERRCHK.
!         ERRGET ASSIGNS TO NFATAL AND NTRACE RESPECTIVELY THE VALUES
!         OF NF AND NT IN COMMON BLOCK MLBLK0 THEREBY ASCERTAINING THE
!         STATE OF THE OPTIONS WHICH CONTROL THE EXECUTION OF ERRCHK.

!     DESCRIPTION OF ARGUMENTS
!     DESCRIPTION OF ARGUMENTS
!         BOTH ARGUMENTS ARE OUTPUT ARGUMENTS OF DATA TYPE INTEGER.
!         NFATAL - CURRENT VALUE OF NF (SEE DESCRIPTION OF ERXSET.)
!         NTRACE - CURRENT VALUE OF NT (SEE DESCRIPTION OF ERXSET.)

CALL erstgt(1,nfatal,ntrace)
RETURN
END SUBROUTINE errget

SUBROUTINE erstgt(k,nfatal,ntrace)                                00004530

!     THIS ROUTINE IS A SLAVE TO ERRGET AND ERRSET WHICH KEEPS
!     THE FLAGS AS LOCAL VARIABLES.

!     *** IF LOCAL VARIABLES ARE NOT NORMALLY RETAINED BETWEEN
!     CALLS ON THIS SYSTEM, THE VARIABLES LNF AND LNT CAN BE
!     PLACED IN A COMMON BLOCK AND PRESET TO THE FOLLOWING
!     VALUES IN THE MAIN PROGRAM.



INTEGER, INTENT(IN OUT)                  :: k
INTEGER, INTENT(IN OUT)                  :: nfatal
INTEGER, INTENT(IN OUT)                  :: ntrace
DATA lnf/-1/,lnt/0/
IF (k <= 0) lnf = nfatal
IF (k <= 0) lnt = ntrace
IF (k > 0) nfatal = lnf
IF (k > 0) ntrace = lnt
RETURN
END SUBROUTINE erstgt
