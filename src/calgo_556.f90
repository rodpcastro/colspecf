module calgo_556
  
  implicit none

contains

! Code converted using TO_F90 by Alan Miller
! Date: 2025-05-26  Time: 03:12:12

!--**--CH775--556--Fix--2:8:1999
!--**--CH774--556--A:1--2:8:1999
!--**--CH566--556--A:H--29:7:1999
!--**--CH565--556--U:D--29:7:1999

SUBROUTINE expint(x,n,kode,m,tol,en,ierr)

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

!     .. Scalar Arguments ..

REAL, INTENT(IN)                         :: x
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: kode
INTEGER, INTENT(IN)                      :: m
REAL, INTENT(IN)                         :: tol
REAL, INTENT(OUT)                        :: en(m)
INTEGER, INTENT(OUT)                     :: ierr


!     ..
!     .. Array Arguments ..

!     ..
!     .. Local Scalars ..
REAL :: aa,aams,ah,ak,at,bk,bt,cc,cnorm,ct,em,emx,etol,euler,fnm,fx,  &
    p1,p2,pt,s,tx,xcut,xlim,xtol,y1,y2,yt
INTEGER :: i,ic,icase,ict,ik,ind,ix,jset,k,kk,kn,ks,lun,ml,mu,nd,nm
!     ..
!     .. Local Arrays ..
REAL :: a(99),b(99),y(2)
!     ..
!     .. External Functions ..
REAL :: digam
INTEGER :: i1mach
EXTERNAL digam,i1mach
!     ..
!     .. Intrinsic Functions ..
INTRINSIC ABS,ALOG,EXP,FLOAT,INT,LOG,REAL,SQRT
!     ..
!     .. Data statements ..

!  For epcf90 on SUN Ultra 10  set
!     XLIM = 87 (sp) and 708 (dp) and XCUT = 1.0

DATA xcut,xlim/1.0E0,87.0E0/
DATA euler/-5.77215664901533E-01/
!     ..
lun = i1mach(2)

etol = 10.0** (-INT((i1mach(11)*LOG(REAL(i1mach(10))))/LOG(10.0)))
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
DO  i = 1,35
  aa = -aa*x/ak
  IF (i == nm) GO TO 30
  s = s - aa/ (ak-fnm)
  IF (ABS(aa) <= xtol*ABS(s)) GO TO 20
  ak = ak + 1.0E0
  CYCLE
  
  20     CONTINUE
  IF (i < 2) GO TO 40
  IF (nd-2 > i .OR. i > nd-1) GO TO 60
  ak = ak + 1.0E0
  CYCLE
  
  30     s = s + aa* (-ALOG(x)+digam(nd))
  xtol = 3.0E0*tol
  40     ak = ak + 1.0E0
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
    GO TO 240
END SELECT

70 IF (icase == 2) RETURN
IF (kode == 1) emx = EXP(-x)
en(1) = (emx-s)/x
RETURN

80 CONTINUE
DO  i = 1,m
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
DO  i = 1,m
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
em = (ah+1.0E0)/ ((x+aa)*xtol*SQRT(ct))
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
at = bk/ (bk+ak+cc+FLOAT(ic))
bk = bk + ak + ak
a(ic) = at
bt = (ak+ak+x)/ (ak+1.0E0)
b(ic) = bt
pt = p2
p2 = bt*p2 - at*p1
p1 = pt
ct = ct + fx
em = em*at* (1.0E0-tx/ct)
IF (em* (ak+1.0E0) > p1*p1) GO TO 190
ict = ic
kk = ic + 1
bt = tx/ (ct+fx)
y2 = (bk/ (bk+cc+FLOAT(kk)))* (p1/p2)* (1.0E0-bt+0.375E0*bt*bt)
y1 = 1.0E0
!     BACKWARD RECURRENCE FOR
!              Y1=             C*U( A ,A,X)
!              Y2= C*(A/(1+A/2))*U(A+1,A,X)
DO  k = 1,ict
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
pt = y2/y1
cnorm = 1.0E0 - pt* (ah+1.0E0)/aa
y(1) = 1.0E0/ (cnorm*aa+x)
y(2) = cnorm*y(1)
IF (icase == 3) GO TO 210
en(ind) = emx*y(jset)
IF (m == 1) RETURN
aa = FLOAT(ks)
SELECT CASE ( icase )
  CASE (    1)
    GO TO 220
  CASE (    2)
    GO TO 240
END SELECT

!     RECURSION SECTION  N*E(N+1,X) + X*E(N,X)=EMX

210 en(1) = emx* (1.0E0-y(1))/x
RETURN

220 k = ind - 1
DO  i = 1,ml
  aa = aa - 1.0E0
  en(k) = (emx-aa*en(k+1))/x
  k = k - 1
END DO
IF (mu <= 0) RETURN
aa = FLOAT(ks)
240 k = ind
DO  i = 1,mu
  en(k+1) = (emx-x*en(k))/aa
  aa = aa + 1.0E0
  k = k + 1
END DO
RETURN


260 WRITE (lun,FMT=9000)
RETURN

270 WRITE (lun,FMT=9010)
RETURN

280 WRITE (lun,FMT=9020)
RETURN

290 WRITE (lun,FMT=9030)
RETURN

300 WRITE (lun,FMT=9040)
RETURN

310 WRITE (lun,FMT=9050)
RETURN

320 WRITE (lun,FMT=9060)
RETURN

330 WRITE (lun,FMT=9070)
RETURN

340 WRITE (lun,FMT=9080)
RETURN

9000 FORMAT (' IN EXPINT, N NOT GREATER THAN 0')
9010 FORMAT (' IN EXPINT, KODE NOT 1 OR 2')
9020 FORMAT (' IN EXPINT, M NOT GREATER THAN 0')
9030 FORMAT (' IN EXPINT, TOL NOT WITHIN LIMITS')
9040 FORMAT (' IN EXPINT, X IS NOT ZERO OR POSITIVE')
9050 FORMAT (' IN EXPINT, THE EXPONENTIAL INTEGRAL IS NOT DE','FINED ',  &
    'FOR X=0 AND N=1')
9060 FORMAT (' IN EXPINT, RELATIVE ERROR TEST FOR SERIES TER','MINATI',  &
    'ON NOT MET IN 36 TERMS')
9070 FORMAT (' IN EXPINT, TERMINATION TEST FOR MILLER ALGORI','THM NO',  &
    'T MET IN 99 STEPS')
9080 FORMAT (' IN EXPINT, AN ERROR IN PLACING INT(X+0.5) WIT','H RESP',  &
    'ECT TO N AND N+M-1 OCCURRED FOR X.GT.XCUT')
END SUBROUTINE expint

REAL FUNCTION digam(n)

!     THIS SUBROUTINE RETURNS VALUES OF PSI(X)=DERIVATIVE OF LOG
!     GAMMA(X), X.GT.0.0 AT INTEGER ARGUMENTS. A TABLE LOOK-UP IS
!     PERFORMED FOR N.LE.100, AND THE ASYMPTOTIC EXPANSION IS
!     EVALUATED FOR N.GT.100.

!     .. Scalar Arguments ..

INTEGER, INTENT(IN)                      :: n

!     ..
!     .. Local Scalars ..
REAL :: ak,ax,fn,fn2,s
INTEGER :: k
!     ..
!     .. Local Arrays ..
REAL :: b(4),c(100),c1(32),c2(27),c3(22),c4(19)
!     ..
!     .. Intrinsic Functions ..
INTRINSIC ALOG
!     ..
!     .. Equivalences ..
EQUIVALENCE (c(1),c1(1))
EQUIVALENCE (c(33),c2(1))
EQUIVALENCE (c(60),c3(1))
EQUIVALENCE (c(82),c4(1))
!     ..
!     .. Data statements ..


DATA c1/-5.7721566490153E-01,4.22784335098467E-01,  &
    9.22784335098467E-01,1.25611766843180E+00,  &
    1.50611766843180E+00,1.70611766843180E+00,  &
    1.87278433509847E+00,2.01564147795561E+00,  &
    2.14064147795561E+00,2.25175258906672E+00,  &
    2.35175258906672E+00,2.44266167997581E+00,  &
    2.52599501330915E+00,2.60291809023222E+00,  &
    2.67434666166079E+00,2.74101332832746E+00,  &
    2.80351332832746E+00,2.86233685773923E+00,  &
    2.91789241329478E+00,2.97052399224215E+00,  &
    3.02052399224215E+00,3.06814303986120E+00,  &
    3.11359758531574E+00,3.15707584618531E+00,  &
    3.19874251285197E+00,3.23874251285197E+00,  &
    3.27720405131351E+00,3.31424108835055E+00,  &
    3.34995537406484E+00,3.38443813268552E+00,  &
    3.41777146601886E+00,3.45002953053499E+00/
DATA c2/3.48127953053499E+00,3.51158256083802E+00,  &
    3.54099432554390E+00,3.56956575411533E+00,  &
    3.59734353189311E+00,3.62437055892013E+00,  &
    3.65068634839382E+00,3.67632737403484E+00,  &
    3.70132737403484E+00,3.72571761793728E+00,  &
    3.74952714174681E+00,3.77278295570029E+00,  &
    3.79551022842757E+00,3.81773245064979E+00,  &
    3.83947158108457E+00,3.86074817682925E+00,  &
    3.88158151016259E+00,3.90198967342789E+00,  &
    3.92198967342789E+00,3.94159751656515E+00,  &
    3.96082828579592E+00,3.97969621032422E+00,  &
    3.99821472884274E+00,4.01639654702455E+00,  &
    4.03425368988170E+00,4.05179754953082E+00, 4.06903892884117E+00/
DATA c3/4.08598808138354E+00,4.10265474805020E+00,  &
    4.11904819067316E+00,4.13517722293122E+00,  &
    4.15105023880424E+00,4.16667523880424E+00,  &
    4.18205985418885E+00,4.19721136934037E+00,  &
    4.21213674247470E+00,4.22684262482764E+00,  &
    4.24133537845082E+00,4.25562109273654E+00,  &
    4.26970559977879E+00,4.28359448866768E+00,  &
    4.29729311880467E+00,4.31080663231818E+00,  &
    4.32413996565151E+00,4.33729786038836E+00,  &
    4.35028487337537E+00,4.36310538619588E+00,  &
    4.37576361404398E+00,4.38826361404398E+00/
DATA c4/4.40060929305633E+00,4.41280441500755E+00,  &
    4.42485260777863E+00,4.43675736968340E+00,  &
    4.44852207556575E+00,4.46014998254249E+00,  &
    4.47164423541606E+00,4.48300787177969E+00,  &
    4.49424382683587E+00,4.50535493794698E+00,  &
    4.51634394893599E+00,4.52721351415338E+00,  &
    4.53796620232543E+00,4.54860450019777E+00,  &
    4.55913081598724E+00,4.56954748265391E+00,  &
    4.57985676100442E+00,4.59006084263708E+00, 4.60016185273809E+00/
DATA b/1.66666666666667E-01,-3.33333333333333E-02,  &
    2.38095238095238E-02,-3.33333333333333E-02/
!     ..

IF (n > 100) GO TO 10
digam = c(n)
RETURN

10 fn = n
ax = 1.0E0
ak = 2.0E0
s = -0.5E0/fn
IF (fn > 1.e+8) GO TO 30
fn2 = fn*fn
DO  k = 1,3
  ax = ax*fn2
  s = s - b(k)/ (ax*ak)
  ak = ak + 2.0E0
END DO
30 CONTINUE
digam = s + ALOG(fn)
RETURN

END FUNCTION digam

end module calgo_556
