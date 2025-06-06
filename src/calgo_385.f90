! Licensed under the ACM Software License Agreement
! Copyright © 1970–2012 Association for Computing Machinery (ACM).
! https://www.acm.org/publications/policies/software-copyright-notice

module calgo_385
!* # CALGO 385
! Algorithm 385.
!
! Procedures:
!
! - `dei`: Exponential integral \(\mathrm{Ei}(x)\)
!
! ## Author
! Kathleen Paciorek
!
! ## History
! - 1970-07-01 - Kathleen Paciorek
!     - Original code.
! - 2006-10-04 - Kathleen Paciorek
!     - F77 code distributed by ACM: <https://calgo.acm.org/385.zip>
! - 2025-05-26 - Rodrigo Castro (GitHub: rodpcastro)
!     - Converted the code from F77 to F90 using TO_F90 by Alan Miller.
! - 2025-05-26 - Rodrigo Castro (GitHub: rodpcastro)
!     - Fixed `frac = q2(8) + x` to `frac = q2(8) / denm` for x in the interval
!       [6,12] according to the algorithm in the reference publication.
!
! ## References
! 1. Kathleen A. Paciorek. 1970. Algorithm 385: Exponential integral Ei(x). Commun.
!*   ACM 13, 7 (July 1970), 446–447. <https://doi.org/10.1145/362686.362696>

! TODO: Refactor the code to modern Fortran.

  implicit none
  private
  public :: dei

contains

  function dei(x1)
    !! CALGO 385 Exponential integral \(\mathrm{Ei}(x)\).
    !
    !! \(\lbrace x \in \mathbb{R} \mid x \neq 0 \rbrace\)

    DOUBLE PRECISION, INTENT(IN) :: x1  !! x1 ≠ 0
    DOUBLE PRECISION :: a(6)
    DOUBLE PRECISION :: b(6)
    DOUBLE PRECISION :: c(8)
    DOUBLE PRECISION :: d(8)
    DOUBLE PRECISION :: dei
    DOUBLE PRECISION :: denm
    DOUBLE PRECISION :: e(8)
    DOUBLE PRECISION :: f(8)
    DOUBLE PRECISION :: frac
    INTEGER :: i
    INTEGER :: j
    INTEGER :: l
    DOUBLE PRECISION :: p0(6)
    DOUBLE PRECISION :: p1(9)
    DOUBLE PRECISION :: p2(9)
    DOUBLE PRECISION :: p3(10)
    DOUBLE PRECISION :: p4(10)
    DOUBLE PRECISION :: px(10)
    DOUBLE PRECISION :: q0(6)
    DOUBLE PRECISION :: q1(9)
    DOUBLE PRECISION :: q2(8)
    DOUBLE PRECISION :: q3(9)
    DOUBLE PRECISION :: q4(9)
    DOUBLE PRECISION :: qx(10)
    DOUBLE PRECISION :: r
    DOUBLE PRECISION :: sump
    DOUBLE PRECISION :: sumq
    DOUBLE PRECISION :: t
    DOUBLE PRECISION :: w
    DOUBLE PRECISION :: x
    DOUBLE PRECISION :: x0

    DOUBLE PRECISION :: xmx0
    DOUBLE PRECISION :: y, maxexp

    SAVE a
    SAVE b
    SAVE c
    SAVE d
    SAVE e
    SAVE f
    SAVE p0
    SAVE p1
    SAVE p2
    SAVE p3
    SAVE p4
    SAVE q0
    SAVE q1
    SAVE q2
    SAVE q3
    SAVE q4
    SAVE x0

    DATA a / -5.77215664901532863D-01,  &
        7.54164313663016620D-01, 1.29849232927373234D-01,  &
        2.40681355683977413D-02, 1.32084309209609371D-03,  &
        6.57739399753264501D-05 /
    DATA b / 1.0D+00,  &
        4.25899193811589822D-01, 7.9779471841022822D-02,  &
        8.30208476098771677D-03, 4.86427138393016416D-04,  &
        1.30655195822848878D-05 /
    DATA c / 8.67745954838443744D-08,  &
        9.99995519301390302D-01, 1.18483105554945844D+01,  &
        4.55930644253389823D+01, 6.99279451291003023D+01,  &
        4.25202034768840779D+01, 8.83671808803843939D+00,  &
        4.01377664940664720D-01 /
    DATA d / 1.0D+00,  &
        1.28481935379156650D+01, 5.64433569561803199D+01,  &
        1.06645183769913883D+02, 8.97311097125289802D+01,  &
        3.14971849170440750D+01, 3.79559003762122243D+00,  &
        9.08804569188869219D-02 /
    DATA e / -9.99999999999973414D-01,  &
        -3.44061995006684895D+01, -4.27532671201988539D+02,  &
        -2.39601943247490540D+03, -6.16885210055476351D+03,  &
        -6.57609698748021179D+03, -2.10607737142633289D+03,  &
        -1.48990849972948169D+01 /
    DATA f / 1.0D+00,  &
        3.64061995006459804D+01, 4.94345070209903645D+02,  &
        3.19027237489543304D+03, 1.03370753085840977D+04,  &
        1.63241453557783503D+04, 1.11497752871096620D+04,  &
        2.37813899102160221D+03 /
    DATA p0 / 1.0D+00,  &
        2.23069937666899751D+00, 1.70277059606809295D+00,  &
        5.10499279623219400D-01, 4.89089253789279154D-02,  &
        3.65462224132368429D-04 /
    DATA p1 / 5.99569946892370010D+09,  &
        -2.50389994886351362D+08, 7.05921609590056747D+08,  &
        -3.36899564201591901D+06, 8.98683291643758313D+06,  &
        7.37147790184657443D+04, 2.85446881813647015D+04,  &
        4.12626667248911939D+02, 1.10639547241639580D+01 /
    DATA p2 / 9.98957666516551704D-01,  &
        5.73116705744508018D+00, 4.18102422562856622D+00,  &
        5.88658240753281111D+00, -1.94132967514430702D+01,  &
        7.89472209294457221D+00, 2.32730233839039141D+01,  &
        -3.67783113478311458D+01, -2.46940983448361265D+00 /
    DATA p3 / 9.99993310616056874D-01,  &
        -1.84508623239127867D+00, 2.65257581845279982D+01,  &
        2.49548773040205944D+01, -3.32361257934396228D+01,  &
        -9.13483569999874255D-01, -2.10574079954804045D+01,  &
        -1.00064191398928483D+01, -1.86009212172643758D+01,  &
        -1.64772117246346314D+00 /
    DATA p4 / 1.00000000000000486D+00,  &
        -3.00000000320981266D+00, -5.00006640413131002D+00,  &
        -7.06810977895029359D+00, -1.52856623636929637D+01,  &
        -7.63147701620253631D+00, -2.79798528624305389D+01,  &
        -1.81949664929868906D+01, -2.23127670777632410D+02,  &
        1.75338801265465972D+02 /
    DATA q0 / 1.0D+00,  &
        2.73069937666899751D+00, 2.73478695106925836D+00,  &
        1.21765962960151532D+00, 2.28817933990526412D-01,  &
        1.31114151194977706D-02 /
    DATA q1 / 2.55926497607616350D+09,  &
        -2.79673351122984591D+09, 8.02827782946956507D+08,  &
        -1.44980714393023883D+08, 1.77158308010799884D+07,  &
        -1.49575457202559218D+06, 8.53771000180749097D+04,  &
        -3.02523682238227410D+03, 5.12578125D+01 /
    DATA q2 / 1.14625253249016191D+00,  &
        -1.99149600231235164D+02, 3.41365212524375539D+02,  &
        5.23165568734558614D+01, 3.17279489254369328D+02,  &
        -8.38767084189640707D+00, 9.65405217429280303D+02,  &
        2.63983007318024593D+00 /
    DATA q3 / 1.00153385204534270D+00,  &
        -1.09355619539109124D+01, 1.99100447081774247D+02,  &
        1.19283242396860101D+03, 4.42941317833792840D+01,  &
        2.53881931563070803D+02, 5.99493232566740736D+01,  &
        6.40380040535241555D+01, 9.79240359921729030D+01 /
    DATA q4 / 1.99999999999048104D+00,  &
        -2.99999894040324960D+00, -7.99243595776339741D+00,  &
        -1.20187763547154743D+01, 7.04831847180424676D+01,  &
        1.17179220502086455D+02, 1.37790390235747999D+02,  &
        3.97277109100414518D+00, 3.97845977167414721D+04 /
    DATA x0 / 0.372507410781366634D+00 /

    ! MAXEXP needs to be set to the largest argument of exp
    ! that will not cause an overflow. This is computed here
    ! but could be embedded as a constant for efficiency reasons.
    maxexp = (INT(LOG(huge(0.0D0))*100))/100.0D0

    x = x1
    1     IF ( x <= 0.0D+00 ) GO TO 100
    IF ( x >= 12.0D+00 ) GO TO 60
    IF ( x >= 6.0D+00 ) GO TO 40

    !  X IN (0,6).

    t = x + x
    t = t / 3.0D+00 - 2.0D+00
    px(10) = 0.0D+00
    qx(10) = 0.0D+00
    px(9) = p1(9)
    qx(9) = q1(9)

    !  THE RATIONAL FUNCTION IS EXPRESSED AS A RATIO OF FINITE SUMS OF
    !  SHIFTED CHEBYSHEV POLYNOMIALS, AND IS EVALUATED BY NOTING THAT
    !  T*(X) = T(2*X-1) AND USING THE CLENSHAW-RICE ALGORITHM FOUND IN
    !  REFERENCE (4).

    DO  l = 2, 8
      i = 10 - l
      px(i) = t * px(i+1) - px(i+2) + p1(i)
      qx(i) = t * qx(i+1) - qx(i+2) + q1(i)
    END DO

    r = ( 0.5D+00 * t * px(2) - px(3) + p1(1) )  &
        / ( 0.5D+00 * t * qx(2) - qx(3) + q1(1) )

    !  ( X - X0 ) = ( X - X1 ) - X2, WHERE X1 = 409576229586. / 2**40 AND
    !  X2 = -.7671772501993940D-12.

    xmx0 = ( x - 409576229586.0D+00 / 1099511627776.0D+00 )  &
        - 0.7671772501993940D-12
    IF ( DABS ( xmx0 ) < 0.037D+00 ) GO TO 15
    dei = LOG ( x / x0 ) + xmx0 * r
    RETURN
    15    y = xmx0 / x0

    !  A RATIONAL APPROXIMATION TO LOG ( X / X0 ) * LOG ( 1 + Y ),
    !  WHERE Y = ( X - X0 ) / X0, AND DABS ( Y ) IS LESS THAN 0.1,
    !  THAT IS FOR DABS ( X - X0 ) LESS THAN 0.037.

    sump = (((( p0(6) * y + p0(5) )  &
        * y + p0(4) ) * y + p0(3) )  &
        * y + p0(2) ) * y + p0(1)

    sumq = (((( q0(6) * y + q0(5) )  &
        * y + q0(4) ) * y + q0(3) )  &
        * y + q0(2) ) * y + q0(1)

    dei = ( sump / ( sumq * x0 ) + r ) * xmx0
    RETURN

    !  X IN (6,12).

    40    denm = p2(9) + x
    frac = q2(8) / denm

    !  THE RATIONAL FUNCTION IS EXPRESSED AS A J-FRACTION.

    DO  j = 2, 8
      i = 9 - j
      denm = p2(i+1) + x + frac
      frac = q2(i) / denm
    END DO

    dei = EXP ( x ) * ( ( p2(1) + frac ) / x )
    RETURN

    60    IF ( x >= 24.0D+00 ) GO TO 80

    !  X IN (12,24).

    denm = p3(10) + x
    frac = q3(9) / denm

    !  THE RATIONAL FUNCTION IS EXPRESSED AS A J-FRACTION.

    DO  j = 2, 9
      i = 10 - j
      denm = p3(i+1) + x + frac
      frac = q3(i) / denm
    END DO

    dei = EXP ( x ) * ( ( p3(1) + frac ) / x )
    RETURN

    !  X GREATER THAN 24.

    80    IF ( x <= maxexp ) GO TO 90

    !  X IS GREATER THAN MAXEXP AND DEI IS SET TO INFINITY.

    dei = huge(0.0D0)
    RETURN
    90    y = 1.0D+00 / x
    denm = p4(10) + x
    frac = q4(9) / denm

    !  THE RATIONAL FUNCTION IS EXPRESSED AS A J-FRACTION.

    DO  j = 2, 9
      i = 10 - j
      denm = p4(i+1) + x + frac
      frac = q4(i) / denm
    END DO

    dei = EXP ( x ) * ( y + y * y * ( p4(1) + frac ) )
    RETURN

    100   IF ( x /= 0.0D+00 ) GO TO 101

    !  X = 0 AND DEI IS SET TO -INFINITY.

    dei = -huge(0.0D0)
    WRITE(*,500)
    500   FORMAT ( ' DEI CALLED WITH A ZERO ARGUMENT, RESULT SET TO -INFINITY')
    RETURN
    101   y = -x
    110   w = 1.0D+00 / y
    IF ( y > 4.0D+00 ) GO TO 300
    IF ( y > 1.0D+00 ) GO TO 200

    !  X IN (-1,0).

    dei = LOG ( y ) - ((((( a(6)  &
        * y + a(5) ) * y + a(4) )  &
        * y + a(4) ) * y + a(2) )  &
        * y + a(1) ) / ((((( b(6)  &
        * y + b(5) ) * y + b(4) )  &
        * y + b(4) ) * y + b(2) )  &
        * y + b(1) )
    RETURN

    !  X IN (-4,-1).

    200   dei = -EXP ( -y ) * (((((((( c(8)  &
        * w + c(7) ) * w + c(6) )  &
        * w + c(5) ) * w + c(4) )  &
        * w + c(3) ) * w + c(2) )  &
        * w + c(1) ) / ((((((( d(8)  &
        * w + d(7) ) * w + d(6) )  &
        * w + d(5) ) * w + d(4) )  &
        * w + d(3) ) * w + d(2) )  &
        * w + d(1) ) )
    RETURN

    !  X LESS THAN -4.

    300   dei = -EXP ( -y ) * ( w * ( 1.0D+00 + w * ((((((( e(8)  &
        * w + e(7) ) * w + e(6) )  &
        * w + e(5) ) * w + e(4) )  &
        * w + e(3) ) * w + e(2) )  &
        * w + e(1) ) / ((((((( f(8)  &
        * w + f(7) ) * w + f(6) )  &
        * w + f(5) ) * w + f(4) )  &
        * w + f(3) ) * w + f(2) )  &
        * w + f(1) ) ) )

    RETURN
  end function dei

end module calgo_385
