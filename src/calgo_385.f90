! Licensed under the ACM Software License Agreement
! Copyright © 1970–2012 Association for Computing Machinery (ACM)
! See ColSpecF LICENSE file for details.

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
! - 2025-06-06 - Rodrigo Castro (GitHub: rodpcastro)
!     - Replaced `DOUBLE PRECISION` by `real(wp)`.
!     - Replaced constant arrays initialized with `DATA` and stored with `SAVE` by
!       `parameter`.
!
! ## References
! 1. Kathleen A. Paciorek. 1970. Algorithm 385: Exponential integral Ei(x). Commun.
!*   ACM 13, 7 (July 1970), 446–447. <https://doi.org/10.1145/362686.362696>

  use csf_kinds, only: wp

  implicit none
  private
  public :: dei

contains

  real(wp) function dei(x)
    !! CALGO 385 Exponential integral \(\mathrm{Ei}(x)\).
    !
    !! \(\lbrace x \in \mathbb{R} \mid x \neq 0 \rbrace\)

    real(wp), intent(in) :: x  !! x ≠ 0
    real(wp) :: denm, frac
    integer :: i, j, l
    real(wp) :: px(10), qx(10)
    real(wp) :: sump, sumq
    real(wp) :: r, t
    real(wp) :: y, xmx0
    real(wp) :: maxexp
    real(wp) :: w

    real(wp), parameter :: a(6) = [ &
      -5.77215664901532863e-1_wp, 7.54164313663016620e-1_wp, &
      1.29849232927373234e-1_wp, 2.40681355683977413e-2_wp, &
      1.32084309209609371e-3_wp, 6.57739399753264501e-5_wp &
    ]
    real(wp), parameter :: b(6) = [ &
      1.0_wp, 4.25899193811589822e-1_wp, &
      7.9779471841022822e-2_wp, 8.30208476098771677e-3_wp, &
      4.86427138393016416e-4_wp, 1.30655195822848878e-5_wp &
    ]
    real(wp), parameter :: c(8) = [ &
      8.67745954838443744e-8_wp, 9.99995519301390302e-1_wp, &
      1.18483105554945844e1_wp, 4.55930644253389823e1_wp, &
      6.99279451291003023e1_wp, 4.25202034768840779e1_wp, &
      8.83671808803843939_wp, 4.01377664940664720e-1_wp &
    ]
    real(wp), parameter :: d(8) = [ &
      1.0_wp, 1.28481935379156650e1_wp, &
      5.64433569561803199e1_wp, 1.06645183769913883e2_wp, &
      8.97311097125289802e1_wp, 3.14971849170440750e1_wp, &
      3.79559003762122243_wp, 9.08804569188869219e-2_wp &
    ]
    real(wp), parameter :: e(8) = [ &
      -9.99999999999973414e-1_wp, -3.44061995006684895e1_wp, &
      -4.27532671201988539e2_wp, -2.39601943247490540e3_wp, &
      -6.16885210055476351e3_wp, -6.57609698748021179e3_wp, &
      -2.10607737142633289e3_wp, -1.48990849972948169e1_wp &
    ]
    real(wp), parameter :: f(8) = [ &
      1.0_wp, 3.64061995006459804e1_wp, &
      4.94345070209903645e2_wp, 3.19027237489543304e3_wp, &
      1.03370753085840977e4_wp, 1.63241453557783503e4_wp, &
      1.11497752871096620e4_wp, 2.37813899102160221e3_wp &
    ]
    real(wp), parameter :: p0(6) = [ &
      1.0_wp, 2.23069937666899751_wp, &
      1.70277059606809295_wp, 5.10499279623219400e-1_wp, &
      4.89089253789279154e-2_wp, 3.65462224132368429e-4_wp &
    ]
    real(wp), parameter :: p1(9) = [ &
      5.99569946892370010e9_wp, -2.50389994886351362e8_wp, &
      7.05921609590056747e8_wp, -3.36899564201591901e6_wp, &
      8.98683291643758313e6_wp, 7.37147790184657443e4_wp, &
      2.85446881813647015e4_wp, 4.12626667248911939e2_wp, &
      1.10639547241639580e1_wp &
    ]
    real(wp), parameter :: p2(9) = [ &
      9.98957666516551704e-1_wp, 5.73116705744508018_wp, &
      4.18102422562856622_wp, 5.88658240753281111_wp, &
      -1.94132967514430702e1_wp, 7.89472209294457221_wp, &
      2.32730233839039141e1_wp, -3.67783113478311458e1_wp, &
      -2.46940983448361265_wp &
    ]
    real(wp), parameter :: p3(10) = [ &
      9.99993310616056874e-1_wp, -1.84508623239127867_wp, &
      2.65257581845279982e1_wp, 2.49548773040205944e1_wp, &
      -3.32361257934396228e1_wp, -9.13483569999874255e-1_wp, &
      -2.10574079954804045e1_wp, -1.00064191398928483e1_wp, &
      -1.86009212172643758e1_wp, -1.64772117246346314_wp &
    ]
    real(wp), parameter :: p4(10) = [ &
      1.00000000000000486_wp, -3.00000000320981266_wp, &
      -5.00006640413131002_wp, -7.06810977895029359_wp, &
      -1.52856623636929637e1_wp, -7.63147701620253631_wp, &
      -2.79798528624305389e1_wp, -1.81949664929868906e1_wp, &
      -2.23127670777632410e2_wp, 1.75338801265465972e2_wp &
    ]
    real(wp), parameter :: q0(6) = [ &
      1.0_wp, 2.73069937666899751_wp, &
      2.73478695106925836_wp, 1.21765962960151532_wp, &
      2.28817933990526412e-1_wp, 1.31114151194977706e-2_wp &
    ]
    real(wp), parameter :: q1(9) = [ &
      2.55926497607616350e9_wp, -2.79673351122984591e9_wp, &
      8.02827782946956507e8_wp, -1.44980714393023883e8_wp, &
      1.77158308010799884e7_wp, -1.49575457202559218e6_wp, &
      8.53771000180749097e4_wp, -3.02523682238227410e3_wp, &
      5.12578125e1_wp &
    ]
    real(wp), parameter :: q2(8) = [ &
      1.14625253249016191_wp, -1.99149600231235164e2_wp, &
      3.41365212524375539e2_wp, 5.23165568734558614e1_wp, &
      3.17279489254369328e2_wp, -8.38767084189640707_wp, &
      9.65405217429280303e2_wp, 2.63983007318024593_wp &
    ]
    real(wp), parameter :: q3(9) = [ &
      1.00153385204534270_wp, -1.09355619539109124e1_wp, &
      1.99100447081774247e2_wp, 1.19283242396860101e3_wp, &
      4.42941317833792840e1_wp, 2.53881931563070803e2_wp, &
      5.99493232566740736e1_wp, 6.40380040535241555e1_wp, &
      9.79240359921729030e1_wp &
    ]
    real(wp), parameter :: q4(9) = [ &
      1.99999999999048104_wp, -2.99999894040324960_wp, &
      -7.99243595776339741_wp, -1.20187763547154743e1_wp, &
      7.04831847180424676e1_wp, 1.17179220502086455e2_wp, &
      1.37790390235747999e2_wp, 3.97277109100414518_wp, &
      3.97845977167414721e4_wp &
    ]
    real(wp), parameter :: x0 = 0.372507410781366634_wp

    ! MAXEXP needs to be set to the largest argument of exp
    ! that will not cause an overflow. This is computed here
    ! but could be embedded as a constant for efficiency reasons.
    maxexp = (INT(LOG(huge(0.0_wp)) * 100.0_wp)) / 100.0_wp

  1 IF (x <= 0.0_wp) GO TO 100
    IF (x >= 12.0_wp) GO TO 60
    IF (x >= 6.0_wp) GO TO 40

    ! X IN (0,6).

    t = x + x
    t = t / 3.0_wp - 2.0_wp
    px(10) = 0.0_wp
    qx(10) = 0.0_wp
    px(9) = p1(9)
    qx(9) = q1(9)

    ! THE RATIONAL FUNCTION IS EXPRESSED AS A RATIO OF FINITE SUMS OF
    ! SHIFTED CHEBYSHEV POLYNOMIALS, AND IS EVALUATED BY NOTING THAT
    ! T*(X) = T(2*X-1) AND USING THE CLENSHAW-RICE ALGORITHM FOUND IN
    ! REFERENCE (4).

    DO l = 2, 8
      i = 10 - l
      px(i) = t * px(i+1) - px(i+2) + p1(i)
      qx(i) = t * qx(i+1) - qx(i+2) + q1(i)
    END DO

    r = (0.5_wp * t * px(2) - px(3) + p1(1)) / (0.5_wp * t * qx(2) - qx(3) + q1(1))

    ! ( X - X0 ) = ( X - X1 ) - X2, WHERE X1 = 409576229586. / 2**40 AND
    ! X2 = -.7671772501993940e-12_wp.

    xmx0 = (x - 409576229586.0_wp / 1099511627776.0_wp) - 0.7671772501993940e-12_wp
    IF (ABS(xmx0) < 0.037_wp) GO TO 15
    dei = LOG(x / x0) + xmx0 * r
    RETURN

  15 y = xmx0 / x0

    ! A RATIONAL APPROXIMATION TO LOG ( X / X0 ) * LOG ( 1 + Y ),
    ! WHERE Y = ( X - X0 ) / X0, AND DABS ( Y ) IS LESS THAN 0.1,
    ! THAT IS FOR DABS ( X - X0 ) LESS THAN 0.037.

    sump = ((((p0(6) * y + p0(5)) * y + p0(4)) * y + p0(3)) * y + p0(2)) * y + p0(1)
    sumq = ((((q0(6) * y + q0(5)) * y + q0(4)) * y + q0(3)) * y + q0(2)) * y + q0(1)

    dei = (sump / (sumq * x0) + r) * xmx0
    RETURN

    ! X IN (6,12).

  40 denm = p2(9) + x
    frac = q2(8) / denm

    ! THE RATIONAL FUNCTION IS EXPRESSED AS A J-FRACTION.

    DO j = 2, 8
      i = 9 - j
      denm = p2(i+1) + x + frac
      frac = q2(i) / denm
    END DO

    dei = EXP(x) * (p2(1) + frac) / x
    RETURN

  60 IF (x >= 24.0_wp) GO TO 80

    ! X IN (12,24).

    denm = p3(10) + x
    frac = q3(9) / denm

    ! THE RATIONAL FUNCTION IS EXPRESSED AS A J-FRACTION.

    DO j = 2, 9
      i = 10 - j
      denm = p3(i+1) + x + frac
      frac = q3(i) / denm
    END DO

    dei = exp(x) * (p3(1) + frac) / x
    RETURN

    ! X GREATER THAN 24.

  80 IF (x <= maxexp) GO TO 90

    ! X IS GREATER THAN MAXEXP AND DEI IS SET TO INFINITY.

    dei = huge(0.0_wp)
    RETURN
  90 y = 1.0_wp / x
    denm = p4(10) + x
    frac = q4(9) / denm

    ! THE RATIONAL FUNCTION IS EXPRESSED AS A J-FRACTION.

    DO j = 2, 9
      i = 10 - j
      denm = p4(i+1) + x + frac
      frac = q4(i) / denm
    END DO

    dei = exp(x) * (y + y * y * (p4(1) + frac))
    RETURN

  100 IF (x /= 0.0_wp) GO TO 101

    ! X = 0 AND DEI IS SET TO -INFINITY.

    dei = -huge(0.0_wp)
    WRITE(*, 500)
  500 FORMAT('DEI CALLED WITH A ZERO ARGUMENT, RESULT SET TO -INFINITY')
    RETURN

  101 y = -x
  110 w = 1.0_wp / y
    IF (y > 4.0_wp) GO TO 300
    IF (y > 1.0_wp) GO TO 200

    ! X IN (-1, 0).

    dei = LOG(y) - ( &
      ((((a(6) * y + a(5)) * y + a(4)) * y + a(4)) * y + a(2)) * y + a(1) &
    ) / ( &
      ((((b(6) * y + b(5)) * y + b(4)) * y + b(4)) * y + b(2)) * y + b(1) &
    )
    RETURN

    ! X IN (-4, -1).

  200 dei = -exp(-y) * ( &
        ( &
          ((((((c(8) * w + c(7)) * w + c(6)) * w + c(5)) &
          * w + c(4)) * w + c(3)) * w + c(2)) * w + c(1) &
        ) / ( &
          ((((((d(8) * w + d(7)) * w + d(6)) * w + d(5)) &
          * w + d(4)) * w + d(3)) * w + d(2)) * w + d(1) &
        ) &
      )
    RETURN

    ! X LESS THAN -4.

  300 dei = -exp(-y) * ( &
        w * ( &
          1.0_wp + w * ( &
            ((((((e(8) * w + e(7)) * w + e(6)) * w + e(5)) &
            * w + e(4)) * w + e(3)) * w + e(2)) * w + e(1) &
          ) / ( &
            ((((((f(8) * w + f(7)) * w + f(6)) * w + f(5)) &
            * w + f(4)) * w + f(3)) * w + f(2)) * w + f(1) &
          ) &
        ) &
      )
    RETURN
  end function dei

end module calgo_385
