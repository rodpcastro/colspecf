PROGRAM cqccex

! Code converted using TO_F90 by Alan Miller
! Date: 2002-03-26  Time: 19:32:11

!  CQCCEX IS A QUICK CHECK PROGRAM TO COMPARE EXPONENTIAL INTEGRALS
!  E(N,Z) FROM SUBROUTINE CEXINT,

!            CALL CEXINT(Z,N,KODE,TOL,M,CY,IERR)

!  AGAINST EXPONENTIAL INTEGRALS FROM QUADRATURE SUBROUTINE CEXQAD,

!            CALL CEXQAD(Z,N,KODE,TOL,CQ,KERR).

!  Z VALUES ARE TAKEN FROM THE REGION -6.5.LE.X.LE.5.5,-6.LE.Y.LE.6.
!  ORDERS N RUN FROM 3 TO 11 AND THE NUMBER OF MEMBERS M IN THE
!  SEQUENCE E(N+K-1,Z), K=1,M RUNS FROM 1 TO 3. BOTH SCALING OPTIONS

!             CY(K) = E(N+K-1,Z)                K=1,M     KODE=1
!                     E(N+K-1,Z)*CEXP(Z)        K=1,M     KODE=2

!  ARE CHECKED AND THE REQUESTED ACCURACY TOL IS THE LARGER OF
!  UNIT ROUNDOFF AND 1.0E-7. RELATIVE ERRORS ERR1 AND ERR2 FOR THE
!  FIRST AND LAST MEMBERS OF THE SEQUENCE ARE COMPUTED AND COMPARED
!  AGAINST 100.0*TOL. IF A CHECK DOES NOT OCCUR, Z,ERR1,ERR2,N,M,KODE
!  AND ERROR FLAGS IERR FROM CEXINT, KERR1 FROM CEXQAD, AND KERR2
!  FROM CEXQAD ARE PRINTED. VALUES CY(1),CQ1 AND CY(N+M-1),CQ2 WHICH
!  WERE COMPARED IN ERR1 AND ERR2 ARE PRINTED NEXT. KERR1.NE.0 OR
!  KERR2.NE.0 INDICATE A PREMATURE TRUNCATION OF THE INTEGRAL EVAL-
!  UATION IN CEXQAD. THE SUFFIXES 1 AND 2 CORRESPOND TO EVALUATIONS
!  AT ORDERS N AND N+M-1 RESPECTIVELY.

!  CQCCEX CALLS CEXINT,CEXQAD AND LOWER LEVEL ROUTINES CEXENZ,CACEXI,
!  PSIXN,CEXQAD,GAUS8,FQCEX,I1MACH,R1MACH,XERROR,FDUMP

USE calgo_683
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

COMPLEX (dp)  :: z, cy(10), cq1, cq2
INTEGER       :: ierr, iprnt, ix, iy, kerr1, kerr2, kode, m, n, nb
REAL (dp)     :: err1, err2, tol, x, xtol, y

iprnt = 0
tol = MAX(1.0D-15, EPSILON(0.0_dp))
xtol = tol * 100.0_dp
WRITE(*, 5000)
WRITE(*, 5100)
DO  kode = 1, 2
  DO  m = 1, 3
    DO  n = 3, 11, 2
      DO  iy = 1, 13, 3
        y = iy - 7
        WRITE(*, 5200) kode, m, n, y
        DO  ix = 2, 14, 4
          x = -7.5_dp + (ix-1)
          IF (y /= 0.0_dp .OR. x > 0.0_dp) THEN
            z = CMPLX(x, y, KIND=dp)
            CALL cexint(z, n, kode, tol, m, cy, ierr)
            nb = n
            CALL cexqad(z, nb, kode, tol, cq1, kerr1)
            err1 = ABS(cy(1)-cq1) / ABS(cq1)
            nb = n + m - 1
            CALL cexqad(z, nb, kode, tol, cq2, kerr2)
            err2 = ABS(cy(m)-cq2) / ABS(cq2)
            IF (err1 >= xtol .OR. err2 >= xtol) THEN
              IF (iprnt /= 1) THEN
                iprnt = 1
                OPEN (7, FILE='CEXDAT7', STATUS='UNKNOWN')
                WRITE (7,5300)
              END IF
              WRITE (7,5400) z, err1, err2, n, m, kode, ierr, kerr1, kerr2
              WRITE (7,5400) cy(1), cq1
              WRITE (7,5400) cy(m), cq2
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO
END DO
IF (iprnt == 0) THEN
  WRITE(*, 5500)
ELSE
  WRITE(*, 5600)
END IF
STOP

5000 FORMAT (' CEXINT VS QUADRATURE FOR PARAMETERS:'/)
5100 FORMAT ('  KODE   M     N       Y        -6.5.LE.X.LE.5.5')
5200 FORMAT (i5, i5, i5, g13.5)
5300 FORMAT (' Z,ERR1,ERR2,N,M,KODE,IERR,KERR1,KERR2')
5400 FORMAT (4g11.3, 6I5)
5500 FORMAT (/' QUICK CHECKS FOR CEXINT ARE OK.')
5600 FORMAT (/' SEE DATA FILE CEXDAT7 FOR ERRORS')
END PROGRAM cqccex
