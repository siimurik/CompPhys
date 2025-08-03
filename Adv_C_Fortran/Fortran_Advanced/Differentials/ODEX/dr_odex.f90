PROGRAM VanDerPolDriver
    IMPLICIT NONE
    INTEGER, PARAMETER :: NDGL = 2, KM = 9, NRDENS = 2
    INTEGER, PARAMETER :: LWORK = NDGL*(KM+5) + 5*KM + 20 + (2*KM*(KM+2)+5)*NRDENS
    INTEGER, PARAMETER :: LIWORK = 2*KM + 21 + NRDENS

    INTEGER :: N, IOUT, ITOL, IPAR, IDID
    REAL(8) :: X, XEND, H, RPAR, TOL, RTOL, ATOL
    REAL(8), DIMENSION(NDGL) :: Y
    REAL(8), DIMENSION(LWORK) :: WORK
    INTEGER, DIMENSION(LIWORK) :: IWORK

    EXTERNAL SOLOUT, FVPOL

    ! Initial problem setup
    N = NDGL
    RPAR = 1.0D-1
    X = 0.0D0
    Y(1) = 2.0D0
    Y(2) = 0.0D0
    XEND = 2.0D0
    TOL = 1.0D-5
    ITOL = 0
    RTOL = TOL
    ATOL = TOL
    IOUT = 2
    H = 0.01D0

    IWORK(:) = 0
    WORK(:) = 0.0D0
    IWORK(8) = NRDENS
    IWORK(21) = 1
    IWORK(22) = 2

    CALL ODEX(N, FVPOL, X, Y, XEND, H, RTOL, ATOL, ITOL, SOLOUT, IOUT, &
              WORK, LWORK, IWORK, LIWORK, RPAR, IPAR, IDID)

    WRITE(6,'(1X,"X =",F5.2,"    Y =",2E18.10)') X, Y(1), Y(2)
    WRITE(6,'("       tol=",D8.2)') TOL
    WRITE(6,'(" fcn=",I5," step=",I4," accpt=",I4," rejct=",I3)') &
        IWORK(17), IWORK(18), IWORK(19), IWORK(20)

END PROGRAM VanDerPolDriver

SUBROUTINE SOLOUT(NR,XOLD,X,Y,N,CON,NCON,ICOMP,ND,RPAR,IPAR,IRTRN)
    IMPLICIT NONE
    INTEGER :: NR, N, NCON, ND, IPAR, IRTRN
    REAL(8) :: XOLD, X, RPAR
    REAL(8), DIMENSION(N) :: Y
    REAL(8), DIMENSION(NCON) :: CON
    INTEGER, DIMENSION(ND) :: ICOMP
    REAL(8), SAVE :: XOUT
    REAL(8) :: SOL1, SOL2
    REAL(8) :: CONTEX

    IF (.FALSE.) THEN
        WRITE(*,*) XOLD, RPAR, IPAR, IRTRN
    END IF

    IF (NR == 1) THEN
        WRITE(6,'(1X,"X =",F5.2,"    Y =",2E18.10,"    NSTEP =",I4)') X, Y(1), Y(2), NR-1
        XOUT = X + 0.1D0
    ELSE
        DO WHILE (X >= XOUT)
            SOL1 = CONTEX(1, XOUT, CON, NCON, ICOMP, ND)
            SOL2 = CONTEX(2, XOUT, CON, NCON, ICOMP, ND)
            WRITE(6,'(1X,"X =",F5.2,"    Y =",2E18.10,"    NSTEP =",I4)') XOUT, SOL1, SOL2, NR-1
            XOUT = XOUT + 0.1D0
        END DO
    END IF
END SUBROUTINE SOLOUT

SUBROUTINE FVPOL(N,X,Y,F,RPAR,IPAR)
    IMPLICIT NONE
    INTEGER :: N, IPAR
    REAL(8) :: X, RPAR
    REAL(8), DIMENSION(N) :: Y, F
    REAL(8) :: EPS

    IF (.FALSE.) THEN
        WRITE(*,*) X, IPAR
    END IF

    EPS = RPAR
    F(1) = Y(2)
    F(2) = ((1.0D0 - Y(1)**2) * Y(2) - Y(1)) / EPS
END SUBROUTINE FVPOL
