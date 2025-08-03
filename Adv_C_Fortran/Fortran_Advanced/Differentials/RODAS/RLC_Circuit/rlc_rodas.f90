!==============================================================================
! DRIVER FOR RODAS SOLVER ON RLC CIRCUIT ODE
!
! Compile with:
!   $ gfortran -O2 rlc_rodas.f90 rodas.f90 decsol.f90 dc_decsol.f90 -o rlc_solver
!   $ ./rlc_solver > rlc_output.dat
!==============================================================================
PROGRAM RLC_RODAS
    IMPLICIT NONE
    INTEGER, PARAMETER :: DP = KIND(1.0D0)  ! Double precision

    ! --- PROBLEM DIMENSION AND WORKSPACE
    INTEGER, PARAMETER :: ND = 2            ! i(t), v_C(t)
    INTEGER, PARAMETER :: LWORK = 2*ND*ND + 14*ND + 20
    INTEGER, PARAMETER :: LIWORK = ND + 20

    ! --- VARIABLES
    REAL(DP) :: Y(ND), WORK(LWORK), RPAR(3), X, XEND, H, RTOL, ATOL
    INTEGER  :: IWORK(LIWORK), ITOL, IFCN, IJAC, MLJAC, MUJAC, IDFX, &
                IMAS, MLMAS, MUMAS, IOUT, IDID, IPAR, N

    ! --- EXTERNAL SUBROUTINES
    EXTERNAL :: FRLC, JRLC, SOLOUT_RLC

    ! --- RLC PARAMETERS (R, L, C)
    ! Stiff RLC (R=1e4, L=0.1, C=1e-6)
    RPAR(1) = 100.0_DP   ! Resistance (Ohms)
    RPAR(2) = 0.1_DP     ! Inductance (Henry)
    RPAR(3) = 1.0D-4     ! Capacitance (Farad)
    N = ND               ! System dimension

    ! --- SOLVER CONFIGURATION
    IFCN   = 0           ! Autonomous system (no explicit X-dependence)
    IJAC   = 1           ! Analytical Jacobian provided
    MLJAC  = N           ! Full Jacobian
    IMAS   = 0           ! No mass matrix
    IOUT   = 1           ! Enable output routine

    ! --- INITIAL CONDITIONS (i(0) = 0, v_C(0) = 0)
    X      = 0.0_DP
    Y(1)   = 0.0_DP      ! Initial current
    Y(2)   = 0.0_DP      ! Initial capacitor voltage

    ! --- INTEGRATION PARAMETERS
    XEND   = 0.1_DP      ! Simulate for 0.1 seconds
    RTOL   = 1.0D-6      ! Relative tolerance
    ATOL   = 1.0D-6      ! Absolute tolerance
    ITOL   = 0           ! Scalar tolerances
    H      = 1.0D-6      ! Initial step guess

    ! --- CLEAR WORK ARRAYS
    WORK  = 0.0_DP
    IWORK = 0

    ! --- CALL RODAS SOLVER
    CALL RODAS(N, FRLC, IFCN, X, Y, XEND, H, &
               RTOL, ATOL, ITOL, &
               JRLC, IJAC, MLJAC, MUJAC, FRLC, IDFX, &
               FRLC, IMAS, MLMAS, MUMAS, &
               SOLOUT_RLC, IOUT, &
               WORK, LWORK, IWORK, LIWORK, RPAR, IPAR, IDID)

    ! --- PRINT FINAL SOLUTION
    WRITE (*, '(A, F8.5, A, F10.6, A, F10.6, A, I6)') &
                'Time(s)=', XEND, '  I(A)=', Y(1), '  V_C(V)=', Y(2),&
                '    NSTEP =', IWORK(16)

    ! --- PRINT FINAL STATISTICS
    WRITE (*, '(A, D8.2)') '       rtol=', RTOL
    WRITE (*, '(A, I7, A, I6, A, I6, A, I6, A, I4, A, I6, A, I7)') &
        ' fcn=', IWORK(14), ' jac=', IWORK(15), ' step=', IWORK(16), &
        ' accpt=', IWORK(17), ' rejct=', IWORK(18), ' dec=', IWORK(19), &
        ' sol=', IWORK(20)

END PROGRAM RLC_RODAS

!==============================================================================
! RIGHT-HAND SIDE OF RLC ODE SYSTEM
!   dy1/dt = (-R/L)*y1 - (1/L)*y2 + Vin/L
!   dy2/dt = (1/C)*y1
!==============================================================================
SUBROUTINE FRLC(N, X, Y, F, RPAR, IPAR)
    IMPLICIT NONE
    INTEGER, PARAMETER :: DP = KIND(1.0D0)  ! Double precision
    INTEGER, INTENT(IN) :: N, IPAR
    REAL(DP), INTENT(IN) :: X, Y(N), RPAR(3)
    REAL(DP), INTENT(OUT) :: F(N)
    REAL(DP) :: R, L, C, Vin

    R = RPAR(1)
    L = RPAR(2)
    C = RPAR(3)
    Vin = 5.0_DP * SIN(100.0_DP * X)  ! AC voltage source (5V, 100Hz)

    F(1) = (-R/L) * Y(1) - (1.0_DP/L) * Y(2) + Vin/L
    F(2) = (1.0_DP/C) * Y(1)

END SUBROUTINE FRLC

!==============================================================================
! JACOBIAN OF RLC SYSTEM
!   J = [ -R/L    -1/L ]
!       [  1/C      0  ]
!==============================================================================
SUBROUTINE JRLC(N, X, Y, DFY, LDFY, RPAR, IPAR)
    IMPLICIT NONE
    INTEGER, PARAMETER :: DP = KIND(1.0D0)  ! Double precision
    INTEGER, INTENT(IN) :: N, LDFY, IPAR
    REAL(DP), INTENT(IN) :: X, Y(N), RPAR(3)
    REAL(DP), INTENT(OUT) :: DFY(LDFY, N)
    REAL(DP) :: R, L, C

    R = RPAR(1)
    L = RPAR(2)
    C = RPAR(3)

    DFY(1, 1) = -R/L
    DFY(1, 2) = -1.0_DP/L
    DFY(2, 1) = 1.0_DP/C
    DFY(2, 2) = 0.0_DP

END SUBROUTINE JRLC

!==============================================================================
! OUTPUT ROUTINE - PRINTS CURRENT AND VOLTAGE AT FIXED INTERVALS
!==============================================================================
SUBROUTINE SOLOUT_RLC(NR, XOLD, X, Y, CONT, LRC, N, RPAR, IPAR, IRTRN)
    IMPLICIT NONE
    INTEGER, PARAMETER :: DP = KIND(1.0D0)  ! Double precision
    INTEGER, INTENT(IN) :: NR, LRC, N, IPAR
    REAL(DP), INTENT(IN) :: XOLD, X, Y(N), RPAR(3)
    INTEGER, INTENT(INOUT) :: IRTRN
    REAL(DP), INTENT(IN) :: CONT(LRC)
    REAL(DP) ::  CONTRO

    REAL(DP) :: XOUT, Y1, Y2, DT
    COMMON /INTERN_RLC/ XOUT  ! Persists between calls

    DT = 0.001_DP  ! Output every 1ms

    IF (NR == 1) THEN
        WRITE (*, '(A, F8.5, A, F10.6, A, F10.6, A, I6)') &
            'Time(s)=', X, '  I(A)=', Y(1), '  V_C(V)=', Y(2),  '    NSTEP =', 0
        XOUT = X + DT
    ELSE
        DO WHILE (XOUT < X - SPACING(X))
            Y1 = CONTRO(1, XOUT, CONT, LRC)  ! Interpolated current
            Y2 = CONTRO(2, XOUT, CONT, LRC)  ! Interpolated voltage
            WRITE (*, '(A, F8.5, A, F10.6, A, F10.6, A, I6)') &
                'Time(s)=', XOUT, '  I(A)=', Y1, '  V_C(V)=', Y2,'    NSTEP =', NR-1
            XOUT = XOUT + DT
        END DO
    END IF

END SUBROUTINE SOLOUT_RLC