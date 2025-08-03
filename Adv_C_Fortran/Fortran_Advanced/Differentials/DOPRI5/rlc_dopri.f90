!==============================================================================
! DRIVER FOR DOPRI5 ON NON-STIFF RLC CIRCUIT
!
! Compile with:
!   gfortran -O2 dr_rlc_dopri5.f90 dopri5.f90 -o rlc_dopri5
!==============================================================================
PROGRAM RLC_DOPRI5
    IMPLICIT NONE
    INTEGER, PARAMETER :: DP = KIND(1.0D0)  ! Double precision

    ! --- SYSTEM DIMENSION AND WORKSPACE
    INTEGER, PARAMETER :: NDGL = 2          ! i(t), v_C(t)
    INTEGER, PARAMETER :: NRDENS = 2        ! Dense output for both variables
    INTEGER, PARAMETER :: LWORK = 8*NDGL + 5*NRDENS + 21
    INTEGER, PARAMETER :: LIWORK = NRDENS + 21

    ! --- VARIABLES
    REAL(DP) :: Y(NDGL), WORK(LWORK), RPAR(3), T, TEND, RTOL, ATOL
    INTEGER  :: IWORK(LIWORK), ITOL, IOUT, IDID, IPAR, N

    ! --- EXTERNAL SUBROUTINES
    EXTERNAL :: FRLC, SOLOUT_RLC

    ! --- RLC PARAMETERS (Non-stiff: R=10Ω, L=0.1H, C=1e-4F)
    ! Non-stiff RLC (R=10, L=0.1, C=1e-4) (For stiff use RODAS)
    RPAR(1) = 10.0_DP    ! Resistance [Ω]
    RPAR(2) = 0.1_DP     ! Inductance [H]
    RPAR(3) = 1.0D-4     ! Capacitance [F]
    N = NDGL             ! System dimension

    ! --- SOLVER CONFIGURATION
    IOUT  = 2            ! Dense output enabled
    ITOL  = 0            ! Scalar tolerances
    RTOL  = 1.0D-6       ! Relative tolerance
    ATOL  = RTOL         ! Absolute tolerance

    ! --- INITIAL CONDITIONS (i(0)=0, v_C(0)=0)
    T     = 0.0_DP
    Y(1)  = 0.0_DP       ! Initial current [A]
    Y(2)  = 0.0_DP       ! Initial capacitor voltage [V]

    ! --- INTEGRATION PARAMETERS
    TEND  = 0.1_DP       ! Simulate for 0.1 seconds
    WORK  = 0.0_DP       ! Clear work arrays
    IWORK = 0            ! Clear integer options

    ! --- DENSE OUTPUT FOR BOTH VARIABLES (i and v_C)
    IWORK(5) = NRDENS
    IWORK(21) = 1        ! Output Y(1) (current)
    IWORK(22) = 2        ! Output Y(2) (voltage)

    ! --- CALL DOPRI5
    CALL DOPRI5(N, FRLC, T, Y, TEND, &
                RTOL, ATOL, ITOL, &
                SOLOUT_RLC, IOUT, &
                WORK, LWORK, IWORK, LIWORK, RPAR, IPAR, IDID)

    ! --- PRINT FINAL SOLUTION
    WRITE (*, '(A, F6.3, A, E18.10, A, E18.10, A, I6)') &
        't =', TEND, '    I =', Y(1), '    V =', Y(2), '    NSTEP =', IWORK(18)
    WRITE (*, '(A, 2E18.10)') 'Final state: I, V_C =', Y(1), Y(2)
                
    ! --- PRINT STATISTICS (including NSTEP = IWORK(18))
    WRITE (*, '(A, D8.2, A, I5, A, I4, A, I4, A, I3)') &
        'Tol=', RTOL, '  Fcn=', IWORK(17), ' Step=', IWORK(18), &
        ' Accpt=', IWORK(19), ' Rejct=', IWORK(20)

END PROGRAM RLC_DOPRI5

!==============================================================================
! OUTPUT ROUTINE - PRINTS SOLUTION AT 1ms INTERVALS
!==============================================================================
SUBROUTINE SOLOUT_RLC(NR, XOLD, T, Y, N, CON, ICOMP, ND, RPAR, IPAR, IRTRN)
    IMPLICIT NONE
    INTEGER, PARAMETER :: DP = KIND(1.0D0)  ! Double precision
    INTEGER, INTENT(IN) :: NR, N, ND, ICOMP(ND), IPAR
    REAL(DP), INTENT(IN) :: XOLD, T, Y(N), CON(5*ND), RPAR(3)
    INTEGER, INTENT(INOUT) :: IRTRN
    REAL(DP) :: CONTD5, DT

    REAL(DP) :: TOUT, Y1, Y2
    COMMON /INTERN_RLC/ TOUT  ! Persists between calls

    DT = 0.001_DP
    IF (NR == 1) THEN
        ! First call: print initial condition and NSTEP=0
        WRITE (*, '(A, F6.3, A, E18.10, A, E18.10, A, I6)') &
            't =', T, '    I =', Y(1), '    V =', Y(2), '    NSTEP =', 0
        TOUT = DT  ! First output at 1ms
    ELSE
        ! Intermediate steps: interpolate and print with NSTEP=NR-1
        DO WHILE (TOUT < T - SPACING(T))
            Y1 = CONTD5(1, TOUT, CON, ICOMP, ND)  ! Interpolated current
            Y2 = CONTD5(2, TOUT, CON, ICOMP, ND)  ! Interpolated voltage
            WRITE (*, '(A, F6.3, A, E18.10, A, E18.10, A, I6)') &
                't =', TOUT, '    I =', Y1, '    V =', Y2, '    NSTEP =', NR-1
            TOUT = TOUT + DT  ! Increment by 1ms
        END DO
    END IF
END SUBROUTINE SOLOUT_RLC

!==============================================================================
! RIGHT-HAND SIDE OF RLC CIRCUIT ODE
!   dy1/dt = (-R/L)*y1 - (1/L)*y2 + Vin/L
!   dy2/dt = (1/C)*y1
!==============================================================================
SUBROUTINE FRLC(N, T, Y, F, RPAR, IPAR)
    IMPLICIT NONE
    INTEGER, PARAMETER :: DP = KIND(1.0D0)  ! Double precision
    INTEGER, INTENT(IN) :: N, IPAR
    REAL(DP), INTENT(IN) :: T, Y(N), RPAR(3)
    REAL(DP), INTENT(OUT) :: F(N)
    REAL(DP) :: R, L, C, Vin

    R = RPAR(1)
    L = RPAR(2)
    C = RPAR(3)
    Vin = 5.0_DP * SIN(100.0_DP * T)  ! AC source: 5V, 100Hz

    F(1) = (-R/L) * Y(1) - (1.0_DP/L) * Y(2) + Vin/L
    F(2) = (1.0_DP/C) * Y(1)
END SUBROUTINE FRLC