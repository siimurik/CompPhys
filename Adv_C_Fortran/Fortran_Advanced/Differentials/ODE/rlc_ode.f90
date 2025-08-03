!==============================================================================
! DRIVER FOR ODE() SUBROUTINE ON NON-STIFF RLC CIRCUIT
!
! Compile with:
!   gfortran -O2 rlc_ode.f90 ode.f90 -o rlc_ode
!==============================================================================
PROGRAM RLC_ODE_SOLVER
    IMPLICIT NONE
    INTEGER, PARAMETER :: DP = KIND(1.0D0)  ! Double precision

    ! --- SYSTEM DIMENSION AND WORKSPACE
    INTEGER, PARAMETER :: NEQN = 2          ! i(t), v_C(t)
    INTEGER, PARAMETER :: NWORK = 100 + 21*NEQN  ! Work array size per ODE documentation
    INTEGER, PARAMETER :: NIWORK = 5        ! Integer work array size

    ! --- VARIABLES
    REAL(DP) :: Y(NEQN), WORK(NWORK), T, TOUT, RELERR, ABSERR
    INTEGER  :: IWORK(NIWORK), IFLAG, NSTEPS
    REAL(DP) :: DT, TMAX
    
    ! --- RLC PARAMETERS (stored in common block for access by FRLC_ODE)
    REAL(DP) :: R, L, C
    COMMON /RLC_PARAMS/ R, L, C

    ! --- EXTERNAL SUBROUTINES
    EXTERNAL :: FRLC_ODE

    ! --- RLC PARAMETERS (Non-stiff: R=10Ω, L=0.1H, C=1e-4F)
    R = 10.0_DP      ! Resistance [Ω]
    L = 0.1_DP       ! Inductance [H]
    C = 1.0D-4       ! Capacitance [F]

    ! --- SOLVER CONFIGURATION
    RELERR = 1.0D-6  ! Relative tolerance
    ABSERR = 1.0D-9  ! Absolute tolerance (smaller for better accuracy)
    IFLAG  = 1       ! Initialize the solver

    ! --- INITIAL CONDITIONS (i(0)=0, v_C(0)=0)
    T     = 0.0_DP
    Y(1)  = 0.0_DP   ! Initial current [A]
    Y(2)  = 0.0_DP   ! Initial capacitor voltage [V]

    ! --- INTEGRATION PARAMETERS
    TMAX  = 0.1_DP   ! Simulate for 0.1 seconds
    DT    = 0.001_DP ! Output every 1ms
    NSTEPS = 0       ! Step counter for output

    ! --- CLEAR WORK ARRAYS
    WORK  = 0.0_DP
    IWORK = 0

    ! --- PRINT HEADER AND SOLVE SYSTEM
    CALL PRINT_HEADER(R, L, C, RELERR, ABSERR)
    CALL SOLVE_RLC_SYSTEM(Y, T, TMAX, DT, RELERR, ABSERR, IFLAG, WORK, IWORK, NSTEPS)
    CALL PRINT_FINAL_STATS(Y, NSTEPS, TMAX, RELERR, ABSERR)

END PROGRAM RLC_ODE_SOLVER

!==============================================================================
! PRINT PROGRAM HEADER AND PARAMETERS
!==============================================================================
SUBROUTINE PRINT_HEADER(R, L, C, RELERR, ABSERR)
    IMPLICIT NONE
    INTEGER, PARAMETER :: DP = KIND(1.0D0)
    REAL(DP), INTENT(IN) :: R, L, C, RELERR, ABSERR
    
    WRITE (*, '(A)') 'RLC Circuit Solution using ODE() Subroutine'
    WRITE (*, '(A, F4.1, A, F5.3, A, E8.1, A)') &
        'Parameters: R=', R, 'Ω, L=', L, 'H, C=', C, 'F'
    WRITE (*, '(A, E8.1, A, E8.1)') &
        'Tolerances: RELERR=', RELERR, ', ABSERR=', ABSERR
    WRITE (*, '(A)') REPEAT('-', 60)
END SUBROUTINE PRINT_HEADER

!==============================================================================
! SOLVE RLC SYSTEM WITH ODE INTEGRATION
!==============================================================================
SUBROUTINE SOLVE_RLC_SYSTEM(Y, T, TMAX, DT, RELERR, ABSERR, IFLAG, WORK, IWORK, NSTEPS)
    IMPLICIT NONE
    INTEGER, PARAMETER :: DP = KIND(1.0D0)
    INTEGER, PARAMETER :: NEQN = 2
    
    REAL(DP), INTENT(INOUT) :: Y(NEQN), T, RELERR, ABSERR
    REAL(DP), INTENT(IN) :: TMAX, DT
    REAL(DP), INTENT(INOUT) :: WORK(*)
    INTEGER, INTENT(INOUT) :: IFLAG, IWORK(*), NSTEPS
    
    REAL(DP) :: TOUT
    EXTERNAL :: FRLC_ODE
    
    ! Print initial condition
    WRITE (*, '(A, F6.3, A, E18.10, A, E18.10, A, I6)') &
        't =', T, '    I =', Y(1), '    V =', Y(2), '    NSTEP =', NSTEPS

    ! Integration loop with output at regular intervals
    TOUT = DT
    DO WHILE (T < TMAX)
        ! Set target time (don't overshoot TMAX)
        TOUT = MIN(TOUT, TMAX)
        
        ! Call ODE solver
        CALL ODE(FRLC_ODE, NEQN, Y, T, TOUT, RELERR, ABSERR, IFLAG, WORK, IWORK)
        
        ! Check for errors
        CALL CHECK_ODE_STATUS(IFLAG)
        
        ! Count steps and print solution
        NSTEPS = NSTEPS + 1
        WRITE (*, '(A, F6.3, A, E18.10, A, E18.10, A, I6)') &
            't =', T, '    I =', Y(1), '    V =', Y(2), '    NSTEP =', NSTEPS
        
        ! Prepare for next output time
        TOUT = TOUT + DT
        IFLAG = 1  ! Continue integration
    END DO
END SUBROUTINE SOLVE_RLC_SYSTEM

!==============================================================================
! CHECK ODE SOLVER STATUS AND HANDLE ERRORS
!==============================================================================
SUBROUTINE CHECK_ODE_STATUS(IFLAG)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: IFLAG
    
    IF (IFLAG < 0) THEN
        WRITE (*, '(A, I3)') 'ERROR: Integration failed with IFLAG =', IFLAG
        STOP
    ELSE IF (IFLAG == 3) THEN
        ! Tolerances were adjusted, continue
        WRITE (*, '(A)') 'WARNING: Error tolerances adjusted'
        IFLAG = 1  ! Reset flag to continue
    ELSE IF (IFLAG == 4) THEN
        WRITE (*, '(A)') 'ERROR: Too many steps (>500)'
        STOP
    ELSE IF (IFLAG == 5) THEN
        WRITE (*, '(A)') 'ERROR: Equations appear to be stiff'
        STOP
    ELSE IF (IFLAG == 6) THEN
        WRITE (*, '(A)') 'ERROR: Invalid input parameters'
        STOP
    END IF
END SUBROUTINE CHECK_ODE_STATUS

!==============================================================================
! PRINT FINAL STATISTICS
!==============================================================================
SUBROUTINE PRINT_FINAL_STATS(Y, NSTEPS, TMAX, RELERR, ABSERR)
    IMPLICIT NONE
    INTEGER, PARAMETER :: DP = KIND(1.0D0)
    REAL(DP), INTENT(IN) :: Y(2), TMAX, RELERR, ABSERR
    INTEGER, INTENT(IN) :: NSTEPS
    
    WRITE (*, '(A)') REPEAT('-', 60)
    WRITE (*, '(A, E18.10, A, E18.10)') 'Final state: I =', Y(1), ', V_C =', Y(2)
    WRITE (*, '(A, I6, A, F6.3, A)') 'Total output steps:', NSTEPS, ' over ', TMAX, ' seconds'
    WRITE (*, '(A, E8.1, A, E8.1)') 'Final tolerances: RELERR=', RELERR, ', ABSERR=', ABSERR
END SUBROUTINE PRINT_FINAL_STATS

!==============================================================================
! RIGHT-HAND SIDE OF RLC CIRCUIT ODE
!   dy1/dt = (-R/L)*y1 - (1/L)*y2 + Vin/L
!   dy2/dt = (1/C)*y1
!
! This subroutine must match the interface expected by ODE():
!   subroutine f ( t, y, yp )
!==============================================================================
SUBROUTINE FRLC_ODE(T, Y, YP)
    IMPLICIT NONE
    INTEGER, PARAMETER :: DP = KIND(1.0D0)  ! Double precision
    REAL(DP), INTENT(IN) :: T, Y(2)
    REAL(DP), INTENT(OUT) :: YP(2)
    
    ! --- RLC PARAMETERS FROM COMMON BLOCK
    REAL(DP) :: R, L, C
    COMMON /RLC_PARAMS/ R, L, C
    
    REAL(DP) :: Vin

    ! --- AC SOURCE: 5V amplitude, 100 rad/s (≈15.9 Hz)
    Vin = 5.0_DP * SIN(100.0_DP * T)

    ! --- DIFFERENTIAL EQUATIONS
    YP(1) = (-R/L) * Y(1) - (1.0_DP/L) * Y(2) + Vin/L  ! di/dt
    YP(2) = (1.0_DP/C) * Y(1)                           ! dv_C/dt
    
END SUBROUTINE FRLC_ODE
