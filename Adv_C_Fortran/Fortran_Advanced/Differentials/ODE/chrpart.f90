!==============================================================================
! DRIVER FOR ODE() SUBROUTINE ON CHARGED PARTICLE IN EARTH'S MAGNETIC FIELD
!==============================================================================
PROGRAM ODE_SOLVER
    IMPLICIT NONE
    ! --- SYSTEM DIMENSION AND WORKSPACE
    INTEGER, PARAMETER :: NEQN = 6
    INTEGER, PARAMETER :: NWORK = 100 + 21*NEQN
    INTEGER, PARAMETER :: NIWORK = 5

    ! --- VARIABLES
    DOUBLE PRECISION :: Y(NEQN), WORK(NWORK), T, RELERR, ABSERR
    INTEGER  :: IWORK(NIWORK), IFLAG, NSTEPS
    DOUBLE PRECISION :: DT, TMAX
    DOUBLE PRECISION :: q, m
    LOGICAL :: PRINT_ALL
    CHARACTER(LEN=100) :: OUTFILE

    ! --- PHYSICAL CONSTANTS (Common block for func)
    DOUBLE PRECISION :: qm, R_maa, B0
    COMMON /constants/ qm, R_maa, B0

    ! --- EXTERNAL SUBROUTINES
    EXTERNAL :: func

    ! Physical constants (SI units)
    q = 1.602D-19          ! Coulombs
    m = 1.673D-27          ! kg (proton mass)
    qm = q/m               ! C/kg
    B0 = 3.07D-5           ! Tesla (Earth's dipole moment coefficient, approx value)
    R_maa = 6.378137D6     ! Earth's radius in meters

    ! Numerical integration parameters
    RELERR = 1.D-8
    ABSERR = 1.D-8
    T  = 0.0D0
    DT = 0.01D0
    NSTEPS = 25000
    TMAX = DT * NSTEPS
    IFLAG = 1

    ! Initial conditions (SI units)
    Y(1) = 3.0D7            ! x (meters) = 30000 km above Earth's surface
    Y(2) = 0.0D0            ! y (meters)
    Y(3) = 0.0D0            ! z (meters)
    Y(4) = 0.0D0            ! vx (m/s)
    Y(5) = 1.0D7            ! vy (m/s) = 10000 km/s
    Y(6) = 2.0D7            ! vz (m/s) = 20000 km/s

    ! --- Option for full output to file
    PRINT_ALL = .TRUE.                ! Set to .FALSE. for no full file output
    OUTFILE = 'timesteps_output.txt'  ! Output filename

    ! Print header
    CALL PRINT_HEADER(RELERR, ABSERR)

    ! Solve system
    CALL SOLVE_ODE_SYSTEM(Y, T, TMAX, DT, RELERR, ABSERR, IFLAG, WORK, IWORK, NSTEPS, PRINT_ALL, OUTFILE)

    ! Print final state
    CALL PRINT_FINAL_STATS(Y, NSTEPS, TMAX, RELERR, ABSERR)

END PROGRAM ODE_SOLVER

!==============================================================================
SUBROUTINE PRINT_HEADER(RELERR, ABSERR)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: RELERR, ABSERR
    WRITE (*, '(A)') 'Charged Particle in Earth''s Dipole Magnetic Field (RKF45)'
    WRITE (*, '(A, E8.1, A, E8.1)') 'Tolerances: RELERR=', RELERR, ', ABSERR=', ABSERR
    WRITE (*, '(A)') REPEAT('-', 60)
END SUBROUTINE PRINT_HEADER

!==============================================================================
SUBROUTINE SOLVE_ODE_SYSTEM(Y, T, TMAX, DT, RELERR, ABSERR, IFLAG, WORK, IWORK, NSTEPS, PRINT_ALL, OUTFILE)
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NEQN = 6

    DOUBLE PRECISION, INTENT(INOUT) :: Y(NEQN), T, RELERR, ABSERR
    DOUBLE PRECISION, INTENT(IN) :: TMAX, DT
    DOUBLE PRECISION, INTENT(INOUT) :: WORK(*)
    INTEGER, INTENT(INOUT) :: IFLAG, IWORK(*)
    INTEGER, INTENT(IN) :: NSTEPS
    LOGICAL, INTENT(IN) :: PRINT_ALL
    CHARACTER(LEN=*), INTENT(IN) :: OUTFILE

    DOUBLE PRECISION :: TOUT
    INTEGER :: ISTEP, OUTUNIT
    EXTERNAL :: func

    ! File output handling
    IF (PRINT_ALL) THEN
        OPEN(NEWUNIT=OUTUNIT, FILE=OUTFILE, STATUS='REPLACE', ACTION='WRITE')
        WRITE(OUTUNIT,'(A)') 't     x      y      z      vx      vy      vz'
    END IF

    ! Print initial condition
    WRITE (*, '(A, F10.3, A, 3E15.6, A, 3E12.4)') &
        't =', T, '   pos (m):', Y(1), Y(2), Y(3), '   vel (m/s):', Y(4), Y(5), Y(6)
    IF (PRINT_ALL) WRITE(OUTUNIT,'(F10.5,1X,3E15.6,1X,3E12.4)') T, Y(1), Y(2), Y(3), Y(4), Y(5), Y(6)

    TOUT = DT
    DO ISTEP = 1, NSTEPS
        TOUT = MIN(T + DT, TMAX)
        CALL ODE(func, NEQN, Y, T, TOUT, RELERR, ABSERR, IFLAG, WORK, IWORK)
        CALL CHECK_ODE_STATUS(IFLAG)

        ! Print every 1000 steps to reduce output
        IF (MOD(ISTEP, 1000) == 0 .OR. ISTEP == 1 .OR. ISTEP == NSTEPS) THEN
            WRITE (*, '(A, F10.3, A, 3E15.6, A, 3E12.4, A, I8)') &
                't =', T, '   pos (m):', Y(1), Y(2), Y(3), '   vel (m/s):', Y(4), Y(5), Y(6), '   step=', ISTEP
        END IF
        IF (PRINT_ALL) WRITE(OUTUNIT,'(F10.5,1X,3E15.6,1X,3E12.4)') T, Y(1), Y(2), Y(3), Y(4), Y(5), Y(6)
    END DO

    IF (PRINT_ALL) CLOSE(OUTUNIT)
END SUBROUTINE SOLVE_ODE_SYSTEM

!==============================================================================
SUBROUTINE CHECK_ODE_STATUS(IFLAG)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: IFLAG
    IF (IFLAG < 0) THEN
        WRITE (*, '(A, I3)') 'ERROR: Integration failed with IFLAG =', IFLAG
        STOP
    ELSE IF (IFLAG == 3) THEN
        WRITE (*, '(A)') 'WARNING: Error tolerances adjusted'
        IFLAG = 1
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
SUBROUTINE PRINT_FINAL_STATS(Y, NSTEPS, TMAX, RELERR, ABSERR)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: Y(6), TMAX, RELERR, ABSERR
    INTEGER, INTENT(IN) :: NSTEPS
    WRITE (*, '(A)') REPEAT('-', 60)
    WRITE (*, '(A, 3E15.6, A, 3E12.4)') 'Final state: pos (m):', Y(1), Y(2), Y(3), '   vel (m/s):', Y(4), Y(5), Y(6)
    WRITE (*, '(A, I8, A, F10.3, A)') 'Total output steps:', NSTEPS, ' over ', TMAX, ' seconds'
    WRITE (*, '(A, E8.1, A, E8.1)') 'Final tolerances: RELERR=', RELERR, ', ABSERR=', ABSERR
END SUBROUTINE PRINT_FINAL_STATS

!==============================================================================
SUBROUTINE calc_magn_field(x, y, z, B)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: x, y, z
    DOUBLE PRECISION, INTENT(OUT) :: B(3)
    DOUBLE PRECISION :: r, mx, my, mz, m, mu0_4pi, r5, r3, rvec(3)
    ! Earth's magnetic moment (z-direction)
    mu0_4pi = 1D-7
    m = 7.94D22         ! Earth's magnetic moment, A·m^2
    mx = 0D0
    my = 0D0
    mz = m

    r = SQRT(x**2 + y**2 + z**2)
    IF (r < 1.D-12) THEN
        B = 0.D0
        RETURN
    END IF

    rvec(1) = x
    rvec(2) = y
    rvec(3) = z

    r3 = r**3
    r5 = r**5

    ! 3(m . r) r - m r^2
    B(1) = mu0_4pi * (3.D0 * (mx*x + my*y + mz*z)*x / r5 - mx / r3)
    B(2) = mu0_4pi * (3.D0 * (mx*x + my*y + mz*z)*y / r5 - my / r3)
    B(3) = mu0_4pi * (3.D0 * (mx*x + my*y + mz*z)*z / r5 - mz / r3)
END SUBROUTINE calc_magn_field

!==============================================================================
SUBROUTINE func(t, y, dy)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: t, y(6)
    DOUBLE PRECISION, INTENT(OUT) :: dy(6)
    DOUBLE PRECISION :: B(3), qm, R_maa, B0
    COMMON /constants/ qm, R_maa, B0

    ! Silence unused variable(s)
    IF (.FALSE.) THEN
        WRITE(*,*) t
    END IF

    ! Position derivatives = velocity
    dy(1) = y(4)
    dy(2) = y(5)
    dy(3) = y(6)

    ! Magnetic field at current position (in meters, not normalized)
    CALL calc_magn_field(y(1), y(2), y(3), B)

    ! Lorentz force: dv/dt = q/m (v x B)
    dy(4) = qm * (y(5)*B(3) - y(6)*B(2))
    dy(5) = qm * (y(6)*B(1) - y(4)*B(3))
    dy(6) = qm * (y(4)*B(2) - y(5)*B(1))
END SUBROUTINE func