! ===============================================================================
! CHARGED PARTICLE MOTION IN EARTH'S MAGNETIC FIELD (RKF45 driver)
! ===============================================================================

PROGRAM charged_particle_motion
    IMPLICIT NONE

    INTEGER, PARAMETER :: neqn = 6, nt = 26000
    DOUBLE PRECISION, PARAMETER :: dt = 0.01D0
    DOUBLE PRECISION, PARAMETER :: relerr = 1.D-8, abserr = 1.D-8

    DOUBLE PRECISION :: y(neqn), yp(neqn)
    DOUBLE PRECISION :: t, tout
    DOUBLE PRECISION :: q, m, qm, B0, R_maa
    INTEGER :: i, flag
    COMMON /constants/ qm, R_maa, B0

    EXTERNAL func
    EXTERNAL rkf45

    ! Physical constants (SI)
    q = 1.602D-19        ! Elementary charge (C)
    m = 1.673D-27        ! Proton mass (kg)
    qm = q/m             ! Charge-to-mass ratio (C/kg)
    R_maa = 6.378137D6   ! Earth radius (m)
    B0 = 3.07D-5         ! Not used (for compatibility with func signature)

    ! Initial conditions (SI)
    y(1) = 3.0D7        ! x (meters) = 30000 km above Earth's surface
    y(2) = 0.0D0        ! y (meters)
    y(3) = 0.0D0        ! z (meters)
    y(4) = 0.0D0        ! vx (m/s)
    y(5) = 1.0D7        ! vy (m/s) = 10000 km/s
    y(6) = 2.0D7        ! vz (m/s) = 20000 km/s

    t = 0.0D0
    flag = 1

    OPEN(10, FILE="proton.dat", status='replace')
    WRITE(10,'(A)') 't     x      y      z      vx      vy      vz'
    WRITE(10,'(F10.5,1X,3E15.7,1X,3E12.4)') t, y(1), y(2), y(3), y(4), y(5), y(6)

    DO i = 1, nt
        tout = t + dt
        CALL rkf45(func, neqn, y, yp, t, tout, relerr, abserr, flag)
        IF (flag < 0 .OR. flag == 6 .OR. flag == 8) THEN
            PRINT *, 'RKF45 error, flag =', flag
            EXIT
        END IF
        WRITE(10,'(F10.5,1X,3E15.7,1X,3E12.4)') t, y(1), y(2), y(3), y(4), y(5), y(6)
        IF (MOD(i, 1000) == 0 .OR. i == 1 .OR. i == nt) THEN
            WRITE (*, '(A, F10.3, A, 3E14.6, A, 3E12.4, A, I8)') &
                't =', t, '   pos (m): [', y(1), y(2), y(3), &
                ' ]  vel (m/s): [', y(4), y(5), y(6), ' ]  step=', i
        END IF
        ! flag instruction options for RKF45:
        ! flag = 1  ! Integrate from the current t to the new tout (target time), 
                    ! and do whatever steps you need to land exactly on tout.
        !-------------------------------------------------------------------------
         flag = 2  ! Continue integrating from where you left off, 
                    ! but you do NOT have to hit tout exactly — 
                    ! just keep going with whatever step size you think is best.
    END DO

    CLOSE(10)
END PROGRAM charged_particle_motion

! ===============================================================================
SUBROUTINE calc_magn_field(x, y, z, B)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: x, y, z
    DOUBLE PRECISION, INTENT(OUT) :: B(3)
    DOUBLE PRECISION :: r, mx, my, mz, m, mu0_4pi, r5, r3
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

    r3 = r**3
    r5 = r**5

    ! 3(m . r) r - m r^2, dipole field
    B(1) = mu0_4pi * (3.D0 * (mx*x + my*y + mz*z)*x / r5 - mx / r3)
    B(2) = mu0_4pi * (3.D0 * (mx*x + my*y + mz*z)*y / r5 - my / r3)
    B(3) = mu0_4pi * (3.D0 * (mx*x + my*y + mz*z)*z / r5 - mz / r3)
END SUBROUTINE calc_magn_field

! ===============================================================================
SUBROUTINE func(t, y, yp)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: t, y(6)
    DOUBLE PRECISION, INTENT(OUT) :: yp(6)
    DOUBLE PRECISION :: B(3), qm, R_maa, B0
    COMMON /constants/ qm, R_maa, B0

    ! Silence unused variable(s)
    IF (.FALSE.) THEN
        WRITE(*,*) t
    END IF

    ! Position derivatives = velocity
    yp(1) = y(4)
    yp(2) = y(5)
    yp(3) = y(6)

    ! Magnetic field at current position (in meters)
    CALL calc_magn_field(y(1), y(2), y(3), B)

    ! Lorentz force: dv/dt = q/m (v x B)
    yp(4) = qm * (y(5)*B(3) - y(6)*B(2))
    yp(5) = qm * (y(6)*B(1) - y(4)*B(3))
    yp(6) = qm * (y(4)*B(2) - y(5)*B(1))
END SUBROUTINE func