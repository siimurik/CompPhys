! Compile and execute with:
!   $ gfortran -O2 pm_airres.f90 rkf45.f90 -o pm
!   $ ./pm
program main
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    call timestamp()
    write(*, '(a)') ""
    write(*, '(a)') "PROJECTILE MOTION WITH AIR RESISTANCE" 
    write(*, '(a)') "USING THE RUNGE-KUTTA-FEHLBERG METHOD" 
    write(*, '(a)') "" 

    call pmwar()

    write (*, '(a)') ""
    write (*, '(a)') "RKF45 SOLVER" 
    write (*, '(a)') "Normal end of execution." 
    write (*, '(a)') "" 
    call timestamp()

    stop 0

end program main

subroutine pmwar()
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: neqn = 4 ! number of equations
    real(kind=dp) :: abserr, relerr
    integer :: iflag, i_step, n_step
    external :: func
    real(kind=dp) :: t, t_out, t_start,t_stop
    real(kind=dp) :: v_0, theta, dt
    real(kind=dp) :: y(neqn), yp(neqn)

    write(*,'(a)') "" 
    write(*,'(a)') "Solve a 2nd order ODE " 
    write(*,'(a)') "  d^2x       dx           " 
    write(*,'(a)') "m ----  +  d --  +  mg = 0" 
    write(*,'(a)') "   dt        dt           " 
    write(*,'(a)') "as a system of 1st order ODE" 
    write(*,'(a)') "[ dx/dt   = v_x" 
    write(*,'(a)') "{ dy/dt   = v_y" 
    write(*,'(a)') "{ dv_x/dt = -d/m*v_x"
    write(*,'(a)') "[ dv_y/dt = -d/m*v_y - g"
    
    ! Error tolerances
    abserr = sqrt(epsilon(abserr))
    relerr = sqrt(epsilon(relerr))

    ! Initialize integration
    iflag = 1

    ! Time span of integration
    t_start =   0.0D0
    t_stop  = 100.0D0
    n_step  = 10000

    t = 0.0D0
    t_out = 0.0D0
    dt = 0.01D0

    ! Initial conditions
    y(1)  = 0.0D0 ! x(1) = 0
    y(2)  = 0.0D0 ! y(1) = 0
    v_0   = 30.D0 ! m/s
    theta = 30.D0 ! degrees
    y(3)  = v_0*cosd(theta) ! v_x
    y(4)  = v_0*sind(theta) ! v_y

    call func(t, y, yp)

    write(*, '(a)') ""
    write(*, '(a)') "  FLAG         T            Y(1)            Y(2)            Y(3)            Y(4)"
    write(*, '(a)') ""
    write ( *, '(i4, 2x, 4g16.8, 4g16.8, 4g16.8, 4g16.8)' ) iflag, t, y(1), y(2), y(3), y(4)

    !do i_step = 1, 10
    !    t_out = t + dt
    !    call rkf45(func, neqn, y, yp, t, t_out, relerr, abserr, iflag)
    !    if (iflag.ne.2) then
    !        iflag = 2
    !    end if
    !    if (y(2).lt.0) stop
    !    write ( *, '(i4, 2x, 4g16.8, 4g16.8, 4g16.8, 4g16.8)' ) iflag, t, y(1), y(2), y(3), y(4)
    !end do

    open(unit=11, file='out.dat')
    write ( 11, '(i4, 2x, 4g16.8, 4g16.8, 4g16.8, 4g16.8)' ) iflag, t, y(1), y(2), y(3), y(4)

    !do i_step = 1, n_step
    !    t     = ( real(n_step - i_step+1, kind=dp)* t_start &
    !            + real(         i_step-1, kind=dp)* t_stop) &
    !            / real(n_step           , kind=dp)
    !    t_out = ( real(n_step - i_step+1, kind=dp)* t_start &
    !            + real(         i_step-1, kind=dp)* t_stop) &
    !            / real(n_step           , kind=dp)
    !    
    !    call rkf45(func, neqn, y, yp, t, t_out, relerr, abserr, iflag)
    !
    !    write ( 11, '(i4, 2x, 4g14.6, 4g14.6, 4g14.6, 4g14.6)' ) iflag, t, y(1), y(2), y(3), y(4)
    !    if (iflag .ne. 2) then
    !        iflag = 2
    !    end if 
    !    if (y(2) .lt. 0.0) stop 
    !end do
    do i_step = 1, n_step
        t_out = t + dt
        call rkf45(func, neqn, y, yp, t, t_out, relerr, abserr, iflag)
        if (iflag.ne.2) then
            iflag = 2
        end if
        if (y(2).lt.0) stop
        write ( 11, '(i4, 2x, 4g16.8, 4g16.8, 4g16.8, 4g16.8)' ) iflag, t, y(1), y(2), y(3), y(4)
    end do
    close(11)

    return
end subroutine pmwar

subroutine func(t, y, yp)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(kind=dp) :: t, m, d, g
    real(kind=dp) :: y(4), yp(4)

    m = 1.5D0 ! kg
    d = 3.0D0 ! drag coef
    g = 9.81D0! m/s^2

    call r8_fake_use(t)

    ! Differential equations for 
    ! Projectile Motion with Air Resistance
    yp(1) = y(3)            ! v_x = dx/dt = v_0*cos(theta)
    yp(2) = y(4)            ! v_y = dy/dt = v_0*sin(theta)
    yp(3) = -d/m*y(3)       ! a_x = -d/m*v_x
    yp(4) = -d/m*y(4) - g   ! a_y = -d/m*v_y - g
    return
end subroutine func

subroutine r8_fake_use ( x )

    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    real ( kind = dp ) :: x

    if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use: variable is NAN.'
    end if

    return
end subroutine r8_fake_use