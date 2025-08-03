program lotvol
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: neqn = 2
    integer iwork(5)
    real(kind=dp) work(100+21*neqn)
    integer iflag
    integer i_step
    integer nmax
    external f2
    real(kind = dp) abserr
    real(kind = dp) relerr
    real(kind = dp) t, dt
    real(kind = dp) t_out
    real(kind = dp) t_start
    real(kind = dp) t_stop
    real(kind = dp) y(neqn)
    real(kind = dp) yp(neqn)

    write ( *, '(a)' ) ""
    !write ( *, '(a)' ) 'TEST05'
    write ( *, '(a)' ) "  Solve a Lotka-Volterra predator-prey model:"
    write ( *, '(a)' ) " "
    write ( *, '(a)' ) "  y'(1) =  y(1)*(alpha - beta*y(2))"
    write ( *, '(a)' ) "  y'(2) =  y(2)*(y(1)*delta - gamma))"

    abserr = 1.0D-10
    relerr = 1.0D-10

    iflag = 1

    t_start =  0.0D+00
    t_stop  = 15.0D+00
    nmax    = int(t_stop*1.0D+03)

    t     = 0.0D+00
    t_out = 0.0D+00
    dt    = (t_stop - t_start)/real(nmax, kind=dp)

    ! Initial conditions
    y(1) = 10.0D+00  ! initial prey population
    y(2) =  5.0D+00  ! initial predator population

    call func_lv(t, y, yp)
    
    write ( *, '(a)' ) ""
    write ( *, '(a)' ) "Initial conditions:"
    write ( *, '(a, 5F9.1)' ) "Initial prey population:", y(1)
    write ( *, '(a, 5F5.1)' ) "Initial predator population:", y(2)

    !write ( *, '(a)' ) '  FLAG       T              Y(1)            Y(2)'
    
    open(unit=11, file='lv.dat')
    write ( 11, '(i4, 1x, 4E20.10, 4E20.10)' ) iflag, t, y(1), y(2)
    do i_step = 1, nmax
        t_out = t + dt
        call ode  (func_lv, neqn, y, t, t_out, relerr, abserr, iflag, work, iwork )
        !call rkf45(func_lv, neqn, y, yp, t, t_out, relerr, abserr, iflag)

        if (iflag.ne.2) then
            iflag = 2
        end if

        write ( 11, '(i4, 1x, 4E20.10, 4E20.10)' ) iflag, t, y(1), y(2)
    end do
    close(11)

    write ( *, '(a)' ) ""
    write ( *, '(a)' ) "Results exported to output file."

    return
end program

subroutine func_lv(t, y, yp)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(kind = dp) t
    real(kind = dp) y(2)
    real(kind = dp) yp(2)
    real(kind = dp) :: alpha, beta, delta, gamma
    
    alpha = 1.5D0 ! prey population growth parameter
    beta  = 1.0D0 ! predator population extinction parameter
    gamma = 3.0D0 ! species interaction parameter 2
    delta = 1.0D0 ! species interaction parameter 1

    call r8_fake_use ( t )

    ! Lotka-Volterra 1st order ODEs
    yp(1) = y(1)*(alpha - beta*y(2))
    yp(2) = y(2)*(y(1)*delta - gamma) 
    return
end subroutine func_lv
    
subroutine r8_fake_use ( x )
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(kind = dp) x

    if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use: variable is NAN.'
    end if

    return
end subroutine r8_fake_use