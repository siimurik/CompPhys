! Compile and execute with:
!   $ gfortran -O2 lotka_volterra.f90 rkf45.f90 -o lv
!   $ ./lv
program main

!*****************************************************************************80
!
!! RKF45_TEST tests the RKF45 library.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 June 2006
!
!  Author:
!
!    John Burkardt
!
    implicit none

    integer, parameter :: rk = kind ( 1.0D+00 )

    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RKF45_TEST'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) '  Test the RKF45 library.'

!  call test04 ( )
    call test05 ( )
!  call test06 ( )
!
!  Terminate.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RKF45_TEST'
    write ( *, '(a)' ) '  Normal end of execution.'
    write ( *, '(a)' ) ' '
    call timestamp ( )

    stop 0
end

subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 solves a vector ODE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 June 2006
!
!  Author:
!
!    John Burkardt
!
    implicit none

    integer, parameter :: rk = kind ( 1.0D+00 )

    integer, parameter :: neqn = 2

    real ( kind = rk ) abserr
    integer flag
    integer i_step
    integer n_step
    external f2
    real ( kind = rk ) relerr
    real ( kind = rk ) t
    real ( kind = rk ) t_out
    real ( kind = rk ) t_start
    real ( kind = rk ) t_stop
    real ( kind = rk ) y(neqn)
    real ( kind = rk ) yp(neqn)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05'
    write ( *, '(a)' ) '  Solve a vector equation using rkf45:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Y''(1) =  Y(2)'
    write ( *, '(a)' ) '  Y''(2) = -Y(1)'

    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    flag = 1

    t_start = 0.0D+00
    !t_stop = 2.0D+00 * 3.14159265D+00
    t_stop = 2000.0D0

    n_step = 10000

    t     = 0.0D+00
    t_out = 0.0D+00

    ! Initial conditions
    y(1) = 6.0D+00  ! initial prey population
    y(2) = 8.0D+00  ! initial predator population

    call f2 ( t, y, yp )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  FLAG       T          Y(1)          Y(2)'
    write ( *, '(a)' ) ' '
    write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

    do i_step = 1, 12
        t = ( real ( n_step - i_step + 1, kind = rk ) * t_start &
            + real (          i_step - 1, kind = rk ) * t_stop ) &
            / real ( n_step,              kind = rk )

        t_out = ( real ( n_step - i_step, kind = rk ) * t_start &
                + real (          i_step, kind = rk ) * t_stop ) &
                / real ( n_step,          kind = rk )

        call rkf45 ( f2, neqn, y, yp, t, t_out, relerr, abserr, flag )

        write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)
    end do

    open(unit=11, file='lv4.dat')
    do i_step = 1, n_step
        t = ( real ( n_step - i_step + 1, kind = rk ) * t_start &
            + real (          i_step - 1, kind = rk ) * t_stop ) &
            / real ( n_step,              kind = rk )

        t_out = ( real ( n_step - i_step, kind = rk ) * t_start &
                + real (          i_step, kind = rk ) * t_stop ) &
                / real ( n_step,          kind = rk )

        call rkf45 ( f2, neqn, y, yp, t, t_out, relerr, abserr, flag )

        write ( 11, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)
    end do
    close(11)

    return
end


subroutine f2 ( t, y, yp )
!*****************************************************************************80
!
!! F2 evaluates the derivative for the ODE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T, the value of the independent variable.
!
!    Input, real ( kind = rk ) Y(NEQN), the value of the dependent variable.
!
!    Output, real ( kind = rk ) YP(NEQN), the value of the derivative
!    dY(1:NEQN)/dT.
!
    implicit none

    integer, parameter :: rk = kind ( 1.0D+00 )

    real ( kind = rk ) t
    real ( kind = rk ) y(2)
    real ( kind = rk ) yp(2)
    real ( kind = rk ) :: alpha, beta, delta, gamma
    
    alpha = 2.0D0 ! prey population growth parameter
    beta  = 1.0D0 ! predator population extinction parameter
    delta = 1.0D0 ! species interaction parameter 1
    gamma = 3.0D0 ! species interaction parameter 2


    call r8_fake_use ( t )

    yp(1) = y(1)*(alpha - delta*y(2))
    yp(2) = y(2)*(y(1)*gamma - beta) 
    return
end

subroutine r8_fake_use ( x )
!*****************************************************************************80
!
!! r8_fake_use() pretends to use an R8 variable.
!
!  Discussion:
!
!    Some compilers will issue a warning if a variable is unused.
!    Sometimes there's a good reason to include a variable in a program,
!    but not to use it.  Calling this function with that variable as
!    the argument will shut the compiler up.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the variable to be "used".
!
    implicit none

    integer, parameter :: rk = kind ( 1.0D+00 )

    real ( kind = rk ) x

    if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use: variable is NAN.'
    end if

    return
end