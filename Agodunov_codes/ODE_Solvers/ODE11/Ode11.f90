program main
!====================================================================
! Initial value problem for a single 1st-order ODE
! Call Methods: Euler, Predictor-Corrector, Runge-Kutta 4th order
!====================================================================
    implicit none
    real(8)             :: dt, xi, ti, xf, tf, tmax
    real(8)             :: eps, t, told, ft, ftold, tnew
    real(8)             :: t1, t2, time
    integer(8)          :: key, count
    real(8), external   :: dx

! initial conditions
    xi = 1.00
    ti = 0.00

! "time" step and max time
    dt  = 0.01
! If you are using a function with a non-analytical solution
! use the NewtonRaphson(...) subroutine. Also determine the 
! desired accuracy and input your function into 'dx(t,x)'
! as well as into the subroutine down below.
    eps = 1E-8 
!    call NewtonRaphson(eps,t,told,ft,ftold,tnew)
!    tmax = t
! For the first function
    tmax = 20.0
! select a method:
! key = 1 Euler, key = 2 predictor-corrector, key = 3 RK 4th order
! Let user choose the method for solving
    print *, 'Input the solver method:'
    print *, '[1] Euler method'
    print *, '[2] Predictor-Corrector method'
    print *, '[3] Runge-Kutta 4th order'
    read *, key

! open files (if needed)


! print the header and initial conditions  
    if (key==1) then
    write (*,*) ('   1st order single ODE will be solved by:')
    write (*,*)  '   Simple Euler method'
    open (unit=6, file = 'result1_euler.dat')
    end if
    if (key==2) then
    write (*,*) ('   1st order single ODE will be solved by:')
    write (*,*)  '   Predictor-Corrector'
    open (unit=6, file = 'result1_pred_cor.dat')
    end if
    if (key==3) then
    write (*,*) ('   1st order single ODE will be solved by:')
    write (*,*)  '   Runge-Kutta 4th order'
    open (unit=6, file = 'result1_rk.dat')
    end if
    write (*,*) ''
    write (*,*) '   t       x(t)    '
    write(6,100) ti, xi

    count = 0
    call cpu_time(t1)
    do while (ti <= tmax)
        tf = ti + dt
        if (key==1) call euler1 (dx,ti,xi,tf,xf)
        if (key==2) call euler1m(dx,ti,xi,tf,xf)
        if (key==3) call RK4D11 (dx,ti,xi,tf,xf)
        write (6,100) tf, xf
        ti = tf
        xi = xf
        count = count + 1
    end do
    call cpu_time(t2)
    time = t2 - t1

    100 format (2f10.5)
    close(6)

    open (unit = 7, file = "loops_n_time1.txt")
    write(7,97)
    write(7,98) count
    write(7,99) time
    97 format (/'Initial value problem for a single 2nd-order ODE', /, &
                'Call Methods: Runge-Kutta 4th order.')
    98 format (/'Number of necessary loops:', i8)
    99 format ( 'Calculation time is ', 1pe10.4, ' seconds.'/,/)
    close(7)
end program main

function dx(t,x)
!----------------------------------------------
! dx(t,x)- function dx/dt in the 1st order ODE
!----------------------------------------------
    implicit none
    real(8) :: dx, x, t

    dx = (-1.0)*x   ! solution x = exp(-t)
! the carrying capacity problem
!   dx = 0.1*x - 2.0E-4*x**2

end function dx

subroutine RK4D11(dx,ti,xi,tf,xf)
!==================================================
! Solution of a single 1st order ODE dx/dt=f(x,t)
! Method: 4th-order Runge-Kutta method
! Alex G. February 2010
!--------------------------------------------------
! input ...
! dx(t,x)- function dx/dt (supplied by a user)
! ti  - initial time
! xi  - initial position
! tf  - time for a solution
! output ...
! xf  - solution at point tf
!==================================================
    implicit none
    real(8) :: dx, ti, xi, tf, xf
    real(8) :: h, k1, k2, k3, k4

    h  = tf - ti

    k1 = h*dx(ti      , xi       )
    k2 = h*dx(ti+h/2.0, xi+k1/2.0)
    k3 = h*dx(ti+h/2.0, xi+k2/2.0)
    k4 = h*dx(ti+h    , xi+k3    )

    xf = xi + (k1 + 2.0*(k2+k3) + k4)/6.0

end subroutine RK4D11

subroutine euler1m(dx,ti,xi,tf,xf)
!==================================================
! Solution of a single 1st order ODE dx/dt=f(x,t)
! Method:  Modified Euler (predictor-corrector)
! Alex G. February 2010
!--------------------------------------------------
! input ...
! dx(t,x)- function dx/dt (supplied by a user)
! ti - initial time
! xi  - initial position
! tf  - time for a solution
! output ...
! xf  - solution at point tf
!==================================================
    implicit none
    real(8) :: dx, xi, ti, xf, tf
    
    xf = xi + dx(ti,xi)*(tf-ti)
    xf = xi + (dx(ti,xi) + dx(tf,xf))*0.5*(tf-ti)
    
end subroutine euler1m
    
subroutine euler1(dx,ti,xi,tf,xf)
!==================================================
! Solution of a single 1st order ODE dx/dt=f(x,t)
! Method:  Simple Euler (predictor-corrector)
! Alex G. February 2010
!--------------------------------------------------
! input ...
! dx(t,x)- function dx/dt (supplied by a user)
! ti - initial time
! xi  - initial position
! tf  - time for a solution
! output ...
! xf  - solution at point tf
!==================================================
    implicit none
    real(8) :: dx, xi, ti, xf, tf
    
    xf = xi + dx(ti,xi)*(tf-ti)

end subroutine euler1

! A specific solution for the problem without 
! the nonanalytical solution. Specific for this prob.
subroutine NewtonRaphson(eps, t, told, ft, ftold, tnew)
!==================================================
! The Newton-Raphson method
! Author: Siim Erik Pugal, July 2022
!--------------------------------------------------
! input ...
! eps   - desired accuracy for calculations
! told  - a random real number (NB! told .ne. t)
! tnew  - next iterative step
!
! NB! must be modified specifically ...
! f(t)    - function's output with the new output
! f(told) - function's output with the old output
!
! output ...
! t     - final time instance
!==================================================
    real(8) :: eps, t, told, ft, ftold, tnew
    integer(8) :: count

! Rough position where the function crosses 
! the horizontal axis the second time.
    t = 490     ! seconds
    told = t - 1

    count   = 0
    do while (abs(t-told) >= eps)
        count   = count + 1
        ft      = 0.1*t    - 2.0E-4*t**2
        ftold   = 0.1*told - 2.0E-4*told**2
        tnew    = t - ft*(t-told)/(ft-ftold)  ! The modified Newton–Raphson method
        told    = t
        t       = tnew
    end do
    
end subroutine NewtonRaphson