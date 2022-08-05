program main
!====================================================================
! Initial value problem for a single 2nd-order ODE
! Call Methods: Runge-Kutta 4th order
! comment: can be used for a system of two first order ODEs
!====================================================================
    implicit none 
    real(8)    :: ti, xi, vi, tf, xf, vf
    real(8)    :: ei, dt
    real(8)    :: tmax
    real(8)    :: t1, t2, time
    real(8), external :: d1x, d2x
    integer(8) :: count

    open (unit = 6, file = "result12.dat")

! initial conditions
    ti =  0.0
    xi = -2.5
    vi =  6.0

! time step and max time
    dt   = 0.01
    tmax = 20.0
    
    ei = vi**2.0/2.0 + xi**2.0/2.0 ! energy for a harmonic oscillator

! print the header and initial conditions
    write (6,*) '       2-nd order single ODE '
    write (6,*) '    Method: Runge-Kutta 4th order'
    write (6,*) '   t       x(t)        v(t)        energy'
    write (6,'(4f10.5)') ti, xi, vi, ei

    count = 0
    call cpu_time(t1)
    do while (ti .le. tmax)
        tf = ti + dt
        call RK4D12(d1x,d2x,ti,xi,vi,tf,xf,vf)
        ei = vf**2/2 + xf**2/2
        write (6,100) tf, xf, vf, ei
        ti = tf
        xi = xf
        vi = vf
        count = count + 1
    end do
    call cpu_time(t2)
    time = t2 - t1
    
    100 format (4f10.5)
    close(6)
    
    open (unit = 7, file = "loops_n_time12.txt")
    write(7,97)
    write(7,98) count
    write(7,99) time
    97 format (/'Initial value problem for a single 2nd-order ODE', /, &
                'Call Methods: Runge-Kutta 4th order.')
    98 format (/'Number of necessary loops:', i6)
    99 format ( 'Calculation time is ', 1pe10.4, ' seconds.'/,/)
    close(7)
! Try plotting with this command in gnuplot
! > plot 'result12.dat' using 0:2 with lines
end program main

function d1x(t,x,v)
!--------------------------------------
! function dx/dt
!--------------------------------------
    implicit none
    real(8) :: d1x, t, x, v

    d1x = v

end function d1x

function d2x(t,x,v)
!--------------------------------------
! function d2x/dt2
!--------------------------------------
    implicit none
    real(8) :: d2x, t, x, v
! simple harmonic oscillator
    real(8), parameter :: k = 1.0, m = 1.0

    d2x = (-1.0*k/m)*x

end function d2x

Subroutine RK4D12(d1x,d2x,ti,xi,vi,tf,xf,vf)
!===========================================================
! Solution of a single second-order ODE d2x/dt2=f(t,x,dx/dt)
! Method:  Runge-Kutta 4th-order
! Alex G. February 2010
!-----------------------------------------------------------
! input ...
! d1x(t,x,v)- function dx/dt   (supplied by a user)
! d2x(t,x,v)- function d2x/dt2 (supplied by a user)
! ti 	- initial time
! xi  - initial position
! vi  - initial velocity (first order derivative)
! tf  - final time (find solution at this time)
! output ...
! xf  - position at point tf
! vf  - velocity at point tf
!===========================================================
    implicit none
    real(8) :: d1x, d2x, ti, xi, vi, tf, xf, vf
    real(8) :: h, t, k1x, k2x, k3x, k4x, k1v, k2v, k3v, k4v
    
    h = tf - ti
    t = ti
    
    k1x = h*d1x(t,xi,vi)
    k1v = h*d2x(t,xi,vi)
    
    k2x = h*d1x(t+h/2.0,xi+k1x/2.0,vi+k1v/2.0)
    k2v = h*d2x(t+h/2.0,xi+k1x/2.0,vi+k1v/2.0)
    
    k3x = h*d1x(t+h/2.0,xi+k2x/2.0,vi+k2v/2.0)
    k3v = h*d2x(t+h/2.0,xi+k2x/2.0,vi+k2v/2.0)
    
    k4x = h*d1x(t+h,xi+k3x,vi+k3v)
    k4v = h*d2x(t+h,xi+k3x,vi+k3v)
    
    xf = xi + (k1x + 2.0*(k2x+k3x) + k4x)/6.0
    vf = vi + (k1v + 2.0*(k2v+k3v) + k4v)/6.0
    
end subroutine RK4D12