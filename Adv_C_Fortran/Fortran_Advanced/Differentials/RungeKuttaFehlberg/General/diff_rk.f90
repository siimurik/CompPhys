PROGRAM harmonic_oscillator
    ! Example program using rk4 and rkdumb subroutines
    ! Solves the simple harmonic oscillator: d²x/dt² + ω²x = 0
    ! This is converted to a system of first-order ODEs:
    ! dy1/dt = y2 (where y1 = x, y2 = dx/dt)
    ! dy2/dt = -ω²*y1

    IMPLICIT NONE
    INTEGER, PARAMETER :: NMAX=50, NSTPMX=200
    INTEGER :: nvar, nstep, k
    REAL :: x1, x2, omega, period
    REAL :: vstart(2)  ! Initial conditions [position, velocity]
    REAL :: xx(NSTPMX), y(NMAX,NSTPMX)
    EXTERNAL derivs
    COMMON /path/ xx, y  ! Storage of results as required by rkdumb
    COMMON /params/ omega  ! Parameters for the differential equation

    ! Set up the problem parameters
    omega = 2.0  ! Angular frequency
    period = 2.0 * 3.14159 / omega  ! One complete period
    nvar = 2     ! Number of variables (position and velocity)
    nstep = 100  ! Number of integration steps
    x1 = 0.0     ! Start time
    x2 = period  ! End time (one complete period)
    ! Initial conditions: start at x=1.0, v=0.0 (classic SHO)
    vstart(1) = 1.0  ! Initial position
    vstart(2) = 0.0  ! Initial velocity

    WRITE(*,*) 'Simple Harmonic Oscillator Solution'
    WRITE(*,*) 'Omega =', omega
    WRITE(*,*) 'Period =', period
    WRITE(*,*) 'Initial position =', vstart(1)
    WRITE(*,*) 'Initial velocity =', vstart(2)
    WRITE(*,*) 'Number of steps =', nstep
    WRITE(*,*)

    ! Call the integration routine
    CALL rkdumb(vstart, nvar, x1, x2, nstep, derivs)

    ! Print results
    WRITE(*,*) 'Results:'
    WRITE(*,*) 'Time      Position    Velocity'
    WRITE(*,*) '--------------------------------'
    DO k = 1, nstep+1, 10  ! Print every 10th point
        WRITE(*,'(F8.4, 2X, F10.6, 2X, F10.6)') xx(k), y(1,k), y(2,k)
    END DO

    ! Verify conservation of energy
    WRITE(*,*)
    WRITE(*,*) 'Energy Conservation Check:'
    WRITE(*,*) 'Initial energy =', 0.5 * omega**2 * vstart(1)**2 + 0.5 * vstart(2)**2
    WRITE(*,*) 'Final energy   =', 0.5 * omega**2 * y(1,nstep+1)**2 + 0.5 * y(2,nstep+1)**2

END PROGRAM harmonic_oscillator

!-----------------------------------------------------------------------
SUBROUTINE derivs(x, y, dydx)
    ! Subroutine to calculate derivatives for the harmonic oscillator
    ! y(1) = position (x)
    ! y(2) = velocity (dx/dt)
    ! dydx(1) = dy_vec(1)/dt = velocity = y(2)
    ! dydx(2) = dy_vec(2)/dt = acceleration = -ω²*y(1)

    IMPLICIT NONE
    REAL :: x, y(2), dydx(2)
    REAL :: omega
    COMMON /params/ omega

    dydx(1) = y(2)              ! dx/dt = velocity
    dydx(2) = -omega**2 * y(1)  ! dv/dt = -ω²*x

END SUBROUTINE derivs

!-----------------------------------------------------------------------
SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs)
    INTEGER, PARAMETER :: NMAX=50 ! Set to the maximum number of functions.
    INTEGER n
    REAL h,x,dydx(n),y(n),yout(n)
    EXTERNAL derivs
    ! Given values for the variables y(1:n) and their derivatives dydx(1:n) known at x, use
    ! the fourth-order Runge-Kutta method to advance the solution over an interval h and return
    ! the incremented variables as yout(1:n), which need not be a distinct array from y. The
    ! user supplies the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
    INTEGER i
    REAL h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)
    hh=h*0.5
    h6=h/6.
    xh=x+hh
    do i=1,n ! First step.
        yt(i)=y(i)+hh*dydx(i)
    end do
    call derivs(xh,yt,dyt) ! Second step.
    do i=1,n
        yt(i)=y(i)+hh*dyt(i)
    end do
    call derivs(xh,yt,dym) !Third step.
    do i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
    end do
    call derivs(x+h,yt,dyt) ! Fourth step.
    do i=1,n ! Accumulate increments with proper weights.
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
    end do
    return
END SUBROUTINE

!-----------------------------------------------------------------------
SUBROUTINE rkdumb(vstart,nvar,x1,x2,nstep,derivs)
    INTEGER nstep,nvar
    INTEGER, PARAMETER :: NMAX=50, NSTPMX=200 ! Maximum number of functions and
    !maximum number of values to be stored.
    REAL x1,x2,vstart(nvar),xx(NSTPMX),y(NMAX,NSTPMX)
    EXTERNAL derivs
    COMMON /path/ xx, y ! Storage of results.
    ! C USES rk4
    ! Starting from initial values vstart(1:nvar) known at x1 use fourth-order Runge-Kutta to
    ! advance nstep equal increments to x2. The user-supplied subroutine derivs(x,v,dvdx)
    ! evaluates derivatives. Results are stored in the common block path. Be sure to dimension
    ! the common block appropriately.
    INTEGER i,k
    REAL h,x,dv(NMAX),v(NMAX)
    do i=1,nvar !Load starting values.
        v(i)=vstart(i)
        y(i,1)=v(i)
    end do
    xx(1)=x1
    x=x1
    h=(x2-x1)/nstep
    do k=1,nstep ! Take nstep steps.
        call derivs(x,v,dv)
        call rk4(v,dv,nvar,x,h,v,derivs)
        if(x+h.eq.x) stop 'stepsize not significant in rkdumb'
        x=x+h
        xx(k+1)=x ! Store intermediate steps.
        do i=1,nvar
            y(i,k+1)=v(i)
        end do
    end do
    return
END SUBROUTINE