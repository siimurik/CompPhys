PROGRAM ODE_Solver
    INTEGER, PARAMETER :: N = 2  ! Number of variables
    INTEGER, PARAMETER :: NSTEP = 100  ! Number of steps
    REAL :: x1, x2, vstart(N)
    REAL :: xx(NSTEP), y(N, NSTEP)
    integer :: i
    COMMON /path/ xx, y !Storage of results.
    
    ! Set initial values
    vstart(1) = 1.0
    vstart(2) = 0.0
    
    ! Set initial and final values of x
    x1 = 0.0
    x2 = 10.0
    
    ! Call the derivative evaluation subroutine
    CALL derivs(x, v, dv)
    
    ! Call the rkdumb subroutine to solve the ODEs
    CALL rkdumb(vstart, N, x1, x2, NSTEP, derivs)
    
    ! Retrieve results from the common block path
    print *, y
    
END PROGRAM ODE_Solver


SUBROUTINE derivs(x, y, dydx)
    REAL, INTENT(IN) :: x
    REAL, INTENT(IN) :: y(:)
    REAL, INTENT(OUT) :: dydx(size(y))
    
    dydx(1) = -0.5*y(1) + 0.1*y(2) ! Derivative of y1
    dydx(2) = 0.5*y(1) - 0.1*y(2)  ! Derivative of y2

END SUBROUTINE derivs


SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs)
    INTEGER, PARAMETER :: NMAX=50 !Set to the maximum number of functions.
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
    do i=1,n             ! First step.
        yt(i)=y(i)+hh*dydx(i)
    end do
    
    call derivs(xh,yt,dyt)  ! Second step.
    
    do i=1,n
        yt(i)=y(i)+hh*dyt(i)
    end do
    
    call derivs(xh,yt,dym)  !Third step.
    
    do i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
    end do
    
    call derivs(x+h,yt,dyt) ! Fourth step.
    
    do i=1,n     ! Accumulate increments with proper weights.
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
    end do
    return
END

SUBROUTINE rkdumb(vstart,nvar,x1,x2,nstep,derivs)
    INTEGER nstep,nvar
    INTEGER, PARAMETER :: NMAX=50, NSTPMX=200 ! Maximum number of functions and
    !maximum number of values to    be stored.
    REAL x1,x2,vstart(nvar),xx(NSTPMX),y(NMAX,NSTPMX)
    EXTERNAL derivs
    COMMON /path/ xx, y !Storage of results.
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
END