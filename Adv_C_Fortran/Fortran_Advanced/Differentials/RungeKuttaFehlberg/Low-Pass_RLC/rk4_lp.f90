PROGRAM low_pass
    ! Circuit equations:
    ! V = i_R * R
    ! i_C = C * dV/dt  
    ! L * di_L/dt = V_b - V
    ! i_L = i_R + i_C
    !
    ! System of ODEs:
    ! dV/dt = (1/C) * (i_L - V/R)
    ! di_L/dt = (1/L) * (V_b - V)

    IMPLICIT NONE
    INTEGER, PARAMETER :: neq=2!, NSTPMX=2000
    INTEGER :: nvar, nstep, k
    REAL :: tstart, tstop, dt
    REAL :: dfdt(neq)  ! Initial conditions [position, velocity]
    REAL, ALLOCATABLE, DIMENSION(:) :: tout!, y(NMAX,NSTPMX)
    REAL, ALLOCATABLE, DIMENSION(:) :: V, i_R, i_L, i_C!, y(NMAX,NSTPMX)
    REAL, ALLOCATABLE, DIMENSION(:,:) :: y
    REAL ::  Vb, C, R, L
    common /circuit/ Vb, C, R, L  ! Share values with derivs() function.
    EXTERNAL derivs
    !COMMON /path/ xx, y ! Storage of results as required by rkdumb

    nvar = neq      ! Number of variables (position and velocity)
    tstart = 0.0    ! Start time
    tstop = 2.0     ! End time 
    dt = 0.001      ! Time-step size
    nstep = int((tstop - tstart) / dt) + 1  ! Total number of steps

    ! Circuit parameters
    Vb = 24.0
    C = 0.001
    R = 100.0
    L = 1.0
    ! P.S. These must be set before rkdumb() can run

    ! +1 necessary due to how rkdumb() works
    allocate(tout(nstep + 1))
    allocate(y(nvar, nstep + 1))

    allocate(V(nstep + 1))
    allocate(i_R(nstep + 1))
    allocate(i_L(nstep + 1))
    allocate(i_C(nstep + 1))

    ! Initial conditions: 
    dfdt(1) = 0.0  ! Voltage
    dfdt(2) = 0.0  ! Current in inductor

    WRITE(*,'(A, F6.3, A)') 'Initial voltage =', dfdt(1), ' V'
    WRITE(*,'(A, F6.3, A)') 'Initial current =', dfdt(2), ' A'
    WRITE(*,'(A, ES12.4)') 'Step size dt:', dt
    WRITE(*,'(A, i5)') 'Number of steps =', nstep
    WRITE(*,*)

    ! Call the integration routine
    CALL rkdumb(dfdt, neq, tstart, tstop, nstep, derivs, tout, y)

    ! Calculate missing values
    do k = 1, nstep+1
        V = y(1,k)         ! Determined by rk4
        i_L(k) = y(2,k)    ! Determined by rk4

        i_R(k) = V(k) / R           ! From the original circuit
        i_C(k) = i_L(k) - i_R(k)    ! From the original circuit
    end do

    ! Print results
    WRITE(*,*) 'Results:'
    WRITE(*,*) '  Time        V           i_L         i_R        i_C'
    WRITE(*,*) '-------------------------------------------------------'
    DO k = 1, nstep+1, 200  ! Print every 200th point
        WRITE(*,'(F8.4, 2X, F10.6, 2X, F10.6, 2X, F10.6,2X, F10.6)') &
            tout(k), y(1,k), y(2,k), i_R(k), i_C(k)
    END DO

    ! output data into a file 
    open(1, file = 'data1.dat', status = 'replace')  
    do k=1, nstep+1
        WRITE(1,'(F8.4, 2X, F10.6, 2X, F10.6)') tout(k), y(1,k), y(2,k)
    end do  
    close(1)

    deallocate(tout, y, V, i_R, i_L, i_C)

END PROGRAM low_pass

!-----------------------------------------------------------------------
SUBROUTINE derivs(x, y, dydx, neq)
    ! Subroutine to calculate derivatives for the harmonic oscillator
    ! y(1) = position (x)
    ! y(2) = velocity (dx/dt)
    ! dydx(1) = dy_vec(1)/dt = velocity = y(2)
    ! dydx(2) = dy_vec(2)/dt = acceleration = -ω²*y(1)

    IMPLICIT NONE
    INTEGER :: neq
    REAL :: x, y(neq), dydx(neq)
    REAL :: Vb, V, iL, C, R, L
    common /circuit/ Vb, C, R, L

    if (.false.) then
        write (*,*) x
    end if

    V = y(1)      ! Voltage across capacitor
    iL = y(2)     ! Current through inductor

    dydx(1) = (1.0 / C) * (iL - V / R)
    dydx(2) = (1.0 / L) * (Vb - V)

END SUBROUTINE derivs

!-----------------------------------------------------------------------
SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs)
    INTEGER n
    !INTEGER :: NMAX=neq ! Set to the maximum number of functions.
    REAL h,x,dydx(n),y(n),yout(n)
    EXTERNAL derivs
    ! Given values for the variables y(1:n) and their derivatives dydx(1:n) known at x, use
    ! the fourth-order Runge-Kutta method to advance the solution over an interval h and return
    ! the incremented variables as yout(1:n), which need not be a distinct array from y. The
    ! user supplies the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
    INTEGER i
    REAL h6,hh,xh,dym(n),dyt(n),yt(n)
    hh=h*0.5
    h6=h/6.
    xh=x+hh
    do i=1,n ! First step.
        yt(i)=y(i)+hh*dydx(i)
    end do
    call derivs(xh,yt,dyt,n) ! Second step.
    do i=1,n
        yt(i)=y(i)+hh*dyt(i)
    end do
    call derivs(xh,yt,dym,n) !Third step.
    do i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
    end do
    call derivs(x+h,yt,dyt,n) ! Fourth step.
    do i=1,n ! Accumulate increments with proper weights.
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
    end do
    return
END SUBROUTINE

!-----------------------------------------------------------------------
SUBROUTINE rkdumb(dfdt,neq,tstart,tstop,nstep,derivs,xx,y)
    IMPLICIT NONE
    INTEGER nstep
    integer neq
    !INTEGER, PARAMETER :: NMAX=50, NSTPMX=2000 ! Maximum number of functions and
    !maximum number of values to be stored.
    REAL tstart,tstop,dfdt(neq),xx(nstep),y(neq,nstep)
    EXTERNAL derivs
    !COMMON /path/ xx, y ! Storage of results.
    ! C USES rk4
    ! Starting from initial values dfdt(1:neq) known at tstart use fourth-order Runge-Kutta to
    ! advance nstep equal increments to tstop. The user-supplied subroutine derivs(x,v,dvdx)
    ! evaluates derivatives. Results are stored in the common block path. Be sure to dimension
    ! the common block appropriately.
    INTEGER i,k
    REAL h,x,dv(neq),v(neq)
    do i=1,neq !Load starting values.
        v(i)=dfdt(i)
        y(i,1)=v(i)
    end do
    xx(1)=tstart
    x=tstart
    h=(tstop-tstart)/nstep
    do k=1,nstep ! Take nstep steps.
        call derivs(x,v,dv,neq)
        call rk4(v,dv,neq,x,h,v,derivs)
        if(x+h.eq.x) stop 'stepsize not significant in rkdumb'
        x=x+h
        xx(k+1)=x ! Store intermediate steps.
        do i=1,neq
            y(i,k+1)=v(i)
        end do
    end do
    return
END SUBROUTINE