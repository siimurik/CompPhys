!   gfortran vanderpol_main.f90 SLATPACK.f90 -o vander

program vanderpol_solver
    implicit none
    
    ! Problem parameters
    integer, parameter :: neq = 2
    integer, parameter :: lrw = 130 + 21*neq    ! Required work array size
    integer, parameter :: liw = 51              ! Required integer work array size
    integer, parameter :: npts = 1000           ! Number of output points
    
    ! Variables for DDEABM
    double precision :: t, tout, y(neq)
    double precision :: rtol, atol
    double precision :: rwork(lrw), rpar(1)
    double precision :: energy
    integer :: info(15), iwork(liw), ipar(1)
    integer :: idid, i
    
    ! Problem-specific parameters
    double precision :: mu                      ! Van der Pol parameter
    double precision :: tstart, tend, dt
    
    ! Output file
    integer, parameter :: unit_out = 10

    external vanderpol_deriv

    
    write(*,*) '================================================='
    write(*,*) 'Van der Pol Oscillator using DDEABM Solver'
    write(*,*) '================================================='
    write(*,*) 'Equation: x'' - mu*(1-x²)*x'' + x = 0'
    write(*,*) 'System:   dx/dt = y'
    write(*,*) '          dy/dt = mu*(1-x²)*y - x'
    write(*,*)
    
    ! Get user input
    write(*,*) 'Enter Van der Pol parameter mu (try 0.5, 2.0, or 10.0): '
    read(*,*) mu
    write(*,*) 'Enter initial conditions:'
    write(*,*) 'x(0) = '
    read(*,*) y(1)
    write(*,*) 'y(0) = x''(0) = '
    read(*,*) y(2)
    write(*,*) 'Enter time span:'
    write(*,*) 'Start time: '
    read(*,*) tstart
    write(*,*) 'End time: '
    read(*,*) tend
    
    ! Store mu in parameter array for passing to derivative function
    rpar(1) = mu
    
    ! Initialize problem
    t = tstart
    dt = (tend - tstart) / dble(npts-1)
    
    ! Set tolerances
    rtol = 1.0d-8
    atol = 1.0d-10
    
    ! Initialize INFO array (all zeros for simplest use)
    do i = 1, 15
        info(i) = 0
    end do
    info(1) = 0  ! First call
    
    ! Initialize work arrays
    do i = 1, lrw
        rwork(i) = 0.0d0
    end do
    do i = 1, liw
        iwork(i) = 0
    end do
    
    ! Open output file
    open(unit_out, file='vanderpol_solution.dat', status='replace')
    write(unit_out,'(A)') '# Van der Pol Oscillator Solution'
    write(unit_out,'(A,F8.3)') '# mu = ', mu
    write(unit_out,'(A)') '# t, x, y, energy'
    
    write(*,*)
    write(*,*) 'Solving Van der Pol equation...'
    write(*,*) 'mu =', mu
    write(*,*) 'Initial conditions: x(0) =', y(1), ', y(0) =', y(2)
    write(*,*) 'Time span: [', tstart, ',', tend, ']'
    write(*,*)
    
    ! Output initial condition
    write(unit_out,'(4ES16.8)') t, y(1), y(2), 0.5d0*(y(1)**2 + y(2)**2)
    
    ! Integration loop
    do i = 1, npts-1
        tout = tstart + i * dt
        
        ! Call DDEABM solver
        call ddeabm(vanderpol_deriv, neq, t, y, tout, info, rtol, atol, &
                   idid, rwork, lrw, iwork, liw, rpar, ipar)
        
        ! Check for successful integration
        if (idid < 0) then
            write(*,*) 'Integration failed with IDID =', idid
            select case(idid)
                case(-1)
                    write(*,*) 'Too many steps attempted (500+)'
                case(-2)
                    write(*,*) 'Error tolerances too stringent'
                case(-3)
                    write(*,*) 'Local error test failed'
                case(-4)
                    write(*,*) 'Problem appears to be stiff'
                    write(*,*) 'Consider using DDEBDF instead'
                case(-33)
                    write(*,*) 'Invalid input detected'
                    exit
                case default
                    write(*,*) 'Unknown error'
            end select
            
            if (idid == -4) then
                write(*,*) 'Continuing with DDEABM (inefficient for stiff problems)...'
                info(1) = 1  ! Continue integration
                cycle
            else if (idid == -1 .or. idid == -2) then
                info(1) = 1  ! Continue integration
                cycle
            else
                exit
            end if
        end if
        
        ! Calculate energy-like quantity for analysis
        ! For Van der Pol, energy is not conserved, but we can track it
        energy = 0.5d0 * (y(1)**2 + y(2)**2)
        
        ! Output solution
        write(unit_out,'(4ES16.8)') t, y(1), y(2), energy
        
        ! Print progress
        if (mod(i, npts/10) == 0) then
            write(*,'(A,F8.3,A,F8.3,A,F8.3,A,F8.3)') &
                'Progress: t=', t, ', x=', y(1), ', y=', y(2), ', E=', energy
        end if
        
        ! Reset INFO(1) for continuation
        info(1) = 1
    end do
    
    close(unit_out)
    
    write(*,*)
    write(*,*) 'Integration completed successfully!'
    write(*,*) 'Final time: t =', t
    write(*,*) 'Final solution: x =', y(1), ', y =', y(2)
    write(*,*) 'Solution saved to: vanderpol_solution.dat'
    write(*,*)
    write(*,*) 'Integration statistics:'
    write(*,*) 'Steps taken:', iwork(26)
    write(*,*) 'Function evaluations:', iwork(27)
    write(*,*) 'Last step size:', rwork(11)
    
    
end program vanderpol_solver

!===============================================================================
! Van der Pol derivative function
!===============================================================================
subroutine vanderpol_deriv(t, y, yprime, rpar, ipar)
    implicit none
    double precision, intent(in) :: t
    double precision, intent(in) :: y(2)
    double precision, intent(out) :: yprime(2)
    double precision, intent(in) :: rpar(*)
    integer, intent(in) :: ipar(*)
    
    double precision :: mu, x, xdot
    
    ! Get parameters
    mu = rpar(1)
    x = y(1)
    xdot = y(2)
    
    ! Van der Pol equations:
    ! dx/dt = y
    ! dy/dt = mu*(1-x²)*y - x
    yprime(1) = xdot
    yprime(2) = mu * (1.0d0 - x*x) * xdot - x
    
end subroutine vanderpol_deriv