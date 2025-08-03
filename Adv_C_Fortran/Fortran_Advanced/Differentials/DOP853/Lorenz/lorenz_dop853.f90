! Compile and execute with:
! $ gfortran lorenz_dop853.f90 dop853.f90 -o lorenz_dop853
! $ ./lorenz_dop853

! DOP853 solver for chaotic dynamical systems
! Solves the Lorenz system - a classic chaotic benchmark:
! dx/dt = σ(y - x)
! dy/dt = x(ρ - z) - y  
! dz/dt = xy - βz
! Classic parameters: σ=10, ρ=28, β=8/3

program dop853_lorenz
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    external florenz, solout_lorenz, f
    
    ! Problem parameters
    integer, parameter :: neq = 3
    real(kind=dp), dimension(neq) :: y
    real(kind=dp) :: x, xend, rtol, atol, rpar
    integer :: n, iout, itol, ipar, idid
    integer :: lwork, liwork
    real(kind=dp), allocatable, dimension(:) :: work
    integer, allocatable, dimension(:) :: iwork
    
    ! Calculate required array sizes and allocate
    call calculate_dop853_arrays(neq, lwork, liwork)
    allocate(work(lwork))
    allocate(iwork(liwork))
    
    ! Initialize problem parameters
    call setup_lorenz_problem(neq, y, x, xend, rtol, atol, rpar, ipar, &
                             n, iout, itol, work, iwork, lwork, liwork)
    
    ! Run the integration
    call run_dop853_integration(florenz, solout_lorenz, n, y, x, xend, &
                               rtol, atol, itol, iout, work, lwork, &
                               iwork, liwork, rpar, ipar, idid)
    
    ! Print final statistics
    call print_dop853_statistics(n, x, y, rtol, liwork, iwork, idid)
    
    ! Clean up
    deallocate(work)
    deallocate(iwork)
    
end program dop853_lorenz

! Subroutine to calculate DOP853 array lengths
subroutine calculate_dop853_arrays(neq, lwork, liwork)
    implicit none
    integer, intent(in) :: neq
    integer, intent(out) :: lwork, liwork
    
    ! For dense output with all components
    lwork = 11 * neq + 8 * neq + 21
    liwork = neq + 21
    
    write(6, '(A,I0,A,I0)') 'Calculated DOP853 array sizes: LWORK = ', lwork, ', LIWORK = ', liwork
    
end subroutine calculate_dop853_arrays

! Setup the Lorenz chaotic system problem
subroutine setup_lorenz_problem(neq, y, x, xend, rtol, atol, rpar, ipar, &
                               n, iout, itol, work, iwork, lwork, liwork)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: neq, lwork, liwork
    real(kind=dp), dimension(neq), intent(out) :: y
    real(kind=dp), intent(out) :: x, xend, rtol, atol, rpar
    integer, intent(out) :: n, iout, itol, ipar
    real(kind=dp), dimension(lwork), intent(inout) :: work
    integer, dimension(liwork), intent(inout) :: iwork
    integer :: i
    
    write(6, '(A)') 'Setting up Lorenz chaotic system...'
    write(6, '(A)') 'Equations:'
    write(6, '(A)') '  dx/dt = σ(y - x)'
    write(6, '(A)') '  dy/dt = x(ρ - z) - y'
    write(6, '(A)') '  dz/dt = xy - βz'
    write(6, '(A)') 'Parameters: σ = 10, ρ = 28, β = 8/3 (classic chaotic regime)'
    write(6, '(A)') ''
    
    ! Problem dimension
    n = neq
    
    ! DOP853 options
    iout = 3      ! Dense output for all components
    itol = 0      ! Scalar tolerances
    
    ! Initial conditions: slightly perturbed from origin
    x = 0.0_dp               ! Initial time
    y(1) = 1.0_dp            ! x = 1.0 (slight perturbation)
    y(2) = 1.0_dp            ! y = 1.0 (slight perturbation)
    y(3) = 1.0_dp            ! z = 1.0 (slight perturbation)
    
    ! Integration parameters
    xend = 20.0_dp           ! Final time (enough to see chaotic behavior)
    
    ! Tolerance settings for chaotic system
    rtol = 1.0d-8           ! Relative tolerance (tight for chaos sensitivity)
    atol = 1.0d-10          ! Absolute tolerance (very tight for accurate trajectories)
    
    ! User parameters (not used in this problem)
    rpar = 0.0_dp
    ipar = 0
    
    ! Initialize work arrays
    do i = 1, min(10, lwork)
        work(i) = 0.0_dp
    end do
    
    do i = 1, min(10, liwork)
        iwork(i) = 0
    end do
    
    ! Set dense output for all components
    iwork(5) = n            ! Number of components for dense output
    iwork(4) = 100          ! Stiffness test frequency (higher for chaotic systems)
    
    write(6, '(A)') 'Problem setup completed.'
    write(6, '(A,3F10.6)') 'Initial conditions [x,y,z]: ', y(1), y(2), y(3)
    write(6, '(A,F8.3)') 'Integration from t = 0 to t = ', xend
    write(6, '(A,ES10.3,A,ES10.3)') 'Tolerances: RTOL = ', rtol, ', ATOL = ', atol
    write(6, '(A)') ''
    
end subroutine setup_lorenz_problem

! Main DOP853 integration routine
subroutine run_dop853_integration(f, solout, n, y, x, xend, rtol, atol, itol, &
                                 iout, work, lwork, iwork, liwork, rpar, ipar, idid)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    external f, solout
    integer, intent(in) :: n, itol, iout, lwork, liwork, ipar
    real(kind=dp), intent(inout) :: x, rpar
    real(kind=dp), intent(in) :: xend, rtol, atol
    real(kind=dp), dimension(n), intent(inout) :: y
    real(kind=dp), dimension(lwork), intent(inout) :: work
    integer, dimension(liwork), intent(inout) :: iwork
    integer, intent(out) :: idid
    
    write(6, '(A)') 'Starting DOP853 integration for Lorenz chaotic system...'
    write(6, '(A)') ' Time        X             Y             Z           NSTEP'
    write(6, '(A)') '----------------------------------------------------------------'
    
    ! Call DOP853 solver
    call dop853(n, f, x, y, xend, &
                rtol, atol, itol, &
                solout, iout, &
                work, lwork, iwork, liwork, rpar, ipar, idid)
    
    ! Check integration status
    select case (idid)
    case (1)
        write(6, '(/A)') 'Integration completed successfully.'
    case (2)
        write(6, '(/A)') 'Integration completed - interrupted by solout.'
    case (-1)
        write(6, '(/A)') 'Error: Input is not consistent.'
    case (-2)
        write(6, '(/A)') 'Error: Larger NMAX is needed.'
    case (-3)
        write(6, '(/A)') 'Error: Step size becomes too small.'
    case (-4)
        write(6, '(/A)') 'Error: Problem is probably stiff (use stiff solver).'
    case default
        write(6, '(/A,I0)') 'Integration finished with status: ', idid
    end select
    
end subroutine run_dop853_integration

! Print final statistics
subroutine print_dop853_statistics(n, x, y, rtol, liwork, iwork, idid)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: n, liwork
    real(kind=dp), intent(in) :: x, rtol
    real(kind=dp), dimension(n), intent(in) :: y
    integer, dimension(liwork), intent(in) :: iwork
    integer, intent(in) :: idid
    real(kind=dp) :: energy_like
    
    write(6, '(/A)') 'Final Integration Statistics:'
    write(6, '(A)') '============================'
    write(6, '(A,F10.6,A,3F12.6)') 'Final time: ', x, ', Final state [x,y,z]: ', y(1), y(2), y(3)
    
    ! Calculate a Lorenz "energy-like" quantity for monitoring
    energy_like = y(1)**2 + y(2)**2 + y(3)**2
    write(6, '(A,F12.6)') 'Final ||r||² (x² + y² + z²): ', energy_like
    
    write(6, '(A,ES10.3)') 'Relative tolerance used: ', rtol
    write(6, '(A,I0)') 'Integration status (IDID): ', idid
    write(6, '(/A)') 'Detailed statistics:'
    write(6, '(A,I6)') ' Function evaluations: ', iwork(17)
    write(6, '(A,I6)') ' Integration steps: ', iwork(18)
    write(6, '(A,I6)') ' Accepted steps: ', iwork(19)
    write(6, '(A,I6)') ' Rejected steps: ', iwork(20)
    
    if (iwork(18) > 0) then
        write(6, '(A,F8.4)') ' Average step size: ', x / real(iwork(18), dp)
        write(6, '(A,F6.2,A)') ' Acceptance ratio: ', &
                               100.0_dp * real(iwork(19), dp) / real(iwork(18), dp), '%'
    endif
    
    ! Analysis for chaotic behavior
    write(6, '(/A)') 'Chaotic System Analysis:'
    write(6, '(A)') '========================'
    if (abs(y(1)) > 20.0_dp .or. abs(y(2)) > 30.0_dp .or. abs(y(3)) > 50.0_dp) then
        write(6, '(A)') ' System has explored the chaotic attractor extensively.'
    else
        write(6, '(A)') ' System may still be in transient phase or bounded region.'
    endif
    
    write(6, '(A,F8.2)') ' Final distance from origin: ', sqrt(energy_like)
    
end subroutine print_dop853_statistics

! Output routine for continuous solution monitoring
subroutine solout_lorenz(nr, xold, x, y, n, con, icomp, nd, rpar, ipar, irtrn, xout)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: nr, n, nd, ipar
    integer, dimension(nd), intent(in) :: icomp
    real(kind=dp), intent(in) :: xold, rpar
    real(kind=dp), intent(in) :: x
    real(kind=dp), dimension(n), intent(in) :: y
    real(kind=dp), dimension(8*nd), intent(in) :: con
    integer, intent(inout) :: irtrn
    real(kind=dp), intent(inout) :: xout
    real(kind=dp), save :: dt = 0.5_dp  ! Output interval for oscillations
    
    ! External function for dense output
    real(kind=dp), external :: contd8
    
    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) xold, rpar, ipar, irtrn, icomp(1)
    endif
    
    if (nr == 1) then
        ! Print initial condition
        write(6, '(1X,F10.6,3F14.6,I8)') x, y(1), y(2), y(3), nr-1
        xout = dt  ! Output interval (smaller for chaotic systems)
    else
        ! Print at regular intervals
        do while (x >= xout)
            write(6, '(1X,F10.6,3F14.6,I8)') xout, &
                   contd8(1, xout, con, icomp, nd), &
                   contd8(2, xout, con, icomp, nd), &
                   contd8(3, xout, con, icomp, nd), nr-1
            xout = xout + dt
        end do
        
        ! Check if we need to print the final time point (only if not already printed)
        if (abs(x - (xout - dt)) < 1.0d-10 .and. x > xout - dt) then
            write(6, '(1X,F10.6,3F14.6,I8)') x, y(1), y(2), y(3), nr-1
        endif
    endif
    
end subroutine solout_lorenz

! Right-hand side: Lorenz system equations
subroutine florenz(n, x, y, f, rpar, ipar)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: n, ipar
    real(kind=dp), intent(in) :: x, rpar
    real(kind=dp), dimension(n), intent(in) :: y
    real(kind=dp), dimension(n), intent(out) :: f
    
    ! Lorenz parameters (classic chaotic regime)
    real(kind=dp), parameter :: sigma = 10.0_dp     ! Prandtl number
    real(kind=dp), parameter :: rho = 28.0_dp       ! Rayleigh number
    real(kind=dp), parameter :: beta = 8.0_dp/3.0_dp ! Geometric factor
    
    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) x, rpar, ipar
    endif
    
    ! Lorenz equations:
    ! dx/dt = σ(y - x)
    ! dy/dt = x(ρ - z) - y
    ! dz/dt = xy - βz
    
    f(1) = sigma * (y(2) - y(1))
    f(2) = y(1) * (rho - y(3)) - y(2)
    f(3) = y(1) * y(2) - beta * y(3)
    
end subroutine florenz