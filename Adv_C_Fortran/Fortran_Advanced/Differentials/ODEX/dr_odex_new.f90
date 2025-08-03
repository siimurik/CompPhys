! Compile and execute with:
! $ gfortran vanderpol_odex.f90 odex.f -o vanderpol_odex
! $ ./vanderpol_odex

! ODEX solver for stiff differential equations
! Solves the Van der Pol equation - a classic stiff oscillator:
! x'' + ε(x² - 1)x' + x = 0
! Rewritten as system:
! dx/dt = y
! dy/dt = (1 - x²)y - x) / ε
! For small ε, this becomes very stiff

program odex_vanderpol
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    external fvpol, solout_vanderpol, f
    
    ! Problem parameters
    integer, parameter :: neq = 2, km = 9, nrdens = 2
    real(kind=dp), dimension(neq) :: y
    real(kind=dp) :: x, xend, rtol, atol, rpar, h
    integer :: n, iout, itol, ipar, idid
    integer :: lwork, liwork
    real(kind=dp), allocatable, dimension(:) :: work
    integer, allocatable, dimension(:) :: iwork
    
    ! Calculate required array sizes and allocate
    call calculate_odex_arrays(neq, km, nrdens, lwork, liwork)
    allocate(work(lwork))
    allocate(iwork(liwork))
    
    ! Initialize problem parameters
    call setup_vanderpol_problem(neq, y, x, xend, rtol, atol, rpar, ipar, &
                                 n, iout, itol, h, work, iwork, lwork, liwork, &
                                 nrdens)
    
    ! Run the integration
    call run_odex_integration(fvpol, solout_vanderpol, n, y, x, xend, h, &
                             rtol, atol, itol, iout, work, lwork, &
                             iwork, liwork, rpar, ipar, idid)
    
    ! Print final statistics
    call print_odex_statistics(n, x, y, rtol, rpar, liwork, iwork, idid)
    
    ! Clean up
    deallocate(work)
    deallocate(iwork)
    
end program odex_vanderpol

! Subroutine to calculate ODEX array lengths
subroutine calculate_odex_arrays(neq, km, nrdens, lwork, liwork)
    implicit none
    integer, intent(in) :: neq, km, nrdens
    integer, intent(out) :: lwork, liwork
    
    ! ODEX workspace requirements
    lwork = neq * (km + 5) + 5 * km + 20 + (2 * km * (km + 2) + 5) * nrdens
    liwork = 2 * km + 21 + nrdens
    
    write(6, '(A,I0,A,I0)') 'Calculated ODEX array sizes: LWORK = ', lwork, ', LIWORK = ', liwork
    write(6, '(A,I0,A,I0)') 'Using KM = ', km, ', NRDENS = ', nrdens
    
end subroutine calculate_odex_arrays

! Setup the Van der Pol stiff oscillator problem
subroutine setup_vanderpol_problem(neq, y, x, xend, rtol, atol, rpar, ipar, &
                                   n, iout, itol, h, work, iwork, lwork, liwork, &
                                   nrdens)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: neq, lwork, liwork, nrdens
    real(kind=dp), dimension(neq), intent(out) :: y
    real(kind=dp), intent(out) :: x, xend, rtol, atol, rpar, h
    integer, intent(out) :: n, iout, itol, ipar
    real(kind=dp), dimension(lwork), intent(inout) :: work
    integer, dimension(liwork), intent(inout) :: iwork
    integer :: i
    
    write(6, '(A)') 'Setting up Van der Pol stiff oscillator...'
    write(6, '(A)') 'Equation: x'' + ε(x² - 1)x'' + x = 0'
    write(6, '(A)') 'System form:'
    write(6, '(A)') '  dx/dt = y'
    write(6, '(A)') '  dy/dt = ((1 - x²)y - x) / ε'
    write(6, '(A)') ''
    
    ! Problem dimension
    n = neq
    
    ! ODEX options
    iout = 2      ! Dense output is used during integration
    itol = 0      ! Scalar tolerances
    
    ! Initial conditions: limit cycle approach
    x = 0.0_dp               ! Initial time
    y(1) = 2.0_dp            ! x = 2.0 (displaced from equilibrium)
    y(2) = 0.0_dp            ! y = 0.0 (initially at rest)
    
    ! Integration parameters
    xend = 20.0_dp           ! Final time (multiple periods for limit cycle)
    
    ! Van der Pol parameter (stiffness parameter)
    rpar = 0.01_dp           ! ε = 0.01 (stiff case - smaller values more stiff)
    
    ! Tolerance settings for stiff system
    rtol = 1.0d-6           ! Relative tolerance
    atol = 1.0d-8           ! Absolute tolerance
    
    ! Initial step size
    h = 0.001_dp            ! Small initial step for stiff problem
    
    ! User parameters
    ipar = 0
    
    ! Initialize work arrays
    do i = 1, min(10, lwork)
        work(i) = 0.0_dp
    end do
    
    do i = 1, min(21, liwork)
        iwork(i) = 0
    end do
    
    ! Set dense output for both components
    iwork(8) = nrdens       ! Number of components for dense output
    iwork(21) = 1           ! First component for dense output
    iwork(22) = 2           ! Second component for dense output
    
    write(6, '(A)') 'Problem setup completed.'
    write(6, '(A,2F10.6)') 'Initial conditions [x,y]: ', y(1), y(2)
    write(6, '(A,F8.3)') 'Integration from t = 0 to t = ', xend
    write(6, '(A,ES10.3)') 'Stiffness parameter ε = ', rpar
    write(6, '(A,ES10.3,A,ES10.3)') 'Tolerances: RTOL = ', rtol, ', ATOL = ', atol
    write(6, '(A,ES10.3)') 'Initial step size: ', h
    write(6, '(A)') ''
    
    if (rpar < 0.1_dp) then
        write(6, '(A)') 'Note: Small ε indicates stiff problem - ODEX is well-suited for this.'
    else
        write(6, '(A)') 'Note: Moderate ε - problem may not be very stiff.'
    endif
    write(6, '(A)') ''
    
end subroutine setup_vanderpol_problem

! Main ODEX integration routine
subroutine run_odex_integration(f, solout, n, y, x, xend, h, rtol, atol, itol, &
                               iout, work, lwork, iwork, liwork, rpar, ipar, idid)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    external f, solout
    integer, intent(in) :: n, itol, iout, lwork, liwork, ipar
    real(kind=dp), intent(inout) :: x, h, rpar
    real(kind=dp), intent(in) :: xend, rtol, atol
    real(kind=dp), dimension(n), intent(inout) :: y
    real(kind=dp), dimension(lwork), intent(inout) :: work
    integer, dimension(liwork), intent(inout) :: iwork
    integer, intent(out) :: idid
    
    write(6, '(A)') 'Starting ODEX integration for Van der Pol equation...'
    write(6, '(A)') ' Time         X             Y           NSTEP'
    write(6, '(A)') '-------------------------------------------------------'
    
    ! Call ODEX solver
    call odex(n, f, x, y, xend, h, &
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
        write(6, '(/A)') 'Error: Matrix is repeatedly singular.'
    case default
        write(6, '(/A,I0)') 'Integration finished with status: ', idid
    end select
    
end subroutine run_odex_integration

! Print final statistics
subroutine print_odex_statistics(n, x, y, rtol, rpar, liwork, iwork, idid)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: n, liwork
    real(kind=dp), intent(in) :: x, rtol, rpar
    real(kind=dp), dimension(n), intent(in) :: y
    integer, dimension(liwork), intent(in) :: iwork
    integer, intent(in) :: idid
    real(kind=dp) :: energy, amplitude, period_est
    
    write(6, '(/A)') 'Final Integration Statistics:'
    write(6, '(A)') '============================'
    write(6, '(A,F10.6,A,2F12.6)') 'Final time: ', x, ', Final state [x,y]: ', y(1), y(2)
    
    ! Calculate energy-like quantity for Van der Pol oscillator
    energy = 0.5_dp * (y(1)**2 + y(2)**2)
    amplitude = sqrt(y(1)**2 + y(2)**2)
    write(6, '(A,F12.6)') 'Final energy-like quantity (½(x² + y²)): ', energy
    write(6, '(A,F12.6)') 'Final amplitude ||(x,y)||: ', amplitude
    
    write(6, '(A,ES10.3)') 'Relative tolerance used: ', rtol
    write(6, '(A,ES10.3)') 'Stiffness parameter ε: ', rpar
    write(6, '(A,I0)') 'Integration status (IDID): ', idid
    write(6, '(/A)') 'Detailed statistics:'
    write(6, '(A,I6)') ' Function evaluations: ', iwork(17)
    write(6, '(A,I6)') ' Integration steps: ', iwork(18)
    write(6, '(A,I6)') ' Accepted steps: ', iwork(19)
    write(6, '(A,I6)') ' Rejected steps: ', iwork(20)
    
    if (iwork(18) > 0) then
        write(6, '(A,F8.6)') ' Average step size: ', x / real(iwork(18), dp)
        write(6, '(A,F6.2,A)') ' Acceptance ratio: ', &
                               100.0_dp * real(iwork(19), dp) / real(iwork(18), dp), '%'
    endif
    
    ! Analysis for Van der Pol behavior
    write(6, '(/A)') 'Van der Pol Analysis:'
    write(6, '(A)') '===================='
    
    ! Estimate period for limit cycle (approximate for small ε)
    if (rpar > 0.0_dp .and. rpar < 1.0_dp) then
        period_est = 2.0_dp * (3.0_dp - 2.0_dp * log(2.0_dp)) / rpar**0.5_dp
        write(6, '(A,F8.3)') ' Estimated limit cycle period: ', period_est
        write(6, '(A,F6.2)') ' Number of periods integrated: ', x / period_est
    endif
    
    if (amplitude > 1.5_dp .and. amplitude < 2.5_dp) then
        write(6, '(A)') ' System appears to be on or near the limit cycle.'
    elseif (amplitude < 0.5_dp) then
        write(6, '(A)') ' System is near the unstable equilibrium at origin.'
    else
        write(6, '(A)') ' System may be in transient phase toward limit cycle.'
    endif
    
    write(6, '(A,F8.4)') ' Final phase space radius: ', amplitude
    
    ! Stiffness analysis
    if (rpar < 0.01_dp) then
        write(6, '(A)') ' Very stiff problem - ODEX should handle efficiently.'
    elseif (rpar < 0.1_dp) then
        write(6, '(A)') ' Moderately stiff problem.'
    else
        write(6, '(A)') ' Mildly stiff or non-stiff problem.'
    endif
    
end subroutine print_odex_statistics

! Output routine for continuous solution monitoring
subroutine solout_vanderpol(nr, xold, x, y, n, con, ncon, icomp, nd, rpar, ipar, irtrn)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: nr, n, ncon, nd, ipar
    integer, dimension(nd), intent(in) :: icomp
    real(kind=dp), intent(in) :: xold, rpar
    real(kind=dp), intent(in) :: x
    real(kind=dp), dimension(n), intent(in) :: y
    real(kind=dp), dimension(ncon), intent(in) :: con
    integer, intent(inout) :: irtrn
    real(kind=dp), save :: xout = 0.0_dp
    real(kind=dp), parameter :: dt = 0.5_dp  ! Output interval
    
    ! External function for dense output
    real(kind=dp), external :: contex
    
    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) xold, ipar, irtrn, rpar
    endif
    
    if (nr == 1) then
        ! Print initial condition
        write(6, '(1X,F12.6,2F14.6,I8)') x, y(1), y(2), nr-1
        xout = dt  ! Set first output time
    else
        ! Print at regular intervals using dense output
        do while (x >= xout)
            write(6, '(1X,F12.6,2F14.6,I8)') xout, &
                   contex(1, xout, con, ncon, icomp, nd), &
                   contex(2, xout, con, ncon, icomp, nd), nr-1
            xout = xout + dt
        end do
        
        ! Print final point if we haven't already
        if (abs(x - (xout - dt)) > 1.0d-10) then
            write(6, '(1X,F12.6,2F14.6,I8)') x, y(1), y(2), nr-1
        endif
    endif
    
end subroutine solout_vanderpol

! Right-hand side: Van der Pol equation
subroutine fvpol(n, x, y, f, rpar, ipar)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: n, ipar
    real(kind=dp), intent(in) :: x, rpar
    real(kind=dp), dimension(n), intent(in) :: y
    real(kind=dp), dimension(n), intent(out) :: f
    
    real(kind=dp) :: eps
    
    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) x, ipar
    endif
    
    ! Van der Pol parameter (stiffness parameter)
    eps = rpar
    
    ! Van der Pol equations:
    ! dx/dt = y
    ! dy/dt = ((1 - x²)y - x) / ε
    
    f(1) = y(2)
    f(2) = ((1.0_dp - y(1)**2) * y(2) - y(1)) / eps
    
end subroutine fvpol