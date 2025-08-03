! Compile and execute with:
! $ gfortran brusselator_dopri5.f90 dopri5.f -o brusselator_dopri5
! $ ./brusselator_dopri5

! DOPRI5 solver for oscillatory chemical reaction systems
! Solves the Brusselator system - a classic autocatalytic oscillator:
! A -> X (rate k1*A)
! 2X + Y -> 3X (rate k2*X²*Y)  
! B + X -> Y + D (rate k3*B*X)
! X -> E (rate k4*X)
! Leading to: dx/dt = A - (B+1)x + x²y, dy/dt = Bx - x²y

program dopri5_brusselator
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    external fbrussel, solout_brussels, f
    
    ! Problem parameters
    integer, parameter :: neq = 2, nrdens = 2
    real(kind=dp), dimension(neq) :: y
    real(kind=dp), dimension(2) :: rpar
    real(kind=dp) :: x, xend, rtol, atol
    integer :: n, iout, itol, ipar, idid
    integer :: lwork, liwork
    real(kind=dp), allocatable, dimension(:) :: work
    integer, allocatable, dimension(:) :: iwork
    
    ! Calculate required array sizes and allocate
    call calculate_dopri5_arrays(neq, nrdens, lwork, liwork)
    allocate(work(lwork))
    allocate(iwork(liwork))
    
    ! Initialize problem parameters
    call setup_brusselator_problem(neq, nrdens, y, x, xend, rtol, atol, rpar, ipar, &
                                  n, iout, itol, work, iwork, lwork, liwork)
    
    ! Run the integration
    call run_dopri5_integration(fbrussel, solout_brussels, n, y, x, xend, &
                               rtol, atol, itol, iout, work, lwork, &
                               iwork, liwork, rpar, ipar, idid)
    
    ! Print final statistics
    call print_dopri5_statistics(n, x, y, rtol, liwork, iwork, idid)
    
    ! Clean up
    deallocate(work)
    deallocate(iwork)
    
end program dopri5_brusselator

! Subroutine to calculate DOPRI5 array lengths
subroutine calculate_dopri5_arrays(neq, nrdens, lwork, liwork)
    implicit none
    integer, intent(in) :: neq, nrdens
    integer, intent(out) :: lwork, liwork
    
    ! For dense output
    lwork = 8 * neq + 5 * nrdens + 21
    liwork = nrdens + 21
    
    write(6, '(A,I0,A,I0)') 'Calculated DOPRI5 array sizes: LWORK = ', lwork, ', LIWORK = ', liwork
    
end subroutine calculate_dopri5_arrays

! Setup the Brusselator oscillatory system problem
subroutine setup_brusselator_problem(neq, nrdens, y, x, xend, rtol, atol, rpar, ipar, &
                                    n, iout, itol, work, iwork, lwork, liwork)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: neq, nrdens, lwork, liwork
    real(kind=dp), dimension(neq), intent(out) :: y
    real(kind=dp), dimension(2), intent(out) :: rpar
    real(kind=dp), intent(out) :: x, xend, rtol, atol
    integer, intent(out) :: n, iout, itol, ipar
    real(kind=dp), dimension(lwork), intent(inout) :: work
    integer, dimension(liwork), intent(inout) :: iwork
    integer :: i
    
    write(6, '(A)') 'Setting up Brusselator oscillatory chemical system...'
    write(6, '(A)') 'Chemical reactions:'
    write(6, '(A)') '  A -> X         (rate k₁A)'
    write(6, '(A)') '  2X + Y -> 3X   (rate k₂X²Y)'
    write(6, '(A)') '  B + X -> Y + D (rate k₃BX)'
    write(6, '(A)') '  X -> E         (rate k₄X)'
    write(6, '(A)') ''
    write(6, '(A)') 'Dimensionless equations:'
    write(6, '(A)') '  dx/dt = A - (B+1)x + x²y'
    write(6, '(A)') '  dy/dt = Bx - x²y'
    write(6, '(A)') ''
    
    ! Problem dimension
    n = neq
    
    ! DOPRI5 options
    iout = 2      ! Dense output used during integration
    itol = 0      ! Scalar tolerances
    
    ! Brusselator parameters (oscillatory regime)
    rpar(1) = 1.0_dp     ! A parameter
    rpar(2) = 3.0_dp     ! B parameter (B > 1 + A² for oscillations)
    
    write(6, '(A,F6.3,A,F6.3)') 'Parameters: A = ', rpar(1), ', B = ', rpar(2)
    write(6, '(A)') 'This gives oscillatory behavior (B > 1 + A² = 2.0)'
    write(6, '(A)') ''
    
    ! Initial conditions: near steady state with small perturbation
    x = 0.0_dp                ! Initial time
    y(1) = 1.5_dp            ! x concentration (slightly perturbed from steady state)
    y(2) = 3.0_dp            ! y concentration (slightly perturbed from steady state)
    
    ! Integration parameters (multiple oscillation periods)
    xend = 20.0_dp           ! Final time (several oscillation periods)
    
    ! Tolerance settings for oscillatory system
    rtol = 1.0d-7           ! Relative tolerance
    atol = 1.0d-9           ! Absolute tolerance
    
    ! User integer parameter (not used)
    ipar = 0
    
    ! Initialize work arrays
    do i = 1, min(20, lwork)
        work(i) = 0.0_dp
    end do
    
    do i = 1, min(20, liwork)
        iwork(i) = 0
    end do
    
    ! Set dense output for both components
    iwork(5) = nrdens       ! Number of components for dense output
    iwork(21) = 1          ! Component 1 (x) for dense output
    iwork(22) = 2          ! Component 2 (y) for dense output
    
    write(6, '(A)') 'Problem setup completed.'
    write(6, '(A,2F10.6)') 'Initial concentrations [x,y]: ', y(1), y(2)
    write(6, '(A,F8.3)') 'Integration from t = 0 to t = ', xend
    write(6, '(A,ES10.3,A,ES10.3)') 'Tolerances: RTOL = ', rtol, ', ATOL = ', atol
    write(6, '(A)') ''
    
end subroutine setup_brusselator_problem

! Main DOPRI5 integration routine
subroutine run_dopri5_integration(f, solout, n, y, x, xend, rtol, atol, itol, &
                                 iout, work, lwork, iwork, liwork, rpar, ipar, idid)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    external f, solout
    integer, intent(in) :: n, itol, iout, lwork, liwork, ipar
    real(kind=dp), intent(inout) :: x
    real(kind=dp), intent(in) :: xend, rtol, atol
    real(kind=dp), dimension(n), intent(inout) :: y
    real(kind=dp), dimension(lwork), intent(inout) :: work
    integer, dimension(liwork), intent(inout) :: iwork
    real(kind=dp), dimension(2), intent(inout) :: rpar
    integer, intent(out) :: idid
    
    write(6, '(A)') 'Starting DOPRI5 integration for Brusselator oscillatory system...'
    write(6, '(A)') '     Time           X             Y         NSTEP'
    write(6, '(A)') '-------------------------------------------------'
    
    ! Call DOPRI5 solver
    call dopri5(n, f, x, y, xend, &
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
    
end subroutine run_dopri5_integration

! Print final statistics
subroutine print_dopri5_statistics(n, x, y, rtol, liwork, iwork, idid)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: n, liwork
    real(kind=dp), intent(in) :: x, rtol
    real(kind=dp), dimension(n), intent(in) :: y
    integer, dimension(liwork), intent(in) :: iwork
    integer, intent(in) :: idid
    real(kind=dp) :: period_estimate
    
    write(6, '(/A)') 'Final Integration Statistics:'
    write(6, '(A)') '============================'
    write(6, '(A,F10.6,A,2F10.6)') 'Final time: ', x, ', Final concentrations [x,y]:', y(1), y(2)
    
    write(6, '(A,ES10.3)') 'Relative tolerance used: ', rtol
    write(6, '(A,I0)') 'Integration status (IDID): ', idid
    write(6, '(/A)') 'Detailed statistics:'
    write(6, '(A,I5)') ' Function evaluations: ', iwork(17)
    write(6, '(A,I5)') ' Integration steps:    ', iwork(18)
    write(6, '(A,I5)') ' Accepted steps:       ', iwork(19)
    write(6, '(A,I5)') ' Rejected steps:       ', iwork(20)
    
    if (iwork(18) > 0) then
        write(6, '(A,F8.4)') ' Average step size: ', x / real(iwork(18), dp)
        write(6, '(A,F8.2,A)') ' Acceptance ratio: ', &
                               100.0_dp * real(iwork(19), dp) / real(iwork(18), dp), '%'
    endif
    
    ! Analysis for oscillatory behavior
    write(6, '(/A)') 'Oscillatory System Analysis:'
    write(6, '(A)') '============================'
    
    ! Rough period estimate for Brusselator (theoretical ≈ 2π/√(A) for A=1, B=3)
    period_estimate = 2.0_dp * 4.0_dp * atan(1.0_dp)  ! 2π approximation
    write(6, '(A,F6.3)') ' Theoretical period estimate: ', period_estimate
    write(6, '(A,I0)') ' Number of periods integrated: ', int(x / period_estimate)
    
    ! Check if system shows expected oscillatory range
    if (y(1) > 0.5_dp .and. y(1) < 4.0_dp .and. y(2) > 1.0_dp .and. y(2) < 6.0_dp) then
        write(6, '(A)') ' Final state within expected oscillatory range.'
    else
        write(6, '(A)') ' Final state may be outside typical oscillatory range.'
    endif
    
end subroutine print_dopri5_statistics

! Output routine for continuous solution monitoring
subroutine solout_brussels(nr, xold, x, y, n, con, icomp, nd, rpar, ipar, irtrn)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: nr, n, nd, ipar
    integer, dimension(nd), intent(in) :: icomp
    real(kind=dp), intent(in) :: xold
    real(kind=dp), intent(in) :: x
    real(kind=dp), dimension(n), intent(in) :: y
    real(kind=dp), dimension(5*nd), intent(in) :: con
    real(kind=dp), dimension(2), intent(in) :: rpar
    integer, intent(inout) :: irtrn
    
    ! External function for dense output
    real(kind=dp), external :: contd5
    
    ! Storage for output control
    real(kind=dp), save :: xout = 0.0_dp
    logical, save :: first_call = .true.
    
    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) xold, rpar(1), ipar, irtrn
    endif
    
    if (first_call) then
        ! Print initial condition
        write(6, '(1X,F10.6,2F14.6,I8)') x, y(1), y(2), nr-1
        xout = 1.0_dp  ! Output interval (good for oscillations)
        first_call = .false.
    else
        ! Print at regular intervals
        do while (x >= xout)
            write(6, '(1X,F10.6,2F14.6,I8)') xout, &
                   contd5(1, xout, con, icomp, nd), &
                   contd5(2, xout, con, icomp, nd), nr-1
            xout = xout + 1.0_dp
        end do
        
        ! Check if we need to print the final time point (only if not already printed)
        if (abs(x - (xout - 1.0_dp)) < 1.0d-10 .and. x > xout - 1.0_dp) then
            write(6, '(1X,F10.6,2F14.6,I8)') x, y(1), y(2), nr-1
        endif
    endif
    
end subroutine solout_brussels

! Right-hand side: Brusselator system equations
subroutine fbrussel(n, x, y, f, rpar, ipar)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: n, ipar
    real(kind=dp), intent(in) :: x
    real(kind=dp), dimension(n), intent(in) :: y
    real(kind=dp), dimension(n), intent(out) :: f
    real(kind=dp), dimension(2), intent(in) :: rpar
    
    ! Brusselator parameters
    real(kind=dp) :: A, B
    
    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) x, ipar
    endif
    
    A = rpar(1)  ! A parameter
    B = rpar(2)  ! B parameter
    
    ! Brusselator equations (dimensionless form):
    ! dx/dt = A - (B+1)x + x²y
    ! dy/dt = Bx - x²y
    
    f(1) = A - (B + 1.0_dp) * y(1) + y(1) * y(1) * y(2)
    f(2) = B * y(1) - y(1) * y(1) * y(2)
    
end subroutine fbrussel