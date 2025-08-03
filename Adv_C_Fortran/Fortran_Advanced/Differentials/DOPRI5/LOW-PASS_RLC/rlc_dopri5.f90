! Compile and execute with:
! $ gfortran rlc_dopri5.f90 dopri5.f90 -o rlc_dopri5
! $ ./rlc_dopri5

! DOPRI5 solver for RLC low-pass filter circuit
! Circuit equations:
! V = i_R * R
! i_C = C * dV/dt  
! L * di_L/dt = V_b - V
! i_L = i_R + i_C
!
! System of ODEs:
! dV/dt = (1/C) * (i_L - V/R)
! di_L/dt = (1/L) * (V_b - V)

program dopri5_rlc_filter
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    external frlc, solout_rlc, f
    
    ! Problem parameters
    integer, parameter :: neq = 2, nrdens = 2
    real(kind=dp), dimension(neq) :: y
    real(kind=dp), dimension(4) :: rpar  ! [Vb, L, R, C]
    real(kind=dp) :: t, tend, rtol, atol
    integer :: n, iout, itol, ipar, idid
    integer :: lwork, liwork
    real(kind=dp), allocatable, dimension(:) :: work
    integer, allocatable, dimension(:) :: iwork
    
    ! Calculate required array sizes and allocate
    call calculate_dopri5_arrays(neq, nrdens, lwork, liwork)
    allocate(work(lwork))
    allocate(iwork(liwork))
    
    ! Initialize problem parameters
    call setup_rlc_problem(neq, nrdens, y, t, tend, rtol, atol, rpar, ipar, &
                          n, iout, itol, work, iwork, lwork, liwork)
    
    ! Run the integration
    call run_dopri5_integration(frlc, solout_rlc, n, y, t, tend, &
                               rtol, atol, itol, iout, work, lwork, &
                               iwork, liwork, rpar, ipar, idid)
    
    ! Print final statistics
    call print_dopri5_statistics(n, t, y, rtol, liwork, iwork, idid, rpar)
    
    ! Clean up
    deallocate(work)
    deallocate(iwork)
    
    ! Close output file
    close(10)
    
    write(6, '(/A)') 'Data written to "rlc_data.csv" for plotting.'
    write(6, '(A)') 'Run "python3 plot_rlc.py" to generate plots.'
    
end program dopri5_rlc_filter

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

! Setup the RLC low-pass filter problem
subroutine setup_rlc_problem(neq, nrdens, y, t, tend, rtol, atol, rpar, ipar, &
                             n, iout, itol, work, iwork, lwork, liwork)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: neq, nrdens, lwork, liwork
    real(kind=dp), dimension(neq), intent(out) :: y
    real(kind=dp), dimension(4), intent(out) :: rpar
    real(kind=dp), intent(out) :: t, tend, rtol, atol
    integer, intent(out) :: n, iout, itol, ipar
    real(kind=dp), dimension(lwork), intent(inout) :: work
    integer, dimension(liwork), intent(inout) :: iwork
    integer :: i
    real(kind=dp) :: omega_n, zeta, tau
    
    write(6, '(A)') 'Setting up RLC Low-Pass Filter Circuit...'
    write(6, '(A)') 'Circuit equations:'
    write(6, '(A)') '  V = i_R * R'
    write(6, '(A)') '  i_C = C * dV/dt'
    write(6, '(A)') '  L * di_L/dt = V_b - V'
    write(6, '(A)') '  i_L = i_R + i_C'
    write(6, '(A)') ''
    write(6, '(A)') 'System of ODEs:'
    write(6, '(A)') '  dV/dt = (1/C) * (i_L - V/R)'
    write(6, '(A)') '  di_L/dt = (1/L) * (V_b - V)'
    write(6, '(A)') ''
    
    ! Problem dimension
    n = neq
    
    ! DOPRI5 options
    iout = 2      ! Dense output used during integration
    itol = 0      ! Scalar tolerances
    
    ! RLC circuit parameters
    rpar(1) = 24.0_dp      ! V_b (input voltage)
    rpar(2) = 1.0_dp       ! L (inductance in H)
    rpar(3) = 100.0_dp     ! R (resistance in Ohms)
    rpar(4) = 1.0e-3_dp    ! C (capacitance in F)
    
    write(6, '(A,F6.1,A)') 'Parameters: V_b = ', rpar(1), ' V'
    write(6, '(A,F6.1,A)') '           L = ', rpar(2), ' H'
    write(6, '(A,F6.1,A)') '           R = ', rpar(3), ' Ω'
    write(6, '(A,ES8.1,A)') '           C = ', rpar(4), ' F'
    
    ! Calculate system characteristics
    omega_n = 1.0_dp / sqrt(rpar(2) * rpar(4))  ! Natural frequency
    zeta = 0.5_dp * rpar(3) * sqrt(rpar(4) / rpar(2))  ! Damping ratio
    tau = rpar(2) / rpar(3)  ! L/R time constant
    
    write(6, '(A,F8.2,A)') '           ω_n = ', omega_n, ' rad/s (natural frequency)'
    write(6, '(A,F6.3)') '           ζ = ', zeta, ' (damping ratio)'
    write(6, '(A,ES8.1,A)') '           τ = L/R = ', tau, ' s (time constant)'
    
    if (zeta < 1.0_dp) then
        write(6, '(A)') '           System is UNDERDAMPED (oscillatory)'
    elseif (abs(zeta - 1.0_dp) < 1.0e-10_dp) then
        write(6, '(A)') '           System is CRITICALLY DAMPED'
    else
        write(6, '(A)') '           System is OVERDAMPED'
    endif
    write(6, '(A)') ''
    
    ! Initial conditions (circuit starts at rest)
    t = 0.0_dp              ! Initial time
    y(1) = 0.0_dp           ! V (initial voltage across capacitor)
    y(2) = 0.0_dp           ! i_L (initial current through inductor)
    
    ! Integration time (several time constants to see settling)
    !tend = 5.0_dp * tau     ! Integrate for 5 time constants
    tend = 2.0_dp
    if (tend < 0.01_dp) tend = 0.01_dp  ! Minimum integration time
    
    ! Tolerance settings
    rtol = 1.0d-8           ! Relative tolerance
    atol = 1.0d-10          ! Absolute tolerance
    
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
    iwork(21) = 1          ! Component 1 (V) for dense output
    iwork(22) = 2          ! Component 2 (i_L) for dense output
    
    write(6, '(A)') 'Problem setup completed.'
    write(6, '(A,2F10.6)') 'Initial conditions [V, i_L]: ', y(1), y(2)
    write(6, '(A,F8.5)') 'Integration from t = 0 to t = ', tend
    write(6, '(A,ES10.3,A,ES10.3)') 'Tolerances: RTOL = ', rtol, ', ATOL = ', atol
    write(6, '(A)') ''
    
end subroutine setup_rlc_problem

! Main DOPRI5 integration routine
subroutine run_dopri5_integration(f, solout, n, y, t, tend, rtol, atol, itol, &
                                 iout, work, lwork, iwork, liwork, rpar, ipar, idid)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    external f, solout
    integer, intent(in) :: n, itol, iout, lwork, liwork, ipar
    real(kind=dp), intent(inout) :: t
    real(kind=dp), intent(in) :: tend, rtol, atol
    real(kind=dp), dimension(n), intent(inout) :: y
    real(kind=dp), dimension(lwork), intent(inout) :: work
    integer, dimension(liwork), intent(inout) :: iwork
    real(kind=dp), dimension(4), intent(inout) :: rpar
    integer, intent(out) :: idid
    
    write(6, '(A)') 'Starting DOPRI5 integration for RLC Low-Pass Filter...'
    write(6, '(A)') '     Time (s)      Voltage (V)    Current (A)     NSTEP'
    write(6, '(A)') '--------------------------------------------------------'
    
    ! Call DOPRI5 solver
    call dopri5(n, f, t, y, tend, &
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
subroutine print_dopri5_statistics(n, t, y, rtol, liwork, iwork, idid, rpar)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: n, liwork
    real(kind=dp), intent(in) :: t, rtol
    real(kind=dp), dimension(n), intent(in) :: y
    integer, dimension(liwork), intent(in) :: iwork
    integer, intent(in) :: idid
    real(kind=dp), dimension(4), intent(in) :: rpar
    real(kind=dp) :: settling_error, final_voltage, expected_voltage
    real(kind=dp) :: power_dissipated, energy_stored_L, energy_stored_C
    
    write(6, '(/A)') 'Final Integration Statistics:'
    write(6, '(A)') '============================'
    write(6, '(A,F10.6,A)') 'Final time: ', t, ' s'
    write(6, '(A,F10.6,A)') 'Final voltage: ', y(1), ' V'
    write(6, '(A,F10.6,A)') 'Final current: ', y(2), ' A'
    
    write(6, '(A,ES10.3)') 'Relative tolerance used: ', rtol
    write(6, '(A,I0)') 'Integration status (IDID): ', idid
    write(6, '(/A)') 'Detailed statistics:'
    write(6, '(A,I5)') ' Function evaluations: ', iwork(17)
    write(6, '(A,I5)') ' Integration steps:    ', iwork(18)
    write(6, '(A,I5)') ' Accepted steps:       ', iwork(19)
    write(6, '(A,I5)') ' Rejected steps:       ', iwork(20)
    
    if (iwork(18) > 0) then
        write(6, '(A,F10.6)') ' Average step size: ', t / real(iwork(18), dp)
        write(6, '(A,F8.2,A)') ' Acceptance ratio: ', &
                               100.0_dp * real(iwork(19), dp) / real(iwork(18), dp), '%'
    endif
    
    ! Circuit analysis
    write(6, '(/A)') 'Circuit Analysis:'
    write(6, '(A)') '=================='
    
    ! Expected steady-state voltage (should equal V_b for low-pass at DC)
    expected_voltage = rpar(1)  ! V_b
    final_voltage = y(1)
    settling_error = abs(final_voltage - expected_voltage) / expected_voltage * 100.0_dp
    
    write(6, '(A,F8.3,A)') ' Expected steady-state voltage: ', expected_voltage, ' V'
    write(6, '(A,F8.3,A)') ' Actual final voltage: ', final_voltage, ' V'
    write(6, '(A,F6.2,A)') ' Settling error: ', settling_error, '%'
    
    ! Power and energy analysis
    power_dissipated = y(1)**2 / rpar(3)  ! P = V²/R
    energy_stored_L = 0.5_dp * rpar(2) * y(2)**2  ! E_L = ½Li²
    energy_stored_C = 0.5_dp * rpar(4) * y(1)**2  ! E_C = ½CV²
    
    write(6, '(A,F8.3,A)') ' Power dissipated in R: ', power_dissipated, ' W'
    write(6, '(A,F8.6,A)') ' Energy stored in L: ', energy_stored_L, ' J'
    write(6, '(A,F8.6,A)') ' Energy stored in C: ', energy_stored_C, ' J'
    
    ! System behavior assessment
    if (settling_error < 5.0_dp) then
        write(6, '(A)') ' System has reached steady state (< 5% error).'
    else
        write(6, '(A)') ' System may need more time to settle.'
    endif
    
end subroutine print_dopri5_statistics

! Output routine for continuous solution monitoring
subroutine solout_rlc(nr, xold, t, y, n, con, icomp, nd, rpar, ipar, irtrn)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: nr, n, nd, ipar
    integer, dimension(nd), intent(in) :: icomp
    real(kind=dp), intent(in) :: xold
    real(kind=dp), intent(in) :: t
    real(kind=dp), dimension(n), intent(in) :: y
    real(kind=dp), dimension(5*nd), intent(in) :: con
    real(kind=dp), dimension(4), intent(in) :: rpar
    integer, intent(inout) :: irtrn
    
    ! External function for dense output
    real(kind=dp), external :: contd5
    
    ! Storage for output control
    real(kind=dp), save :: tout = 0.0_dp
    real(kind=dp), save :: dt_out = 0.0_dp
    logical, save :: first_call = .true.
    integer, parameter :: file_unit = 10
    real(kind=dp) :: V_current, iL_current, iR_current, iC_current, power_R, power_total
    
    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) xold, rpar(1), ipar, irtrn
    endif
    
    if (first_call) then
        ! Open output file for plotting data
        open(unit=file_unit, file='rlc_data.csv', status='replace')
        write(file_unit, '(A)') 'Time,Voltage,Current_L,Current_R,Current_C,Power_R,Power_Total'
        
        ! Print initial condition to screen
        write(6, '(1X,F12.6,2F14.6,I8)') t, y(1), y(2), nr-1
        
        ! Calculate derived quantities at initial time
        V_current = y(1)
        iL_current = y(2)
        iR_current = V_current / rpar(3)
        iC_current = iL_current - iR_current
        power_R = V_current * iR_current
        power_total = rpar(1) * iL_current
        
        ! Write initial condition to file
        write(file_unit, '(7(F16.8,","))') t, V_current, iL_current, iR_current, iC_current, power_R, power_total
        
        ! Set output interval based on time constant
        dt_out = rpar(2) / rpar(3) / 50.0_dp  ! τ/50 for high resolution plotting
        if (dt_out < 1.0e-6_dp) dt_out = 1.0e-6_dp  ! Minimum output interval
        
        tout = dt_out
        first_call = .false.
    else
        ! Print at regular intervals
        do while (t >= tout)
            V_current = contd5(1, tout, con, icomp, nd)
            iL_current = contd5(2, tout, con, icomp, nd)
            iR_current = V_current / rpar(3)
            iC_current = iL_current - iR_current
            power_R = V_current * iR_current
            power_total = rpar(1) * iL_current
            
            write(6, '(1X,F12.6,2F14.6,I8)') tout, V_current, iL_current, nr-1
            write(file_unit, '(7(F16.8,","))') tout, V_current, iL_current, iR_current, iC_current, power_R, power_total
            tout = tout + dt_out
        end do
    endif
    
end subroutine solout_rlc

! Right-hand side: RLC low-pass filter equations
subroutine frlc(n, t, y, f, rpar, ipar)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: n, ipar
    real(kind=dp), intent(in) :: t
    real(kind=dp), dimension(n), intent(in) :: y
    real(kind=dp), dimension(n), intent(out) :: f
    real(kind=dp), dimension(4), intent(in) :: rpar
    
    ! Circuit parameters
    real(kind=dp) :: Vb, L, R, C
    real(kind=dp) :: V, iL
    
    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) t, ipar
    endif
    
    Vb = rpar(1)  ! Input voltage
    L = rpar(2)   ! Inductance
    R = rpar(3)   ! Resistance
    C = rpar(4)   ! Capacitance
    
    V = y(1)      ! Voltage across capacitor
    iL = y(2)     ! Current through inductor
    
    ! RLC low-pass filter equations:
    ! dV/dt = (1/C) * (i_L - V/R)
    ! di_L/dt = (1/L) * (V_b - V)
    
    f(1) = (1.0_dp / C) * (iL - V / R)
    f(2) = (1.0_dp / L) * (Vb - V)
    
end subroutine frlc