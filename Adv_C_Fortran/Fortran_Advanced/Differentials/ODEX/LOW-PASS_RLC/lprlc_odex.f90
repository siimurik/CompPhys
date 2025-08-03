! Compile and execute with:
! $ gfortran lprlc_odex.f90 odex.f90 -o lprlc_odex
! $ ./lprlc_odex

! ODEX solver for RLC low-pass filter circuit
! Circuit equations:
! V = i_R * R
! i_C = C * dV/dt  
! L * di_L/dt = V_b - V
! i_L = i_R + i_C
!
! System of ODEs:
! dV/dt = (1/C) * (i_L - V/R)
! di_L/dt = (1/L) * (V_b - V)

program odex_rlc_filter
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    external frlc, solout_rlc, f
    
    ! Problem parameters
    integer, parameter :: neq = 2, km = 9, nrdens = 2
    real(kind=dp), dimension(neq) :: y
    real(kind=dp), dimension(4) :: rpar  ! [Vb, L, R, C]
    real(kind=dp) :: t, tend, rtol, atol, h
    integer :: n, iout, itol, ipar, idid
    integer :: lwork, liwork
    real(kind=dp), allocatable, dimension(:) :: work
    integer, allocatable, dimension(:) :: iwork
    
    ! Calculate required array sizes and allocate
    call calculate_odex_arrays(neq, km, nrdens, lwork, liwork)
    allocate(work(lwork))
    allocate(iwork(liwork))
    
    ! Initialize problem parameters
    call setup_rlc_problem(neq, nrdens, y, t, tend, rtol, atol, h, rpar, ipar, &
                           n, iout, itol, work, iwork, lwork, liwork)
    
    ! Run the integration
    call run_odex_integration(frlc, solout_rlc, n, y, t, tend, h, &
                             rtol, atol, itol, iout, work, lwork, &
                             iwork, liwork, rpar, ipar, idid)
    
    ! Print final statistics
    call print_odex_statistics(n, t, y, rtol, liwork, iwork, idid, rpar)
    
    ! Clean up
    deallocate(work)
    deallocate(iwork)
    
    ! Close output file
    close(10)
    
    write(6, '(/A)') 'Data written to "rlc_odex_data.csv" for plotting.'
    write(6, '(A)') 'Run "python3 plot_Vt.py" to generate plots.'
    
end program odex_rlc_filter

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

! Setup the RLC low-pass filter problem
subroutine setup_rlc_problem(neq, nrdens, y, t, tend, rtol, atol, h, rpar, ipar, &
                             n, iout, itol, work, iwork, lwork, liwork)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: neq, nrdens, lwork, liwork
    real(kind=dp), dimension(neq), intent(out) :: y
    real(kind=dp), dimension(4), intent(out) :: rpar
    real(kind=dp), intent(out) :: t, tend, rtol, atol, h
    integer, intent(out) :: n, iout, itol, ipar
    real(kind=dp), dimension(lwork), intent(inout) :: work
    integer, dimension(liwork), intent(inout) :: iwork
    integer :: i
    real(kind=dp) :: omega_n, zeta, tau, Q_factor, stiffness_ratio
    
    write(6, '(A)') 'Setting up RLC Low-Pass Filter Circuit with ODEX...'
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
    
    ! ODEX options
    iout = 2      ! Dense output used during integration
    itol = 0      ! Scalar tolerances
    
    ! RLC circuit parameters - choose values that might be stiff
    rpar(1) = 24.0_dp       ! V_b (input voltage)
    rpar(2) = 1.0_dp    ! L (inductance in H) - smaller L for higher frequencies
    rpar(3) = 100.0_dp       ! R (resistance in Ohms)
    rpar(4) = 1.0e-3_dp   ! C (capacitance in F) - smaller C for higher frequencies
    
    write(6, '(A,F6.1,A)') 'Parameters: V_b = ', rpar(1), ' V'
    write(6, '(A,ES8.1,A)') '           L = ', rpar(2), ' H'
    write(6, '(A,F6.1,A)')  '           R = ', rpar(3), ' Ω'
    write(6, '(A,ES8.1,A)') '           C = ', rpar(4), ' F'
    
    ! Calculate system characteristics
    omega_n = 1.0_dp / sqrt(rpar(2) * rpar(4))  ! Natural frequency
    zeta = 0.5_dp * rpar(3) * sqrt(rpar(4) / rpar(2))  ! Damping ratio
    tau = rpar(2) / rpar(3)  ! L/R time constant
    Q_factor = omega_n * rpar(2) / rpar(3)  ! Quality factor
    
    write(6, '(A,F10.1,A)') '         ω_n = ', omega_n, ' rad/s (natural frequency)'
    write(6, '(A,F8.1,A)')  '         f_n = ', omega_n/(2.0_dp*3.14159_dp), ' Hz'
    write(6, '(A,F8.3)')    '           ζ = ', zeta, ' (damping ratio)'
    write(6, '(A,F8.3)')    '           Q = ', Q_factor, ' (quality factor)'
    write(6, '(A,ES8.1,A)') '     τ = L/R = ', tau, ' s (time constant)'
    
    ! Assess stiffness
    stiffness_ratio = omega_n * tau
    write(6, '(A,F8.2)') '           Stiffness ratio (ω_n * τ): ', stiffness_ratio
    
    if (zeta < 1.0_dp) then
        write(6, '(A)') '           System is UNDERDAMPED (oscillatory)'
        if (Q_factor > 10.0_dp) then
            write(6, '(A)') '           High Q - potentially stiff, ODEX is appropriate'
        endif
    elseif (abs(zeta - 1.0_dp) < 1.0e-6_dp) then
        write(6, '(A)') '           System is CRITICALLY DAMPED'
    else
        write(6, '(A)') '           System is OVERDAMPED'
    endif
    
    if (stiffness_ratio > 100.0_dp) then
        write(6, '(A)') '           STIFF SYSTEM detected - ODEX is well-suited'
    elseif (stiffness_ratio > 10.0_dp) then
        write(6, '(A)') '           Moderately stiff system - ODEX advantageous'
    else
        write(6, '(A)') '           Mildly stiff system'
    endif
    write(6, '(A)') ''
    
    ! Initial conditions (circuit starts at rest)
    t = 0.0_dp              ! Initial time
    y(1) = 0.0_dp           ! V (initial voltage across capacitor)
    y(2) = 0.0_dp           ! i_L (initial current through inductor)
    
    ! Integration time (enough time constants to see settling)
    !tend = 10.0_dp * tau    ! Integrate for 10 time constants
    tend = 2.0_dp           ! End time. Integrate for 2 seconds
    if (tend < 1.0e-4_dp) tend = 1.0e-4_dp  ! Minimum integration time for high frequency
    !if (tend > 1.0_dp) tend = 10.0_dp        ! Maximum integration time for display
    
    ! Tolerance settings (tighter for potentially stiff system)
    rtol = 1.0d-8           ! Relative tolerance
    atol = 1.0d-10          ! Absolute tolerance
    
    ! Initial step size (conservative for stiff system)
    h = tau / 1000.0_dp     ! Start with τ/1000
    if (h > 1.0e-6_dp) h = 1.0e-6_dp      ! Cap maximum initial step
    if (h < 1.0e-12_dp) h = 1.0e-12_dp    ! Minimum reasonable step
    
    ! User integer parameter (not used)
    ipar = 0
    
    ! Initialize work arrays
    do i = 1, min(20, lwork)
        work(i) = 0.0_dp
    end do
    
    do i = 1, min(30, liwork)
        iwork(i) = 0
    end do
    
    ! Set dense output for both components
    iwork(8) = nrdens       ! Number of components for dense output
    iwork(21) = 1          ! Component 1 (V) for dense output
    iwork(22) = 2          ! Component 2 (i_L) for dense output
    
    write(6, '(A)') 'Problem setup completed.'
    write(6, '(A,2F12.6)') 'Initial conditions [V, i_L]: ', y(1), y(2)
    write(6, '(A,ES10.3)') 'Integration from t = 0 to t = ', tend
    write(6, '(A,ES10.3,A,ES10.3)') 'Tolerances: RTOL = ', rtol, ', ATOL = ', atol
    write(6, '(A,ES10.3)') 'Initial step size: ', h
    write(6, '(A)') ''
    
end subroutine setup_rlc_problem

! Main ODEX integration routine
subroutine run_odex_integration(f, solout, n, y, t, tend, h, rtol, atol, itol, &
                               iout, work, lwork, iwork, liwork, rpar, ipar, idid)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    external f, solout
    integer, intent(in) :: n, itol, iout, lwork, liwork, ipar
    real(kind=dp), intent(inout) :: t, h
    real(kind=dp), intent(in) :: tend, rtol, atol
    real(kind=dp), dimension(n), intent(inout) :: y
    real(kind=dp), dimension(lwork), intent(inout) :: work
    integer, dimension(liwork), intent(inout) :: iwork
    real(kind=dp), dimension(4), intent(inout) :: rpar
    integer, intent(out) :: idid
    
    write(6, '(A)') 'Starting ODEX integration for RLC Low-Pass Filter...'
    write(6, '(A)') '     Time (s)      Voltage (V)    Current (A)     NSTEP'
    write(6, '(A)') '--------------------------------------------------------'
    
    ! Call ODEX solver
    call odex(n, f, t, y, tend, h, &
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
subroutine print_odex_statistics(n, t, y, rtol, liwork, iwork, idid, rpar)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: n, liwork
    real(kind=dp), intent(in) :: t, rtol
    real(kind=dp), dimension(n), intent(in) :: y
    integer, dimension(liwork), intent(in) :: iwork
    integer, intent(in) :: idid
    real(kind=dp), dimension(4), intent(in) :: rpar
    real(kind=dp) :: settling_error, final_voltage, expected_voltage
    real(kind=dp) :: power_dissipated, energy_stored_L, energy_stored_C, total_energy
    real(kind=dp) :: final_current, steady_state_current, tau, omega_n
    
    write(6, '(/A)') 'Final Integration Statistics:'
    write(6, '(A)') '============================'
    write(6, '(A,ES12.5,A)') 'Final time: ', t, ' s'
    write(6, '(A,F12.6,A)') 'Final voltage: ', y(1), ' V'
    write(6, '(A,F12.6,A)') 'Final current: ', y(2), ' A'
    
    write(6, '(A,ES10.3)') 'Relative tolerance used: ', rtol
    write(6, '(A,I0)') 'Integration status (IDID): ', idid
    write(6, '(/A)') 'Detailed statistics:'
    write(6, '(A,I6)') ' Function evaluations: ', iwork(17)
    write(6, '(A,I6)') ' Integration steps:    ', iwork(18)
    write(6, '(A,I6)') ' Accepted steps:       ', iwork(19)
    write(6, '(A,I6)') ' Rejected steps:       ', iwork(20)
    
    if (iwork(18) > 0) then
        write(6, '(A,ES10.3)') ' Average step size:', t / real(iwork(18), dp)
        write(6, '(A,F8.2,A)') ' Acceptance ratio: ', &
                               100.0_dp * real(iwork(19), dp) / real(iwork(18), dp), '%'
    endif
    
    ! Circuit analysis
    write(6, '(/A)') 'Circuit Analysis:'
    write(6, '(A)') '=================='
    
    ! Calculate circuit characteristics
    tau = rpar(2) / rpar(3)  ! L/R time constant
    omega_n = 1.0_dp / sqrt(rpar(2) * rpar(4))
    
    ! Expected steady-state values (should equal V_b for low-pass at DC)
    expected_voltage = rpar(1)  ! V_b
    steady_state_current = 0.0_dp  ! No DC current in steady state
    final_voltage = y(1)
    final_current = y(2)
    
    settling_error = abs(final_voltage - expected_voltage) / expected_voltage * 100.0_dp
    
    write(6, '(A,F8.3,A)') ' Expected steady-state voltage:', expected_voltage, ' V'
    write(6, '(A,F8.3,A)') ' Actual final voltage:         ', final_voltage, ' V'
    write(6, '(A,F6.2,A)') ' Voltage settling error:        ', settling_error, '%'
    write(6, '(A,F8.6,A)') ' Final current (should → 0):   ', final_current, ' A'
    
    ! Power and energy analysis
    power_dissipated = y(1)**2 / rpar(3)  ! P = V²/R
    energy_stored_L = 0.5_dp * rpar(2) * y(2)**2  ! E_L = ½Li²
    energy_stored_C = 0.5_dp * rpar(4) * y(1)**2  ! E_C = ½CV²
    total_energy = energy_stored_L + energy_stored_C
    
    write(6, '(A,F10.6,A)')  ' Power dissipated in R: ', power_dissipated, ' W'
    write(6, '(A,ES14.3,A)') ' Energy stored in L: ', energy_stored_L, ' J'
    write(6, '(A,ES14.3,A)') ' Energy stored in C: ', energy_stored_C, ' J'
    write(6, '(A,ES13.3,A)') ' Total stored energy: ', total_energy, ' J'
    
    ! System behavior assessment
    write(6, '(/A)') 'Settling Analysis:'
    write(6, '(A)') '=================='
    if (settling_error < 1.0_dp) then
        write(6, '(A)') ' Voltage has reached steady state (< 1% error).'
    elseif (settling_error < 5.0_dp) then
        write(6, '(A)') ' Voltage is close to steady state (< 5% error).'
    else
        write(6, '(A)') ' System may need more time to settle.'
    endif
    
    if (abs(final_current) < 1.0e-6_dp) then
        write(6, '(A)') ' Current has effectively reached zero (< 1 μA).'
    elseif (abs(final_current) < 1.0e-3_dp) then
        write(6, '(A)') ' Current is small but not yet zero (< 1 mA).'
    else
        write(6, '(A)') ' Current is still significant - system settling.'
    endif
    
    ! Efficiency assessment for ODEX
    write(6, '(/A)') 'ODEX Performance:'
    write(6, '(A)') '================='
    if (iwork(20) > 0) then
        write(6, '(A,I0,A)') ' Rejected steps: ', iwork(20), ' (indicates stiffness handled)'
    else
        write(6, '(A)') ' No rejected steps - very smooth integration.'
    endif
    
    if (iwork(18) < 1000) then
        write(6, '(A)') ' Efficient integration - suitable for stiff system.'
    else
        write(6, '(A)') ' Many steps required - system may be challenging.'
    endif
    
end subroutine print_odex_statistics

! Output routine for continuous solution monitoring
subroutine solout_rlc(nr, told, t, y, n, con, ncon, icomp, nd, rpar, ipar, irtrn)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, intent(in) :: nr, n, ncon, nd, ipar
    integer, dimension(nd), intent(in) :: icomp
    real(kind=dp), intent(in) :: told
    real(kind=dp), intent(in) :: t
    real(kind=dp), dimension(n), intent(in) :: y
    real(kind=dp), dimension(ncon), intent(in) :: con
    real(kind=dp), dimension(4), intent(in) :: rpar
    integer, intent(inout) :: irtrn
    
    ! External function for dense output
    real(kind=dp), external :: contex
    
    ! Storage for output control
    real(kind=dp), save :: tout = 0.0_dp
    real(kind=dp), save :: dt_out = 0.0_dp
    logical, save :: first_call = .true.
    logical, parameter ::  iprint_into_termina = .true. ! TERMINAL PRINTOUT OPTION
    integer, parameter :: file_unit = 10
    real(kind=dp) :: V_current, iL_current, iR_current, iC_current
    real(kind=dp) :: power_R, power_source, power_L, power_C, tau
    
    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) told, ipar, irtrn
    endif
    
    if (first_call) then
        ! Calculate time constant for output interval
        tau = rpar(2) / rpar(3)  ! L/R
        
        ! Open output file for plotting data
        open(unit=file_unit, file='rlc_odex_data.csv', status='replace')
        write(file_unit, '(A)') 'Time,Voltage,Current_L,Current_R,Current_C,Power_R,Power_Source,Power_L,Power_C,Energy_L,Energy_C'
        
        ! Print initial condition to screen
        write(6, '(1X,ES12.5,2F14.6,I8)') t, y(1), y(2), nr-1
        
        ! Calculate derived quantities at initial time
        call calculate_circuit_quantities(y, rpar, V_current, iL_current, iR_current, &
                                         iC_current, power_R, power_source, power_L, power_C)
        
        ! Write initial condition to file
        write(file_unit, '(11(ES16.8,","))') t, V_current, iL_current, iR_current, iC_current, &
                                            power_R, power_source, power_L, power_C, &
                                            0.5_dp*rpar(2)*iL_current**2, 0.5_dp*rpar(4)*V_current**2
        
        ! Set output interval for high-resolution plotting
        !dt_out = tau / 100.0_dp  ! τ/100 for very high resolution
        dt_out = 0.001  ! DT - time stepping size
        !if (dt_out > 1.0e-4_dp) dt_out = 1.0e-4_dp  ! Cap for slow circuits
        !if (dt_out < 1.0e-9_dp) dt_out = 1.0e-9_dp  ! Floor for fast circuits
        
        tout = dt_out
        first_call = .false.
    else
        ! Print at regular intervals using dense output
        do while (t >= tout)
            V_current = contex(1, tout, con, ncon, icomp, nd)
            iL_current = contex(2, tout, con, ncon, icomp, nd)
            
            call calculate_circuit_quantities_from_state(V_current, iL_current, rpar, &
                                                        iR_current, iC_current, power_R, &
                                                        power_source, power_L, power_C)
            if (iprint_into_termina .eqv. .true.) then
                write(6, '(1X,ES12.5,2F14.6,I8)') tout, V_current, iL_current, nr-1
            end if

            write(file_unit, '(11(ES16.8,","))') tout, V_current, iL_current, iR_current, iC_current, &
                                                power_R, power_source, power_L, power_C, &
                                                0.5_dp*rpar(2)*iL_current**2, 0.5_dp*rpar(4)*V_current**2
            tout = tout + dt_out
        end do
    endif

end subroutine solout_rlc

! Helper subroutine to calculate circuit quantities
subroutine calculate_circuit_quantities(y, rpar, V, iL, iR, iC, power_R, power_source, power_L, power_C)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(kind=dp), dimension(2), intent(in) :: y
    real(kind=dp), dimension(4), intent(in) :: rpar
    real(kind=dp), intent(out) :: V, iL, iR, iC, power_R, power_source, power_L, power_C
    
    V = y(1)
    iL = y(2)
    iR = V / rpar(3)
    iC = iL - iR
    power_R = V * iR
    power_source = rpar(1) * iL
    power_L = V * iL  ! Power delivered to inductor
    power_C = V * iC  ! Power delivered to capacitor
    
end subroutine calculate_circuit_quantities

! Helper subroutine for quantities from interpolated state
subroutine calculate_circuit_quantities_from_state(V, iL, rpar, iR, iC, power_R, power_source, power_L, power_C)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(kind=dp), intent(in) :: V, iL
    real(kind=dp), dimension(4), intent(in) :: rpar
    real(kind=dp), intent(out) :: iR, iC, power_R, power_source, power_L, power_C
    
    iR = V / rpar(3)
    iC = iL - iR
    power_R = V * iR
    power_source = rpar(1) * iL
    power_L = V * iL
    power_C = V * iC
    
end subroutine calculate_circuit_quantities_from_state

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
