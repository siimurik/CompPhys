! Compile with:
! gfortran dayum_lsoda.f90 YBER_ODEPACK.f ODEPACK_MODULES.f90 -o dam_lsoda -lopenblas -std=legacy

program dlsoda_damped_spring
    use iso_c_binding
    use odepack_interface
    use odepack_common
    implicit none
    
    ! External subroutines
    external spring_equations, jac_dummy_silent
    
    ! Constants and parameters
    integer, parameter :: neq = 2
    integer, parameter :: nvars = 3
    
    ! Variable declarations
    double precision, dimension(neq) :: y
    double precision, dimension(nvars) :: rpar  ! [mass, damping, stiffness]
    double precision :: dt, tstart, tstop
    
    ! Physical parameters
    double precision :: mass, damping, stiffness
    
    call print_title()

    ! Initialize parameters
    rpar = 0.0D0

    ! Physical parameters
    mass = 1.0d0      ! Mass (kg)
    damping = 0.5d0   ! Damping coefficient (Ns/m) 
    stiffness = 4.0d0 ! Spring constant (N/m)
    
    ! Pack parameters
    rpar(1) = mass
    rpar(2) = damping  
    rpar(3) = stiffness

    ! Initial conditions
    y(1) = 1.0D0    ! Initial position (m)
    y(2) = 0.0D0    ! Initial velocity (m/s)
    tstart = 0.0D0
    tstop  = 20.0D0
    dt = 1.0d0

    call run_lsoda_integration(neq, nvars, y, tstart, tstop, dt, rpar)
    
end program dlsoda_damped_spring

subroutine run_lsoda_integration(neq, nvars, y, tstart, tstop, dt, rpar)
    use iso_c_binding
    use odepack_interface
    use odepack_common
    implicit none

    ! Variable declarations
    integer, intent(in) :: neq, nvars
    double precision, dimension(neq), intent(inout) :: y
    double precision, dimension(neq) :: atol
    double precision, dimension(nvars), intent(in) :: rpar  ! [mass, damping, stiffness]
    integer :: iopt, istate, itask, itol, jt, liw, lrw, ipar
    integer, allocatable :: iwork(:)
    double precision, intent(inout) :: tstart, tstop, dt
    double precision :: rtol, t, tout
    double precision, allocatable :: rwork(:)
    type(odepack_common_data), target :: common_data

    ! Initialize parameters
    ipar = 0

    ! Calculate work array dimensions and allocate
    call calculate_work_dimensions(neq, lrw, liw)
    allocate(rwork(lrw), iwork(liw))
    
    ! Initialize the solver
    call initialize_solver(iwork, rwork, common_data, rtol, atol, itol, itask, &
                          istate, iopt, lrw, liw, jt)
    
    call setup_initial_conditions(neq, y, t, dt, tstart, tstop, nvars, rpar)
        
    ! Open output file and write initial state
    call open_output_file()
    call output_state(t, y, rpar, .true.)
    
    ! Main integration loop
    call integrate_lsoda(y, t, tout, dt, tstart, tstop, iwork, rwork, common_data, &
                         rtol, atol, itol, itask, istate, iopt, lrw, liw, jt)
    
    ! Cleanup and final output
    close(11)
    call print_statistics(liw, iwork, t)
    
    ! Deallocate work arrays
    deallocate(rwork, iwork)

end subroutine run_lsoda_integration

!-----------------------------------------------------------------------
! Print simulation header
!-----------------------------------------------------------------------
subroutine print_title()
    implicit none
    write(*,'(a)') "=================================================="
    write(*,'(a)') "    ____  _    ____  _____  ____  _____  "
    write(*,'(a)') "   |  _ \| |  / ___|/  _  \|  _ \|  ___| "
    write(*,'(a)') "   | | | | |  \___ \| | | || | | | |___  "
    write(*,'(a)') "   | |_| | |_____ ) | |_| || |_| | |___  "
    write(*,'(a)') "   |____/|____|____/\_____/|____/|_____| "
    write(*,'(a)') ""
    write(*,'(a)') "     Damped Spring-Mass System Solver"
    write(*,'(a)') "=================================================="
end subroutine print_title

!-----------------------------------------------------------------------
! Calculate proper work array dimensions for DLSODA
!-----------------------------------------------------------------------
subroutine calculate_work_dimensions(neq, lrw, liw)
    implicit none
    integer, intent(in) :: neq
    integer, intent(out) :: lrw, liw
    
    ! DLSODA automatically switches between stiff and non-stiff methods
    ! We need to size for the larger requirement
    
    ! For DLSODA, we need the maximum of:
    ! - Non-stiff: 20 + 16*NEQ  
    ! - Stiff: 22 + 9*NEQ + NEQ²
    ! The stiff requirement is larger for NEQ > 7
    lrw = max(20 + 16*neq, 22 + 9*neq + neq*neq)
    
    ! For DLSODA, LIW needs to handle both methods
    liw = 20 + neq
    
    write(*,'(a,i0)') 'Calculated LRW: ', lrw
    write(*,'(a,i0)') 'Calculated LIW: ', liw
end subroutine calculate_work_dimensions

!-----------------------------------------------------------------------
! Initialize solver parameters and arrays
!-----------------------------------------------------------------------
subroutine initialize_solver(iwork, rwork, common_data, rtol, atol, itol, &
                            itask, istate, iopt, lrw, liw, jt)
    use odepack_common
    implicit none
    integer, parameter :: neq = 2
    integer, intent(in) :: lrw, liw
    integer, dimension(liw), intent(out) :: iwork
    double precision, dimension(lrw), intent(out) :: rwork
    type(odepack_common_data), intent(out) :: common_data
    double precision, intent(out) :: rtol
    double precision, dimension(neq), intent(out) :: atol
    integer, intent(out) :: itol, itask, istate, iopt, jt
    
    ! Initialize arrays
    iwork = 0
    rwork = 0.0D0
    common_data%ierr = 0
    
    ! Set tolerances
    itol = 1                     ! Scalar tolerances
    rtol = 1.0d-8               ! Relative tolerance
    atol = 1.0d-10              ! Absolute tolerance (array but same value)
    
    ! Set solver options
    itask = 1                   ! Normal computation to TOUT
    istate = 1                  ! First call
    iopt = 0                    ! No optional inputs (start simple)
    jt = 2                      ! Internally generated full Jacobian
    
    ! If you want to use optional inputs, set iopt = 1 and uncomment below:
    ! iopt = 1
    ! Optional inputs (when IOPT = 1)
    ! RWORK optional inputs:
    ! rwork(5) = 0.0D0          ! H0: Initial step size (0 = let solver choose)
    ! rwork(6) = 0.0D0          ! HMAX: Maximum step size (0 = infinite)  
    ! rwork(7) = 0.0D0          ! HMIN: Minimum step size (0 = no minimum)
    ! IWORK optional inputs:
    ! iwork(5) = 0              ! MAXORD: Maximum order (0 = use default)
    ! iwork(6) = 500            ! MXSTEP: Maximum steps per call
    ! iwork(7) = 10             ! MXHNIL: Max warnings for T+H=T
end subroutine initialize_solver

!-----------------------------------------------------------------------
! Set up initial spring-mass conditions
!-----------------------------------------------------------------------
subroutine setup_initial_conditions(neq, y, t, dt, tstart, tstop, nvars, rpar)
    implicit none
    integer, intent(in) :: neq, nvars
    double precision, dimension(neq), intent(out) :: y
    double precision, dimension(nvars), intent(in) :: rpar
    double precision, intent(out) :: t, dt, tstart, tstop
    double precision :: mass, damping, stiffness
    
    mass      = rpar(1)
    damping   = rpar(2)
    stiffness = rpar(3)

    ! Initial conditions
    y(1) = 1.0D0    ! Initial position (m)
    y(2) = 0.0D0    ! Initial velocity (m/s)
    
    tstart = 0.0D0
    tstop = 20.0D0
    t = tstart
    dt = 1.0d0

    ! Print simulation parameters
    call print_header(mass, damping, stiffness)

end subroutine setup_initial_conditions

!-----------------------------------------------------------------------
! Print simulation header
!-----------------------------------------------------------------------
subroutine print_header(mass, damping, stiffness)
    implicit none
    double precision, intent(in) :: mass, damping, stiffness
    
    write(*,'(a)') "DAMPED SPRING-MASS SYSTEM USING DLSODA"
    write(*,'(a)') "Second-order differential equation solver"
    write(*,'(a,f6.3,a)') "Mass: ", mass, " kg"
    write(*,'(a,f6.3,a)') "Damping coefficient: ", damping, " Ns/m"
    write(*,'(a,f6.3,a)') "Spring constant: ", stiffness, " N/m"
    write(*,'(a)') ""
    write(*,'(a)') "    Time         Y(1)           Y(2)        Energy"
    write(*,'(a)') "------------------------------------------------------"
end subroutine print_header

!-----------------------------------------------------------------------
! Open output file for spring data
!-----------------------------------------------------------------------
subroutine open_output_file()
    implicit none
    
    open(unit=11, file='spring.dat')
end subroutine open_output_file

!-----------------------------------------------------------------------
! Output current state to screen and file
!-----------------------------------------------------------------------
subroutine output_state(t, y, rpar, force_output)
    implicit none
    integer, parameter :: neq = 2
    double precision, intent(in) :: t
    double precision, dimension(neq), intent(in) :: y
    double precision, dimension(3), intent(in) :: rpar
    logical, intent(in) :: force_output
    double precision :: kinetic_energy, potential_energy, total_energy
    
    ! Calculate energy components
    kinetic_energy = 0.5d0 * rpar(1) * y(2)**2
    potential_energy = 0.5d0 * rpar(3) * y(1)**2
    total_energy = kinetic_energy + potential_energy
    
    ! Always write to file
    write(11, '(4g16.8)') t, y(1), y(2), total_energy
    
    ! Write to screen if forced or condition met
    if (force_output) then
        write(*, '(1x,f8.4,3d15.6)') t, y(1), y(2), total_energy
    endif
end subroutine output_state

!-----------------------------------------------------------------------
! Check if simulation should terminate
!-----------------------------------------------------------------------
logical function should_terminate(y, t)
    implicit none
    integer, parameter :: neq = 2
    double precision, dimension(neq), intent(in) :: y
    double precision, intent(in) :: t
    
    should_terminate = .false.

    ! Silence unused variable message
    if (.false.) then 
        write (*,*) t
    end if
    
    ! Stop if solution has essentially decayed to zero
    if (abs(y(1)) < 1.0d-6 .and. abs(y(2)) < 1.0d-6) then
        write(*,'(a)') 'Solution has decayed to negligible values - stopping integration.'
        should_terminate = .true.
    endif
end function should_terminate

!-----------------------------------------------------------------------
! Main integration loop
!-----------------------------------------------------------------------
subroutine integrate_lsoda(y, t, tout, dt, tstart, tstop, iwork, rwork, &
                           common_data, rtol, atol, itol, itask, istate, &
                           iopt, lrw, liw, jt)
    use odepack_interface
    use odepack_common
    implicit none
    
    integer, parameter :: neq = 2
    integer, parameter :: output_interval = 1
    integer, intent(in) :: lrw, liw
    
    double precision, dimension(neq), intent(inout) :: y
    double precision, intent(inout) :: t, tout
    double precision, intent(in) :: dt, tstart, tstop, rtol
    double precision, dimension(neq), intent(in) :: atol
    integer, dimension(liw), intent(inout) :: iwork
    double precision, dimension(lrw), intent(inout) :: rwork
    type(odepack_common_data), intent(inout) :: common_data
    integer, intent(in) :: itol, itask, iopt, jt
    integer, intent(inout) :: istate
    
    ! Physical parameters for output (hardcoded to match equations)
    double precision, parameter :: mass = 1.0d0
    double precision, parameter :: damping = 0.5d0
    double precision, parameter :: stiffness = 4.0d0
    double precision, dimension(3) :: rpar_output
    
    external spring_equations, jac_dummy_silent
    logical :: should_terminate
    integer :: i_step, ntotal
    
    ! Set up parameters for output only
    rpar_output(1) = mass
    rpar_output(2) = damping
    rpar_output(3) = stiffness
    
    ! Calculate total number of steps
    ntotal = int((tstop - tstart)/dt)
    
    do i_step = 1, ntotal
        tout = t + dt
        
        call dlsoda(spring_equations, neq, y, t, tout, itol, rtol, atol, &
                   itask, istate, iopt, rwork, lrw, iwork, liw, &
                   jac_dummy_silent, jt, common_data)
        
        if (istate < 0) then
            write(*,'(//a,i3)') ' Error halt: ISTATE =', istate
            select case (istate)
            case (-1)
                write(*,'(a)') 'Excess work done on this call. Try increasing MXSTEP or using larger tolerances.'
                write(*,'(a)') 'Current solution is very small - this may be near the numerical limit.'
            case (-2)
                write(*,'(a)') 'Excess accuracy requested. (Tolerances too small.)'
            case (-3)
                write(*,'(a)') 'Illegal input detected. (See printed message.)'
            case (-4)
                write(*,'(a)') 'Repeated error test failures. (Check all input.)'
            case (-5)
                write(*,'(a)') 'Repeated convergence failures. (Perhaps bad Jacobian or wrong JT/tolerances.)'
            case (-6)
                write(*,'(a)') 'Error weight became zero during problem. (Solution component vanished.)'
            end select
            
            ! If we've integrated far enough and solution is very small, consider it successful
            if (istate == -1 .and. t > 10.0d0 .and. abs(y(1)) < 1.0d-3 .and. abs(y(2)) < 1.0d-3) then
                write(*,'(a)') 'Solution has decayed to near zero - stopping integration.'
                exit
            else
                stop 1
            endif
        endif
        
        ! Output to screen and file
        if (mod(i_step, output_interval) == 0) then
            call output_state(t, y, rpar_output, .true.)
        else
            call output_state(t, y, rpar_output, .false.)
        endif
        
        ! Check termination conditions
        if (should_terminate(y, t)) exit
        
        ! Don't go beyond our intended final time
        if (tout > tstop) exit
    end do
    
    write(*,'(a)') 'Integration completed successfully.'
end subroutine integrate_lsoda

!-----------------------------------------------------------------------
! Print final statistics and optional outputs
!-----------------------------------------------------------------------
subroutine print_statistics(liw, iwork, t)
    implicit none
    integer, intent(in) :: liw
    integer, intent(in) :: iwork(liw)
    double precision, intent(in) :: t
    
    write(*,'(a)') ""
    write(*,'(a)') "=== FINAL INTEGRATION STATISTICS ==="
    write(*,'(a,i6)') 'Total integration steps (NST):', iwork(11)
    write(*,'(a,i6)') 'Function evaluations (NFE):   ', iwork(12)
    write(*,'(a,i6)') 'Jacobian evaluations (NJE):   ', iwork(13)
    write(*,'(a,i6)') 'Method order last used (NQU): ', iwork(14)
    write(*,'(a,i6)') 'Length of RWORK required:     ', iwork(17)
    write(*,'(a,i6)') 'Length of IWORK required:     ', iwork(18)
    write(*,'(a,f8.1,a)') 'Final time: ', t, ' seconds'
    
    ! Additional useful statistics if available
    if (liw >= 22) then
        write(*,'(a,i0)') 'No. nonlinear iterations: ', iwork(19)
        write(*,'(a,i0)') 'No. nonlinear convergence failures: ', iwork(20)
        write(*,'(a,i0)') 'No. error test failures: ', iwork(21)
    endif
end subroutine print_statistics

!-----------------------------------------------------------------------
! Spring-mass system equations of motion
!-----------------------------------------------------------------------
subroutine spring_equations(neq, t, y, ydot, common_data)
    use odepack_common
    implicit none
    integer, intent(in) :: neq
    double precision, intent(in) :: t
    double precision, dimension(neq), intent(in) :: y
    double precision, dimension(neq), intent(out) :: ydot
    type(odepack_common_data), intent(inout) :: common_data
    
    ! Local variables - hardcoded for simplicity
    double precision, parameter :: mass = 1.0d0
    double precision, parameter :: damping = 0.5d0
    double precision, parameter :: stiffness = 4.0d0
    
    ! Equations of motion
    ydot(1) = y(2)  ! dx/dt = v
    ydot(2) = -(damping/mass) * y(2) - (stiffness/mass) * y(1)  ! dv/dt = -(c/m)*v - (k/m)*x
    
    common_data%ierr = 0
    
    ! Avoid unused variable warnings
    if (t /= t) write(*,*) 'Unused variable access'
end subroutine spring_equations

!-----------------------------------------------------------------------
! Dummy Jacobian subroutine (not used when jt=2)
!-----------------------------------------------------------------------
subroutine jac_dummy_silent(neq, t, y, ml, mu, pd, nrowpd, common_data)
    use odepack_common
    implicit none
    integer, intent(in) :: neq, ml, mu, nrowpd
    double precision, intent(in) :: t, y(neq)
    double precision, intent(out) :: pd(nrowpd, neq)
    type(odepack_common_data), intent(inout) :: common_data
    
    ! Initialize outputs
    pd = 0.0d0
    common_data%ierr = 0
    
    ! Silence unused variable warnings with simple conditionals
    if (neq < 0 .or. ml < 0 .or. mu < 0 .or. nrowpd < 0) continue
    if (t < -huge(t) .or. size(y) < 0) continue
end subroutine jac_dummy_silent