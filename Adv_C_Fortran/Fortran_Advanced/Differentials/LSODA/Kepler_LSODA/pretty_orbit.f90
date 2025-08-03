! Compile with:
! gfortran pretty_orbit.f90 YBER_ODEPACK.f ODEPACK_MODULES.f90 -o prop -lopenblas -std=legac

program dlsoda_orbital
    use iso_c_binding
    use odepack_interface
    use odepack_common
    implicit none
    
    ! External subroutines
    external orbital_equations, jac_dummy_silent
    
    ! Constants and parameters
    !integer, parameter :: dp = kind(0.0d0)
    integer, parameter :: neq = 4
    !integer, parameter :: max_steps = 10000
    !integer, parameter :: output_interval = 200
    !integer, parameter :: jt_method = 2  ! Internal generated Jacobian (kept for reference)
    
    ! Variable declarations
    double precision, dimension(neq) :: atol, y
    integer :: iopt, istate, itask, itol, jt, liw, lrw
    integer, allocatable :: iwork(:)
    double precision :: rtol, t, tout, dt
    double precision, allocatable :: rwork(:)
    type(odepack_common_data), target :: common_data
    
    call print_title()

    ! Calculate work array dimensions and allocate
    call calculate_work_dimensions(neq, lrw, liw)
    allocate(rwork(lrw), iwork(liw))
    
    ! Initialize the simulation
    call initialize_solver(iwork, rwork, common_data, rtol, atol, itol, itask, &
                          istate, iopt, lrw, liw, jt)
    
    call setup_initial_conditions(y, t, dt)
        
    ! Open output file and write initial state
    call open_output_file()
    call output_state(t, y, .true.)
    
    ! Main integration loop
    call integrate_orbit(y, t, tout, dt, iwork, rwork, common_data, rtol, atol, &
                        itol, itask, istate, iopt, lrw, liw, jt)
    
    ! Cleanup and final output
    close(11)
    call print_statistics(liw, iwork, t)
    
    ! Deallocate work arrays
    deallocate(rwork, iwork)
    
end program dlsoda_orbital

!-----------------------------------------------------------------------
! Print simulation header version 1
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
    write(*,'(a)') "      Orbital Mechanics Simulator"
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
    integer, parameter :: neq = 4
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
! Set up initial orbital conditions
!-----------------------------------------------------------------------
subroutine setup_initial_conditions(y, t, dt)
    implicit none
    integer, parameter :: neq = 4
    double precision, dimension(neq), intent(out) :: y
    double precision, intent(out) :: t, dt
    
    double precision :: earth_radius, altitude, initial_radius
    double precision :: v_0, theta
    
    ! Orbital parameters
    earth_radius = 6371.0D0
    altitude = 200.0D0
    initial_radius = earth_radius + altitude
    v_0 = 7.8D0
    theta = 90.7D0
    
    ! Initial state vector: [x, y, vx, vy]
    y(1) = initial_radius
    y(2) = 0.0D0
    y(3) = v_0 * cosd(theta)
    y(4) = v_0 * sind(theta)
    
    t = 0.0D0
    dt = 1.0D0

    ! Comment out if uneccessary
    call print_header(altitude, v_0, theta)

end subroutine setup_initial_conditions

!-----------------------------------------------------------------------
! Print simulation header
!-----------------------------------------------------------------------
subroutine print_header(altitude, v0, theta)
    implicit none
    double precision, intent(in) :: altitude, v0, theta
    
    write(*,'(a)') "ORBITAL MECHANICS USING DLSODA"
    write(*,'(a)') "Object orbiting Earth in elliptical orbit"
    write(*,'(a,f6.1,a)') "Initial altitude: ", altitude, " km"
    write(*,'(a,f4.1,a)') "Initial velocity: ", v0, " km/s"
    write(*,'(a,f4.1,a)') "Launch angle: ", theta, " degrees"
    write(*,'(a)') ""
    write(*,'(a)') "   TIME(s)      X(km)       Y(km)       VX(km/s)    VY(km/s)     R(km)"
    write(*,'(a)') "------------------------------------------------------------------------"
end subroutine print_header

!-----------------------------------------------------------------------
! Open output file for orbit data
!-----------------------------------------------------------------------
subroutine open_output_file()
    implicit none
    
    open(unit=11, file='orbit.dat')
end subroutine open_output_file

!-----------------------------------------------------------------------
! Output current state to screen and file
!-----------------------------------------------------------------------
subroutine output_state(t, y, force_output)
    implicit none
    integer, parameter :: neq = 4
    double precision, intent(in) :: t
    double precision, dimension(neq), intent(in) :: y
    logical, intent(in) :: force_output
    double precision :: radius
    
    radius = sqrt(y(1)**2 + y(2)**2)
    
    ! Always write to file
    write(11, '(6g16.8)') t, y(1), y(2), y(3), y(4), radius
    
    ! Write to screen if forced or condition met
    if (force_output) then
        write(*, '(f8.1, 4f12.3, f10.1)') t, y(1), y(2), y(3), y(4), radius
    endif
end subroutine output_state

!-----------------------------------------------------------------------
! Check if orbit should terminate
!-----------------------------------------------------------------------
logical function should_terminate(y, t)
    implicit none
    integer, parameter :: neq = 4
    double precision, dimension(neq), intent(in) :: y
    double precision, intent(in) :: t
    double precision, parameter :: earth_radius = 6371.0D0
    
    should_terminate = .false.
    
    if (sqrt(y(1)**2 + y(2)**2) < earth_radius) then
        write(*,*) 'Object crashed into Earth at t =', t
        should_terminate = .true.
    endif
end function should_terminate

!-----------------------------------------------------------------------
! Main integration loop
!-----------------------------------------------------------------------
subroutine integrate_orbit(y, t, tout, dt, iwork, rwork, common_data, &
                          rtol, atol, itol, itask, istate, iopt, lrw, liw, jt)
    use odepack_interface
    use odepack_common
    implicit none
    
    integer, parameter :: neq = 4
    integer, parameter :: max_steps = 10000
    integer, parameter :: output_interval = 1000
    integer, intent(in) :: lrw, liw
    
    double precision, dimension(neq), intent(inout) :: y
    double precision, intent(inout) :: t, tout
    double precision, intent(in) :: dt, rtol
    double precision, dimension(neq), intent(in) :: atol
    integer, dimension(liw), intent(inout) :: iwork
    double precision, dimension(lrw), intent(inout) :: rwork
    type(odepack_common_data), intent(inout) :: common_data
    integer, intent(in) :: itol, itask, iopt, jt
    integer, intent(inout) :: istate
    
    external orbital_equations, jac_dummy_silent
    logical :: should_terminate
    integer :: i_step
    
    do i_step = 1, max_steps
        tout = t + dt
        
        call dlsoda(orbital_equations, neq, y, t, tout, itol, rtol, atol, &
                   itask, istate, iopt, rwork, lrw, iwork, liw, &
                   jac_dummy_silent, jt, common_data)
        
        if (istate < 0) then
            write(*,*) 'Integration failed with istate =', istate
            write(*,*) 'Last successful time:', rwork(15)
            exit
        endif
        
        ! Output periodically to screen
        if (mod(i_step, output_interval) == 0) then
            call output_state(t, y, .true.)
        else
            call output_state(t, y, .false.)
        endif
        
        ! Check termination conditions
        if (should_terminate(y, t)) exit
    end do
end subroutine integrate_orbit

!-----------------------------------------------------------------------
! Print final statistics and optional outputs
!-----------------------------------------------------------------------
subroutine print_statistics(liw, iwork, t)
    implicit none
    integer, intent(in) :: liw
    integer, intent(in) :: iwork(liw)
    double precision, intent(in) :: t
    
    write(*,'(a)') ""
    write(*,'(a)') "=== INTEGRATION STATISTICS ==="
    write(*,'(a,i6)') 'Total integration steps (NST):', iwork(11)
    write(*,'(a,i6)') 'Function evaluations (NFE):   ', iwork(12)
    write(*,'(a,i6)') 'Jacobian evaluations (NJE):   ', iwork(13)
    write(*,'(a,i6)') 'Method order last used (NQU): ', iwork(14)
    write(*,'(a,i6)') 'Length of RWORK required:     ', iwork(17)
    write(*,'(a,i6)') 'Length of IWORK required:     ', iwork(18)
    write(*,'(a,f8.1,a)') 'Final time: ', t, ' seconds'
    write(*,'(a,f8.1,a)') 'Final time: ', t/60.0, ' minutes'
    if (t > 3600.0) then
        write(*,'(a,f9.2,a)') 'Final time: ', t/3600.0, ' hours'
    endif
end subroutine print_statistics

!-----------------------------------------------------------------------
! Orbital equations of motion
!-----------------------------------------------------------------------
subroutine orbital_equations(neq, t, y, ydot, common_data)
    use odepack_common
    implicit none
    integer, intent(in) :: neq
    double precision, intent(in) :: t
    double precision, dimension(neq), intent(in) :: y
    double precision, dimension(neq), intent(out) :: ydot
    type(odepack_common_data), intent(inout) :: common_data
    
    double precision, parameter :: GM = 3.986004418D5  ! km³/s²
    double precision :: r, r3
    
    ! Calculate distance and its cube
    r = sqrt(y(1)**2 + y(2)**2)
    r3 = r**3
    
    ! Equations of motion
    ydot(1) = y(3)                    ! dx/dt = vx
    ydot(2) = y(4)                    ! dy/dt = vy
    ydot(3) = -GM * y(1) / r3         ! dvx/dt = -GM*x/r³
    ydot(4) = -GM * y(2) / r3         ! dvy/dt = -GM*y/r³
    
    common_data%ierr = 0
    
    ! Avoid unused variable warning
    if (t /= t) write(*,*) 'Time variable is NaN'
end subroutine orbital_equations

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