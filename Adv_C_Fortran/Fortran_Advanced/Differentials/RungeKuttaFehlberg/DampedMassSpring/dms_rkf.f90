! Compile with: gfortran damspr_rkf45.f90 rkf45.f90 -o dam_rkf45
program main
    implicit none
    integer, parameter :: neqn = 2
    double precision :: mass, damping, stiffness
    double precision :: tstart, tstop, dt

    double precision, dimension(neqn) :: y

    ! Physical parameters (same as DVODE version)
    mass = 1.0d0      ! Mass (kg)
    damping = 0.5d0   ! Damping coefficient (Ns/m) 
    stiffness = 4.0d0 ! Spring constant (N/m)

    ! Initial conditions (same as DVODE version)
    y(1) = 1.0d0    ! Initial position (m)
    y(2) = 0.0d0    ! Initial velocity (m/s)
    tstart = 0.0d0  ! Start time
    tstop = 20.0d0  ! End time
    dt = 1.0d0

    ! Store parameters in common block for access by RHS function
    call set_parameters(mass, damping, stiffness)

    call solve_rkf(neqn, y, tstart, tstop, dt)

end program main

! Module to hold parameters (alternative to common block)
module parameters
    implicit none
    double precision :: param_mass, param_damping, param_stiffness
end module parameters

! Subroutine to set parameters
subroutine set_parameters(mass, damping, stiffness)
    use parameters
    implicit none
    double precision, intent(in) :: mass, damping, stiffness
    
    param_mass = mass
    param_damping = damping
    param_stiffness = stiffness
end subroutine set_parameters

subroutine solve_rkf(neqn, y, t_start, tstop, dt)
    implicit none

    integer, intent(in) :: neqn
    double precision, intent(in) :: t_start, tstop, dt
    double precision, dimension(neqn), intent(inout) :: y
    double precision, dimension(neqn) :: yp

    integer :: flag, ntotal, iout
    double precision :: relerr, abserr, tout, t

    external fex_rkf45_wrapper

    ! Initialize time
    t = t_start

    ! Error tolerances (similar to DVODE version)
    relerr = 1.0d-8   ! Relative error tolerance
    abserr = 1.0d-8  ! Absolute error tolerance

    ! Initialize flag for first call
    flag = 1

    ! Set up output parameters
    ntotal = int(tstop/dt)
    tout = t + dt

    write(6, '(A)') ''
    write(6, '(A)') 'Starting RKF45 integration...'
    write(6, '(A)') '    Time         Y(1)           Y(2)        Energy'
    write(6, '(A)') '------------------------------------------------------'

    ! Print initial condition
    call print_state(t, y)

    ! Initialize derivatives for first call
    call fex_rkf45_wrapper(t, y, yp)

    ! Main integration loop
    do iout = 1, ntotal
        ! Call RKF45 integrator (standard interface without rpar)
        call rkf45(fex_rkf45_wrapper, neqn, y, yp, t, tout, relerr, abserr, flag)
        
        ! Check for errors
        if (flag > 2) then
            call handle_rkf45_error(flag, relerr, abserr)
            if (flag == 8) then
                stop 1
            endif
            ! Reset flag and continue for recoverable errors
            flag = 2
        endif
        
        ! Print current state
        call print_state(t, y)
        
        ! Stop if solution has essentially decayed to zero
        if (abs(y(1)) < 1.0d-6 .and. abs(y(2)) < 1.0d-6) then
            write(6, '(A)') 'Solution has decayed to negligible values - stopping integration.'
            exit
        endif
        
        ! Set up next output time
        tout = tout + dt
        if (tout > tstop) exit
    enddo
    
    write(6, '(A)') 'Integration completed successfully.'

end subroutine solve_rkf

! Right-hand side function wrapper (standard RKF45 interface)
subroutine fex_rkf45_wrapper(t, y, yp)
    implicit none
    double precision, intent(in) :: t
    double precision, dimension(*), intent(in) :: y
    double precision, dimension(*), intent(out) :: yp
    
    call fex_rkf45(t, y, yp)
end subroutine fex_rkf45_wrapper

! Right-hand side function (same physics as DVODE version)
subroutine fex_rkf45(t, y, ydot)
    use parameters
    implicit none
    double precision, intent(in) :: t
    double precision, dimension(*), intent(in) :: y
    double precision, dimension(*), intent(out) :: ydot

    ! Silence unused parameter warning
    if (.false.) then
        write(*,*) t
    end if
    
    ! Same equations as DVODE version
    ydot(1) = y(2)  ! dx/dt = v
    ydot(2) = -(param_damping/param_mass) * y(2) - (param_stiffness/param_mass) * y(1)  ! dv/dt = -(c/m)*v - (k/m)*x
    
end subroutine fex_rkf45

! Subroutine to print current state with energy calculation
subroutine print_state(t, y)
    use parameters
    implicit none
    double precision, intent(in) :: t
    double precision, dimension(*), intent(in) :: y
    double precision :: kinetic_energy, potential_energy, total_energy
    
    ! Calculate energies (same as DVODE version)
    kinetic_energy = 0.5d0 * param_mass * y(2)**2
    potential_energy = 0.5d0 * param_stiffness * y(1)**2
    total_energy = kinetic_energy + potential_energy
    
    write(6, 20) t, y(1), y(2), total_energy
20  format(' ', F8.4, 3d15.6)
end subroutine print_state

! Error handling for RKF45
subroutine handle_rkf45_error(flag, relerr, abserr)
    implicit none
    integer, intent(in) :: flag
    double precision, intent(in) :: relerr, abserr
    
    write(6, '(/A,I0)') 'RKF45 returned with FLAG = ', flag
    
    select case (flag)
    case (3)
        write(6, '(A)') 'Integration not completed: RELERR too small.'
        write(6, '(A,E12.4)') 'RELERR increased to: ', relerr
        write(6, '(A)') 'Continuing integration...'
    case (4)
        write(6, '(A)') 'Integration not completed: Too many function evaluations.'
        write(6, '(A)') 'Continuing with reset function counter...'
    case (5)
        write(6, '(A)') 'Integration not completed: Solution vanished.'
        write(6, '(A,E12.4)') 'Consider using non-zero ABSERR: ', abserr
        write(6, '(A)') 'Attempting to continue...'
    case (6)
        write(6, '(A)') 'Integration not completed: Requested accuracy not achievable.'
        write(6, '(A)') 'Consider increasing error tolerances.'
        write(6, '(A)') 'This may indicate a trouble spot or singularity.'
    case (7)
        write(6, '(A)') 'Integration inefficient: Output frequency too high.'
        write(6, '(A)') 'Consider using larger output intervals.'
    case (8)
        write(6, '(A)') 'Invalid input parameters detected!'
        write(6, '(A)') 'Check: NEQN > 0, RELERR >= 0, ABSERR >= 0, valid FLAG'
    case default
        write(6, '(A)') 'Unknown error condition.'
    end select
    write(6, '(A)') ''
end subroutine handle_rkf45_error