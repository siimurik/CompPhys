! gfortran damspr_vode.f90 VODEPACK.f90 -o dam -lopenblas

program main
    implicit none
    integer, parameter :: neq = 2
    integer, parameter :: nvars = 3
    integer :: ipar
    double precision :: mass, damping, stiffness
    double precision :: tstart, tstop, dt

    double precision, dimension(neq) :: y
    double precision, dimension(nvars) :: rpar  ! [mass, damping, stiffness]

    ! Initialize parameters
    rpar = 0.0D0
    ipar = 0

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

    ! Run the main integration
    call run_dvode_integration(neq, nvars, y, tstart, tstop, dt, rpar)

end program main

! Subroutine to automatically calculate VODE array lengths
subroutine calculate_vode_arrays(neq, mf, lrw, liw)
    implicit none
    integer, intent(in) :: neq, mf
    integer, intent(out) :: lrw, liw
    integer, parameter :: ml = 5, mu = 5  ! Default bandwidths for banded case
    
    write(6, '(A)') 'Choosing method based on MF number...'
    write(6, '(A, i4)') 'MF number:', mf
    write(6, '(A)') 'Chosen method:'

    ! Calculate LRW based on method flag
    select case (mf)
    case (10)
        ! Nonstiff (Adams) method
        write(6, '(A)') 'Nonstiff (Adams) method.'
        lrw = 20 + 16 * neq
    case (21, 22)
        write(6, '(A)') 'Stiff (BDF) method with full Jacobian.'
        ! Stiff (BDF) method with full Jacobian
        lrw = 22 + 9 * neq + 2 * neq**2
    case (24, 25)
        ! Stiff (BDF) method with banded Jacobian
        write(6, '(A)') 'Stiff (BDF) method with banded Jacobian.'
        lrw = 22 + 11 * neq + (3 * ml + 2 * mu) * neq
    case default
        ! Default to full Jacobian case
        write(6, '(A)') 'Default case'
        lrw = 22 + 9 * neq + 2 * neq**2
    end select
    
    ! Calculate LIW based on method flag
    select case (mf)
    case (10)
        ! Nonstiff method
        liw = 30
    case (21, 22, 24, 25)
        ! Stiff methods
        liw = 30 + neq
    case default
        ! Default to stiff case
        liw = 30 + neq
    end select
    
    write(6, '(A)') ''
    write(6, '(A,I0,A,I0)') 'Calculated array sizes: LRW = ', lrw, ', LIW = ', liw

end subroutine calculate_vode_arrays

! Main integration subroutine
subroutine run_dvode_integration(neq, nvars, y, tstart, tstop, dt, rpar)
    implicit none
    external fex, jac
    double precision, parameter :: TOL = 1.0d-8
    integer, intent(in) :: neq, nvars
    integer :: mf, liw, lrw, ipar, ntotal
    double precision, dimension(nvars), intent(inout) :: rpar
    double precision, intent(in) :: tstart, tstop, dt
    double precision :: t, tout
    double precision, dimension(neq), intent(inout) :: y
    double precision :: atol, rtol
    double precision, allocatable, dimension(:) :: rwork
    integer, allocatable, dimension(:) :: iwork
    
    ! Local variables
    integer :: itol, itask, istate, iopt, iout
    double precision :: kinetic_energy, potential_energy, total_energy

    ! Initialize time variables
    t = tstart
    
    ! Initialize parameters
    ipar = 0

    ! 10 - Nonstiff (Adams) method
    ! 21, 22 - Stiff (BDF) method with full Jacobian
    ! 24, 25 - Stiff (BDF) method with banded Jacobian
    mf = 10   ! Choose optimal method here

    ! Allocate work arrays
    call calculate_vode_arrays(neq, mf, lrw, liw)
    allocate(rwork(lrw))
    allocate(iwork(liw))

    ! Tolerances - use different relative and absolute tolerances
    rtol = TOL  ! relative tolerance
    atol = 1.0d-12  ! absolute tolerance (smaller for better accuracy near zero)

    ! Solver options
    itol = 1      ! both rtol and atol are scalars
    itask = 1     ! normal computation
    istate = 1    ! first call
    iopt = 1      ! use optional input

    ! Set up output times
    ntotal = int((tstop - tstart)/dt)
    tout = t + dt

    ! Set optional inputs
    rwork(5) = 0.0d0    ! Initial step size (0 = automatic)
    rwork(6) = dt * 2.0d0  ! Maximum step size (reasonable value)
    rwork(7) = 0.0d0    ! Minimum step size (0 = automatic)
    iwork(5) = 0        ! Maximum order (0 = default)
    iwork(6) = 5000     ! Maximum number of steps (increased from default 500)
    iwork(7) = 1        ! Maximum number of messages

    write(6, '(A)') ''
    write(6, '(A)') 'Starting DVODE integration...'
    write(6, '(A)') '    Time         Y(1)           Y(2)        Energy'
    write(6, '(A)') '------------------------------------------------------'

    ! Print initial condition
    kinetic_energy = 0.5d0 * rpar(1) * y(2)**2
    potential_energy = 0.5d0 * rpar(3) * y(1)**2
    total_energy = kinetic_energy + potential_energy
    write(6, 20) t, y(1), y(2), total_energy

    ! Main integration loop
    do iout = 1, ntotal
        call dvode(fex, neq, y, t, tout, itol, rtol, atol, itask, istate, &
                   iopt, rwork, lrw, iwork, liw, jac, mf, rpar, ipar)
        
        ! Calculate total energy for output
        kinetic_energy = 0.5d0 * rpar(1) * y(2)**2
        potential_energy = 0.5d0 * rpar(3) * y(1)**2
        total_energy = kinetic_energy + potential_energy
        
        ! Edit terminal printing here
        write(6, 20) t, y(1), y(2), total_energy
20      format(' ', F8.4, 3d15.6)
        
        if (istate < 0) then
            write(6, 90) istate
90          format(///' Error halt: ISTATE =', i3)
            select case (istate)
            case (-1)
                write(6, '(A)') 'Excess work done on this call. Try increasing MXSTEP or using larger tolerances.'
                write(6, '(A)') 'Current solution is very small - this may be near the numerical limit.'
            case (-2)
                write(6, '(A)') 'Excess accuracy requested. (Tolerances too small.)'
            case (-3)
                write(6, '(A)') 'Illegal input detected. (See printed message.)'
            case (-4)
                write(6, '(A)') 'Repeated error test failures. (Check all input.)'
            case (-5)
                write(6, '(A)') 'Repeated convergence failures. (Perhaps bad Jacobian or wrong MF/tolerances.)'
            case (-6)
                write(6, '(A)') 'Error weight became zero during problem. (Solution component vanished.)'
            end select
            
            ! If we've integrated far enough and solution is very small, consider it successful
            if (istate == -1 .and. t > 10.0d0 .and. abs(y(1)) < 1.0d-3 .and. abs(y(2)) < 1.0d-3) then
                write(6, '(A)') 'Solution has decayed to near zero - stopping integration.'
                exit
            else
                stop 1
            endif
        endif
        
        ! Stop if solution has essentially decayed to zero
        if (abs(y(1)) < 1.0d-6 .and. abs(y(2)) < 1.0d-6) then
            write(6, '(A)') 'Solution has decayed to negligible values - stopping integration.'
            exit
        endif
        
        tout = tout + dt
        if (tout > tstop) exit  ! Don't go beyond our intended final time
    enddo
    
    write(6, '(A)') 'Integration completed successfully.'

    call print_dvode_statistics(iwork, liw)

    ! Clean up
    deallocate(rwork)
    deallocate(iwork)

end subroutine run_dvode_integration

! Right-hand side function
subroutine fex(neq, x, y, ydot, rpar, ipar)
    implicit none
    integer, intent(in) :: neq, ipar
    double precision, intent(in) :: x
    double precision, dimension(*), intent(in) :: rpar
    double precision, intent(in), dimension(neq) :: y
    double precision, intent(out), dimension(neq) :: ydot
    
    ! Local variables
    double precision :: mass, damping, stiffness

    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) x, ipar
    end if
    
    mass = rpar(1)
    damping = rpar(2)
    stiffness = rpar(3)
    
    ydot(1) = y(2)  ! dx/dt = v
    ydot(2) = -(damping/mass) * y(2) - (stiffness/mass) * y(1)  ! dv/dt = -(c/m)*v - (k/m)*x
    
end subroutine fex

! Jacobian function
subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
    implicit none
    integer, intent(in) :: neq, ml, mu, nrpd, ipar
    double precision, intent(in) :: t
    double precision, dimension(*), intent(in) :: rpar
    double precision, intent(in), dimension(neq) :: y
    double precision, intent(out), dimension(nrpd, neq) :: pd
    double precision :: mass, damping, stiffness
    
    ! Silence unused parameter warnings
    if (.false.) then
        write(*,*) t, ml, mu, y(1), ipar
    endif
    
    mass = rpar(1)
    damping = rpar(2)
    stiffness = rpar(3)
    
    ! Jacobian matrix elements
    ! ∂f₁/∂x₁ = 0, ∂f₁/∂x₂ = 1
    pd(1, 1) = 0.0d0
    pd(1, 2) = 1.0d0
    
    ! ∂f₂/∂x₁ = -k/m, ∂f₂/∂x₂ = -c/m
    pd(2, 1) = -stiffness/mass
    pd(2, 2) = -damping/mass
    
end subroutine jac

! Subroutine to print final statistics
subroutine print_dvode_statistics(iwork, liw)
    implicit none
    integer, intent(in) :: liw
    integer, dimension(liw), intent(in) :: iwork
    
    write(6, '(/A)') 'Final Integration Statistics:'
    write(6, '(A)') '============================'
    write(6, 60) iwork(11), iwork(12), iwork(13), iwork(19), &
                 iwork(20), iwork(21), iwork(22)
60  format(' No. steps = ', i4, / &
           ' No. f-s   = ', i4, / &
           ' No. J-s   = ', i4, / &
           ' No. LU-s  = ', i4, / &
           ' No. nonlinear iterations = ', i4, / &
           ' No. nonlinear convergence failures =', i3, / &
           ' No. error test failures =', i4, /)
    
    ! Additional useful statistics
    if (liw >= 22) then
        write(6, '(A,I0)') ' Method order used: ', iwork(14)
        write(6, '(A,I0)') ' Length of real work array used: ', iwork(17)
        write(6, '(A,I0)') ' Length of integer work array used: ', iwork(18)
    endif

end subroutine print_dvode_statistics