!   gfortran lorenz.f90 odepack.f odepack_interface.f90 odepack_common.f90 odepack_sub1.f odepack_sub2.f odepack_mod.f90 -o lorenz -lopenblas -std=legacy

program lorenz_attractor
    use iso_c_binding, only: c_ptr, c_f_pointer, c_loc
    use odepack_interface, only: DINTDY, DROOTS, DSTODA
    use odepack_common
    implicit none
    
    ! Problem parameters
    integer, parameter :: neq_size = 3
    type(odepack_common_data), target :: common_data
    
    ! Variables
    double precision :: y(neq_size), t, t_end
    integer :: istate, nsteps
    
    ! Initialize Lorenz system
    call initialize_lorenz(y, t, t_end)
    istate = 1
    nsteps = 1000
    
    ! Initialize common_data structure
    common_data%ierr = 0
    
    write(*,'(a)') '# Lorenz Attractor using DLSODA'
    write(*,'(a)') '# sigma=10, rho=28, beta=8/3 (chaotic regime)'
    write(*,'(a)') '# t           x           y           z'
    
    ! Integrate and output
    call differentiate_system(y, t, t_end, nsteps, istate, common_data)
    
    ! Print final statistics
    call print_statistics(common_data)
    
end program lorenz_attractor

!-----------------------------------------------------------------------
! Initialize the Lorenz system parameters and initial conditions
!-----------------------------------------------------------------------
subroutine initialize_lorenz(y, t, t_end)
    implicit none
    double precision, intent(out) :: y(3), t, t_end
    
    ! Initial conditions (slightly off equilibrium to trigger chaos)
    y(1) = 1.0d0      ! x(0) = 1
    y(2) = 1.0d0      ! y(0) = 1  
    y(3) = 1.0d0      ! z(0) = 1
    
    ! Time span
    t = 0.0d0         ! Start time
    t_end = 50.0d0    ! End time (long enough to see chaotic behavior)
    
end subroutine initialize_lorenz

!-----------------------------------------------------------------------
! Main integration routine with error handling
!-----------------------------------------------------------------------
subroutine differentiate_system(y, t, t_end, nsteps, istate, common_data)
    use odepack_common
    implicit none
    
    double precision, intent(inout) :: y(3), t
    double precision, intent(in) :: t_end
    integer, intent(in) :: nsteps
    integer, intent(inout) :: istate
    type(odepack_common_data), intent(inout) :: common_data
    
    ! DLSODA parameters
    external :: lorenz_equations, jac_dummy
    integer, parameter :: lrw = 70, liw = 23
    integer :: neq, itol, itask, iopt, jt
    double precision :: rtol, atol, rwork(lrw), tout, dt
    integer :: iwork(liw), i
    
    ! Setup DLSODA parameters
    call setup_dlsoda_params(neq, itol, rtol, atol, itask, iopt, jt)
    
    ! Initialize work arrays
    rwork = 0.0d0
    iwork = 0
    
    ! Integration parameters
    dt = t_end / dble(nsteps)
    
    ! Output initial condition
    write(*,'(f10.4, 3f12.6)') t, y(1), y(2), y(3)
    
    ! Main integration loop
    do i = 1, nsteps
        tout = t + dt
        
        ! Double Livermore Solver for Ordinary Differential Equations 
        call dlsoda(lorenz_equations, neq, y, t, tout, itol, rtol, atol, &
                   itask, istate, iopt, rwork, lrw, iwork, liw, &
                   jac_dummy, jt, common_data)
        
        if (istate < 0) then
            call handle_error(istate, iwork)
            return
        end if
        
        ! Output every 10th step to avoid too much data
        if (mod(i, 10) == 0) then
            write(*,'(f10.4, 3f12.6)') t, y(1), y(2), y(3)
        end if
        
        ! Print method switches
        if (i > 1 .and. mod(i, 100) == 0) then
            call check_method_switch(i, iwork)
        end if
    end do
    
end subroutine differentiate_system

!-----------------------------------------------------------------------
! Setup DLSODA parameters
!-----------------------------------------------------------------------
subroutine setup_dlsoda_params(neq, itol, rtol, atol, itask, iopt, jt)
    implicit none
    integer, intent(out) :: neq, itol, itask, iopt, jt
    double precision, intent(out) :: rtol, atol
    
    neq   = 3           ! Number of equations
    itol  = 1           ! Scalar tolerances
    rtol  = 1.0d-8      ! Relative tolerance (tight for accuracy)
    atol  = 1.0d-10     ! Absolute tolerance
    itask = 1           ! Normal computation
    iopt  = 0           ! No optional inputs
    jt    = 2           ! Internal Jacobian computation
    
end subroutine setup_dlsoda_params

!-----------------------------------------------------------------------
! Lorenz equations: dx/dt = sigma*(y-x), dy/dt = x*(rho-z)-y, dz/dt = x*y-beta*z
!-----------------------------------------------------------------------
subroutine lorenz_equations(neq, t, y, ydot, common_data)
    use odepack_common
    implicit none
    integer, intent(in) :: neq
    double precision, intent(in) :: t, y(neq)
    double precision, intent(out) :: ydot(neq)
    type(odepack_common_data), intent(inout) :: common_data
    
    ! Lorenz parameters (classic chaotic values)
    double precision, parameter :: sigma = 10.0d0
    double precision, parameter :: rho = 28.0d0
    double precision, parameter :: beta = 8.0d0/3.0d0
    
    ! Initialize error flag
    common_data%ierr = 0
    
    ! Check for invalid inputs
    if (any(y /= y)) then  ! Check for NaN
        common_data%ierr = -1
        return
    end if
    
    if (any(abs(y) > 1.0d10)) then  ! Check for overflow
        common_data%ierr = -2
        return
    end if
    
    ! Lorenz equations
    ydot(1) = sigma * (y(2) - y(1))             ! dx/dt
    ydot(2) = y(1) * (rho - y(3)) - y(2)        ! dy/dt  
    ydot(3) = y(1) * y(2) - beta * y(3)         ! dz/dt
    
    ! Check for invalid outputs
    if (any(ydot /= ydot)) then  ! Check for NaN in derivatives
        common_data%ierr = -3
        return
    end if
    
end subroutine lorenz_equations

!-----------------------------------------------------------------------
! Error handling routine
!-----------------------------------------------------------------------
subroutine handle_error(istate, iwork)
    implicit none
    integer, intent(in) :: istate, iwork(*)
    
    write(*,'(/a,i0)') 'Integration failed with ISTATE = ', istate
    
    select case(istate)
    case(-1)
        write(*,'(a)') 'Excess work done - possibly wrong Jacobian type'
    case(-2)
        write(*,'(a)') 'Tolerances too small for machine precision'
    case(-3)
        write(*,'(a)') 'Illegal input parameters detected'
    case(-4)
        write(*,'(a)') 'Repeated error test failures - check inputs'
    case(-5)
        write(*,'(a)') 'Convergence failures - possibly bad Jacobian'
    case(-6)
        write(*,'(a)') 'Error weight became zero'
    case(-7)
        write(*,'(a)') 'Insufficient workspace'
        write(*,'(a,i0)') ' Required RWORK: ', iwork(17)
        write(*,'(a,i0)') ' Required IWORK: ', iwork(18)
    case(-8)
        write(*,'(a)') 'Error in user function (F routine)'
        write(*,'(a)') 'Possible causes: NaN, overflow, or numerical instability'
        write(*,'(a)') 'Try smaller time steps or different initial conditions'
    case default
        write(*,'(a,i0)') 'Unknown error code: ', istate
    end select
    
end subroutine handle_error

!-----------------------------------------------------------------------
! Check for method switches and report
!-----------------------------------------------------------------------
subroutine check_method_switch(step, iwork)
    implicit none
    integer, intent(in) :: step, iwork(*)
    !double precision, intent(in) :: rwork(*)
    
    ! Report current method and statistics
    if (mod(step, 200) == 0) then
        write(*,'(a,i0,a,i0,a,i0,a,i0)') &
            '# Step ', step, ': Method=', iwork(20), &
            ', Steps=', iwork(11), ', F-evals=', iwork(12)
    end if
    
end subroutine check_method_switch

!-----------------------------------------------------------------------
! Print final integration statistics
!-----------------------------------------------------------------------
subroutine print_statistics(common_data)
    use odepack_common
    implicit none
    type(odepack_common_data), intent(in) :: common_data
    
    ! Note: In a real implementation, you would access statistics from common_data
    ! For this example, we'll provide a template
    write(*,'(/a)') '# Final Statistics:'
    write(*,'(a)') '# - Lorenz system exhibits sensitive dependence on initial conditions'
    write(*,'(a)') '# - DLSODA automatically handles the varying stiffness in chaotic regions'
    write(*,'(a)') '# - The attractor has a fractal structure with dimension ~2.06'
    write(*,'(a)') '# - Method switching helps maintain efficiency throughout trajectory'
    
end subroutine print_statistics

!-----------------------------------------------------------------------
! Dummy Jacobian (not used with JT=2)
!-----------------------------------------------------------------------
subroutine jac_dummy(neq, t, y, ml, mu, pd, nrowpd, common_data)
    use odepack_common
    implicit none
    integer, intent(in) :: neq, ml, mu, nrowpd
    double precision, intent(in) :: t, y(neq)
    double precision, intent(out) :: pd(nrowpd, neq)
    type(odepack_common_data), intent(inout) :: common_data

    ! Initialize error flag
    common_data%ierr = 0
    
    ! Dummy routine - not called when JT=2
    
end subroutine jac_dummy

!-----------------------------------------------------------------------
! Fake Use Subroutines to Silence Unused Variable Warnings
!-----------------------------------------------------------------------
! These subroutines do nothing but "use" variables to prevent compiler
! warnings about unused dummy arguments. They are designed to be optimized
! away by the compiler while satisfying the requirement that all dummy 
! arguments appear to be used.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! For real/double precision scalars
!-----------------------------------------------------------------------
subroutine r8_fake_use(x)
    implicit none
    double precision, intent(in) :: x
    
    ! Do absolutely nothing with x, but compiler thinks it's used
    if (x < -huge(x)) then
        continue  ! This will never execute but compiler sees x is "used"
    end if
end subroutine r8_fake_use

!-----------------------------------------------------------------------
! For real/double precision arrays (1D)
!-----------------------------------------------------------------------
subroutine r8_array_fake_use(x)
    implicit none
    double precision, intent(in) :: x(:)
    
    ! Use array in a way that will never execute
    if (size(x) < 0) then
        continue  ! Never executes since size() is always >= 0
    end if
end subroutine r8_array_fake_use

!-----------------------------------------------------------------------
! For real/double precision 2D arrays
!-----------------------------------------------------------------------
subroutine r8_matrix_fake_use(x)
    implicit none
    double precision, intent(in) :: x(:,:)
    
    ! Use matrix in a way that will never execute
    if (size(x,1) < 0) then
        continue  ! Never executes
    end if
end subroutine r8_matrix_fake_use

!-----------------------------------------------------------------------
! For integers
!-----------------------------------------------------------------------
subroutine int_fake_use(i)
    implicit none
    integer, intent(in) :: i
    
    ! Do nothing with i but compiler thinks it's used
    if (i < -huge(i)) then
        continue  ! Never executes
    end if
end subroutine int_fake_use

!-----------------------------------------------------------------------
! For integer arrays
!-----------------------------------------------------------------------
subroutine int_array_fake_use(i)
    implicit none
    integer, intent(in) :: i(:)
    
    ! Use array in a way that will never execute
    if (size(i) < 0) then
        continue  ! Never executes
    end if
end subroutine int_array_fake_use

!-----------------------------------------------------------------------
! For logical variables
!-----------------------------------------------------------------------
subroutine logical_fake_use(flag)
    implicit none
    logical, intent(in) :: flag
    
    ! Use logical in a way that will never execute
    if (flag .and. .not. flag) then
        continue  ! Never executes (contradiction)
    end if
end subroutine logical_fake_use

!-----------------------------------------------------------------------
! For character strings
!-----------------------------------------------------------------------
subroutine char_fake_use(str)
    implicit none
    character(*), intent(in) :: str
    
    ! Use string in a way that will never execute
    if (len(str) < 0) then
        continue  ! Never executes since len() is always >= 0
    end if
end subroutine char_fake_use

!-----------------------------------------------------------------------
! Generic subroutine for multiple variables (mixed types)
!-----------------------------------------------------------------------
subroutine multi_fake_use(r_var, i_var, l_var)
    implicit none
    double precision, intent(in), optional :: r_var
    integer, intent(in), optional :: i_var
    logical, intent(in), optional :: l_var
    
    ! Use all present variables
    if (present(r_var)) then
        if (r_var < -huge(r_var)) continue
    end if
    
    if (present(i_var)) then
        if (i_var < -huge(i_var)) continue
    end if
    
    if (present(l_var)) then
        if (l_var .and. .not. l_var) continue
    end if
end subroutine multi_fake_use

!-----------------------------------------------------------------------
! Example usage in your dummy Jacobian subroutine:
!-----------------------------------------------------------------------
subroutine jac_dummy_example(neq, t, y, ml, mu, pd, nrowpd, common_data)
    use odepack_common
    implicit none
    integer, intent(in) :: neq, ml, mu, nrowpd
    double precision, intent(in) :: t, y(neq)
    double precision, intent(out) :: pd(nrowpd, neq)
    type(odepack_common_data), intent(inout) :: common_data
    
    ! Silence unused variable warnings
    call int_fake_use(neq)
    call r8_fake_use(t)
    call r8_array_fake_use(y)
    call int_fake_use(ml)
    call int_fake_use(mu)
    call int_fake_use(nrowpd)
    
    ! Initialize outputs
    pd = 0.0d0
    common_data%ierr = 0

    contains
    subroutine r8_array_fake_use(x)
    implicit none
    double precision, intent(in) :: x(:)
    
    ! Use array in a way that will never execute
    if (size(x) < 0) then
        continue  ! Never executes since size() is always >= 0
    end if
end subroutine r8_array_fake_use
    
    ! Dummy routine - not called when JT=2
end subroutine jac_dummy_example

!-----------------------------------------------------------------------
! Alternative: Single "universal" fake use subroutine using transfer
!-----------------------------------------------------------------------
subroutine universal_fake_use(var)
    implicit none
    class(*), intent(in) :: var
    
    ! This uses unlimited polymorphic dummy argument (Fortran 2003+)
    ! Works with any type but requires modern Fortran
    continue  ! Do nothing, but var is "used"
end subroutine universal_fake_use