!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Van der Pol Oscillator with ODEPACK LSODA
!------------------------------------------------------------------------------
! Compile and execute:
!    gfortran -O3 -std=legacy -o vdp van_der_pol_lsoda.f90 opkdmain.f
! 
! This program solves the stiff Van der Pol oscillator:
!    x'' - μ(1-x²)x' + x = 0
! Converted to system: y1' = y2, y2' = μ(1-y1²)y2 - y1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Global module for Van der Pol parameter
module vdp_params
    implicit none
    double precision :: mu_global
end module vdp_params

program van_der_pol_solver
    use vdp_params
    implicit none
    
    ! Problem parameters
    integer, parameter :: neq = 2           ! Number of equations
    integer, parameter :: lrw = 22 + 9*neq + neq**2  ! Real work array size
    integer, parameter :: liw = 20 + neq    ! Integer work array size
    integer, parameter :: nout = 1000       ! Number of output points
    
    ! Variables for LSODA
    double precision :: y(neq)              ! Solution vector
    double precision :: t, tout             ! Time variables
    double precision :: rtol, atol          ! Tolerances
    integer :: itol, itask, istate, iopt    ! LSODA control parameters
    double precision :: rwork(lrw)          ! Real work array
    integer :: iwork(liw)                   ! Integer work array
    integer :: jt                           ! Jacobian type
    
    ! Problem-specific variables
    double precision :: t_final             ! Final integration time
    double precision :: dt_out              ! Output time step
    integer :: i                            ! Loop counter
    integer :: nsteps, nfe, nje             ! Statistics variables
    
    ! Output file
    integer, parameter :: output_unit = 10
    
    ! External procedure declarations
    external :: vdp_rhs, jac_dummy
    
    ! Print documentation header
    call print_header()
    
    ! Problem setup
    mu_global = 50.0d0          ! Large μ makes the problem stiff
    t = 0.0d0                   ! Initial time
    t_final = 20.0d0            ! Final time
    dt_out = t_final / dble(nout) ! Output time step
    
    ! Initial conditions
    y(1) = 2.0d0                ! x(0) = 2.0 (displacement)
    y(2) = 0.0d0                ! x'(0) = 0.0 (velocity)
    
    ! LSODA parameters
    rtol = 1.0d-8               ! Relative tolerance
    atol = 1.0d-10              ! Absolute tolerance
    itol = 1                    ! Scalar tolerances
    itask = 1                   ! Normal computation mode
    istate = 1                  ! First call to LSODA
    iopt = 0                    ! No optional inputs
    jt = 2                      ! Internal numerical Jacobian
    
    ! Initialize work arrays
    do i = 1, lrw
        rwork(i) = 0.0d0
    end do
    do i = 1, liw
        iwork(i) = 0
    end do
    
    write(*,'(A)') 'PROBLEM SETUP:'
    write(*,'(A,F6.1)') 'Van der Pol parameter μ = ', mu_global
    write(*,'(A,F4.1,A,F4.1,A)') 'Time domain: [', t, ', ', t_final, ']'
    write(*,'(A,F4.1,A,F4.1,A)') 'Initial conditions: y1(0) = ', y(1), ', y2(0) = ', y(2)
    write(*,'(A,ES9.2,A,ES9.2)') 'Tolerances: rtol = ', rtol, ', atol = ', atol
    write(*,'(A,I0)') 'Output points: ', nout
    write(*,*)
    write(*,'(A)') 'INTEGRATION PROGRESS:'
    write(*,'(A)') '================================================================================'
    
    ! Open output file
    open(unit=output_unit, file='van_der_pol_solution.dat', status='replace')
    write(output_unit, '(A)') '# Van der Pol Oscillator Solution'
    write(output_unit, '(A,F6.1)') '# μ = ', mu_global
    write(output_unit, '(A)') '# t           y1(t)         y2(t)'
    
    ! Write initial condition
    write(output_unit, '(3ES15.6)') t, y(1), y(2)
    
    ! Integration loop
    do i = 1, nout
        tout = i * dt_out
        
        ! Call LSODA
        call dlsoda(vdp_rhs, neq, y, t, tout, itol, rtol, atol, itask, &
                   istate, iopt, rwork, lrw, iwork, liw, jac_dummy, jt)
        
        ! Check for errors
        if (istate < 0) then
            write(*,'(A,I0)') 'ERROR: LSODA failed with istate = ', istate
            call print_error_message(istate)
            exit
        end if
        
        ! Write solution to file
        write(output_unit, '(3ES15.6)') t, y(1), y(2)
        
        ! Progress indicator
        if (mod(i, nout/20) == 0) then
            write(*,'(A,F6.2,A,F8.4,A,ES12.4,A,ES12.4)') &
                'Progress: ', 100.0d0*i/nout, '% | t = ', t, &
                ' | y1 = ', y(1), ' | y2 = ', y(2)
        end if
    end do
    
    close(output_unit)
    
    ! Print final results and statistics
    write(*,*)
    write(*,'(A)') '================================================================================'
    write(*,'(A)') 'INTEGRATION COMPLETED SUCCESSFULLY'
    write(*,'(A)') '================================================================================'
    write(*,*)
    write(*,'(A)') 'FINAL SOLUTION:'
    write(*,'(A,F8.4)') 'Final time:        ', t
    write(*,'(A,ES15.6)') 'y1(t_final):       ', y(1)
    write(*,'(A,ES15.6)') 'y2(t_final):       ', y(2)
    write(*,*)
    
    ! Extract statistics from IWORK array
    nsteps = iwork(11)          ! Number of steps taken
    nfe = iwork(12)             ! Number of function evaluations
    nje = iwork(13)             ! Number of Jacobian evaluations
    
    write(*,'(A)') 'PERFORMANCE STATISTICS:'
    write(*,'(A,I0)') 'Number of time steps:      ', nsteps
    write(*,'(A,I0)') 'Function evaluations:      ', nfe
    write(*,'(A,I0)') 'Jacobian evaluations:      ', nje
    if (nsteps > 0) then
        write(*,'(A,F8.4)') 'Average step size:         ', t_final / dble(nsteps)
        write(*,'(A,F6.2)') 'Function evals per step:   ', dble(nfe) / dble(nsteps)
    end if
    write(*,*)
    write(*,'(A)') 'OUTPUT:'
    write(*,'(A)') 'Solution saved to: van_der_pol_solution.dat'
    write(*,'(A)') 'Use gnuplot: plot "van_der_pol_solution.dat" u 1:2 w l title "y1", "" u 1:3 w l title "y2"'
    write(*,'(A)') 'Phase plot: plot "van_der_pol_solution.dat" u 2:3 w l title "Phase Portrait"'
    write(*,*)
    write(*,'(A)') '================================================================================'

contains

    ! Print header information
    subroutine print_header()
        write(*,'(A)') '================================================================================'
        write(*,'(A)') '           VAN DER POL OSCILLATOR - ODEPACK LSODA SOLVER'
        write(*,'(A)') '================================================================================'
        write(*,'(A)') 'Differential Equation: x" - μ(1-x²)x" + x = 0'
        write(*,'(A)') 'System Form: y1" = y2, y2" = μ(1-y1²)y2 - y1'
        write(*,*)
        write(*,'(A)') 'SOLVER CHARACTERISTICS:'
        write(*,'(A)') '• LSODA: Automatic stiff/non-stiff method switching'
        write(*,'(A)') '• Adams methods for non-stiff regions'
        write(*,'(A)') '• BDF methods for stiff regions'
        write(*,'(A)') '• Adaptive step size and order control'
        write(*,'(A)') '• Self-monitoring integration process'
        write(*,*)
        write(*,'(A)') 'PROBLEM CHARACTERISTICS:'
        write(*,'(A)') '• Large μ creates stiff relaxation oscillations'
        write(*,'(A)') '• Exhibits both fast (stiff) and slow (non-stiff) dynamics'
        write(*,'(A)') '• Excellent test case for automatic stiffness detection'
        write(*,'(A)') '================================================================================'
        write(*,*)
    end subroutine print_header
    
    ! Print error message based on istate
    subroutine print_error_message(istate)
        integer, intent(in) :: istate
        write(*,*)
        write(*,'(A)') 'ERROR DIAGNOSTICS:'
        select case(istate)
            case(-1)
                write(*,'(A)') 'Excess work done (check tolerances or time span)'
            case(-2)
                write(*,'(A)') 'Excess accuracy requested (relax tolerances)'
            case(-3)
                write(*,'(A)') 'Illegal input detected'
            case(-4)
                write(*,'(A)') 'Repeated error test failures (check problem formulation)'
            case(-5)
                write(*,'(A)') 'Repeated convergence failures (check Jacobian)'
            case(-6)
                write(*,'(A)') 'Error weight became zero (check tolerances)'
            case default
                write(*,'(A)') 'Unknown error occurred'
        end select
    end subroutine print_error_message

end program van_der_pol_solver

! Van der Pol right-hand side function (external)
subroutine vdp_rhs(neq, t, y, ydot, rpar, ipar)
    use vdp_params
    implicit none
    integer :: neq, ipar(*)
    double precision :: t, y(neq), ydot(neq), rpar(*)
    
    ! Van der Pol equations: y1' = y2, y2' = μ(1-y1²)y2 - y1
    ydot(1) = y(2)
    ydot(2) = mu_global * (1.0d0 - y(1)**2) * y(2) - y(1)
end subroutine vdp_rhs

! Dummy Jacobian routine (LSODA will compute numerically)
subroutine jac_dummy(neq, t, y, ml, mu, pd, nrowpd, rpar, ipar)
    implicit none
    integer :: neq, ml, mu, nrowpd, ipar(*)
    double precision :: t, y(neq), pd(nrowpd, neq), rpar(*)
    ! This routine is not called when jt = 2
end subroutine jac_dummy