module romberg_utils
    implicit none
    
    ! Module parameters
    integer, parameter :: MAXITER = 10
    private :: MAXITER
    
    ! Public interface
    public :: romberg_method, trapezoidal_rule, print_table
    
    ! Abstract interface for user-defined functions
    abstract interface
        real(8) function func_interface(x)
            real(8), intent(in) :: x
        end function func_interface
    end interface
    
contains

    !---------------------------------------------------------------------------
    ! TRAPEZOIDAL_RULE - Compute definite integral using trapezoidal rule
    !---------------------------------------------------------------------------
    ! PURPOSE:
    !   Approximates the definite integral of f(x) from a to b using the
    !   trapezoidal rule with n equal intervals.
    !
    ! INPUTS:
    !   f     - Function to integrate (must match func_interface)
    !   a     - Lower limit of integration
    !   b     - Upper limit of integration
    !   n     - Number of intervals (must be positive)
    !
    ! OUTPUT:
    !   trapezoidal_rule - Approximation to the integral
    !
    ! FORMULA:
    !   ∫[a,b] f(x)dx ≈ h/2 * [f(a) + f(b)] + h * Σ[i=1 to n-1] f(a + i*h)
    !   where h = (b-a)/n
    !---------------------------------------------------------------------------
    real(8) function trapezoidal_rule(f, a, b, n)
        procedure(func_interface) :: f
        real(8), intent(in) :: a, b
        integer, intent(in) :: n
        
        real(8) :: h, sum_val, x
        integer :: i
        
        ! Calculate step size
        h = (b - a) / real(n, 8)
        
        ! Initialize sum with endpoint contributions
        sum_val = 0.5d0 * (f(a) + f(b))
        
        ! Add interior point contributions
        do i = 1, n - 1
            x = a + real(i, 8) * h
            sum_val = sum_val + f(x)
        end do
        
        ! Apply step size
        trapezoidal_rule = h * sum_val
    end function trapezoidal_rule
    
    !---------------------------------------------------------------------------
    ! ROMBERG_METHOD - Main Romberg integration routine
    !---------------------------------------------------------------------------
    ! PURPOSE:
    !   Computes the definite integral of f(x) from a to b using the Romberg
    !   method with Richardson extrapolation. The method uses successive
    !   halving of the step size and extrapolation to achieve high accuracy.
    !
    ! INPUTS:
    !   f           - Function to integrate (must match func_interface)
    !   a           - Lower limit of integration
    !   b           - Upper limit of integration
    !   tol         - Convergence tolerance (positive)
    !
    ! OUTPUTS:
    !   J           - Romberg table of approximations J(i,j)
    !   final_n     - Final iteration number where convergence occurred
    !   converged   - .true. if method converged, .false. otherwise
    !
    ! ALGORITHM:
    !   1. Compute J(i,1) using trapezoidal rule with 2^(i-1) intervals
    !   2. Apply Richardson extrapolation: J(i,j) = J(i,j-1) + 
    !      [J(i,j-1) - J(i-1,j-1)] / (4^(j-1) - 1)
    !   3. Estimate error: ε(i,i) = [J(i,i) - J(i-1,i-1)] / (4^(i-1) - 1)
    !   4. Check convergence: |ε(i,i)| ≤ tol
    !   5. Final result: J(i,i) with error correction
    !---------------------------------------------------------------------------
    subroutine romberg_method(f, a, b, tol, J, final_n, converged)
        procedure(func_interface) :: f
        real(8), intent(in) :: a, b, tol
        real(8), intent(out) :: J(MAXITER, MAXITER)
        integer, intent(out) :: final_n
        logical, intent(out) :: converged
        
        integer :: i, l, n_points
        real(8) :: epsilon_val, power_of_4
        
        ! Initialize outputs
        J = 0.0d0
        converged = .false.
        final_n = 0
        
        ! Validate inputs
        if (tol <= 0.0d0) then
            write(*,'(A)') 'Error: Tolerance must be positive'
            return
        endif
        
        if (a >= b) then
            write(*,'(A)') 'Error: Lower limit must be less than upper limit'
            return
        endif
        
        ! Step 1: Calculate J_11 with n = 1 (h = b-a)
        n_points = 1
        J(1, 1) = trapezoidal_rule(f, a, b, n_points)
        
        write(*,'(A)') 'Step 1:'
        write(*,'(A, E19.12)') 'J_11 = ', J(1, 1)
        
        ! Step 2: Calculate J_21 with n = 2 (h = (b-a)/2)
        n_points = 2
        J(2, 1) = trapezoidal_rule(f, a, b, n_points)
        
        write(*,*)
        write(*,'(A)') 'Step 2:'
        write(*,'(A, F16.14)') 'J_21 = ', J(2, 1)
        
        ! Calculate first error estimate
        epsilon_val = (J(2, 1) - J(1, 1)) / (4.0d0 - 1.0d0)
        write(*,'(A, E19.12)') 'ε_21 = (J_21 - J_11)/(4-1) = ', epsilon_val
        
        ! Calculate J_22
        J(2, 2) = J(2, 1) + epsilon_val
        write(*,'(A, F16.14)') 'J_22 = J_21 + ε_21 = ', J(2, 2)
        
        ! Check convergence for first iteration
        if (abs(epsilon_val) <= tol) then
            converged = .true.
            final_n = 2
            write(*,'(A, E19.12, A, E10.4)') 'Converged! |ε_21| = ', abs(epsilon_val), ' ≤ ', tol
            call print_table(J, 2)
            return
        else
            write(*,'(A, E19.12, A, E10.4)') 'Not converged: |ε_21| = ', abs(epsilon_val), ' > ', tol
        endif
        
        ! Continue with more iterations
        do i = 3, MAXITER
            n_points = 2**(i-1)
            
            ! Calculate J_i1 using trapezoidal rule
            J(i, 1) = trapezoidal_rule(f, a, b, n_points)
            
            write(*,*)
            write(*,'(A, i1, A)') 'Step ', i, ':'
            write(*,'(A, i1, A, F16.14)') 'J_', i, '1 = ', J(i, 1)
            
            ! Richardson extrapolation to fill the i-th row
            do l = 2, i
                power_of_4 = 4.0d0**(l-1)
                J(i, l) = J(i, l-1) + (J(i, l-1) - J(i-1, l-1)) / (power_of_4 - 1.0d0)
                write(*,'(A, i1, i1, A, F16.14)') 'J_', i, l, ' = ', J(i, l)
            end do
            
            ! Calculate error estimate for the diagonal element
            power_of_4 = 4.0d0**(i-1)
            epsilon_val = (J(i, i) - J(i-1, i-1)) / (power_of_4 - 1.0d0)
            write(*,'(A, i1, i1, A, i1, i1, A, i1, i1, A, i1, A, E19.12)') &
            'ε_', i, i, ' = (J_', i, i, ' - J_', i-1, i-1, ')/(', int(power_of_4), '-1) = ', epsilon_val
            
            ! Update J_ii with error correction
            J(i, i) = J(i, i-1) + epsilon_val
            write(*,'(A, i1, i1, A, i1, i1, A, i1, i1, A, F16.14)') &
            'J_', i, i, ' = J_', i, i-1, ' + ε_', i, i, ' = ', J(i, i)
            
            ! Check convergence
            if (abs(epsilon_val) <= tol) then
                converged = .true.
                final_n = i
                write(*,'(A, i1, i1, A, E19.12, A, E10.4)') &
                'Converged! |ε_', i, i, '| = ', abs(epsilon_val), ' ≤ ', tol
                exit
            else
                write(*,'(A, i1, i1, A, E19.12, A, E10.4)') &
                'Not converged: |ε_', i, i, '| = ', abs(epsilon_val), ' > ', tol
            endif
        end do
        
        ! Check if maximum iterations reached
        if (.not. converged) then
            write(*,'(A)') 'Warning: Did not converge within maximum iterations'
            final_n = MAXITER
        endif
        
        ! Print the final Romberg table
        call print_table(J, final_n)
        
    end subroutine romberg_method
    
    !---------------------------------------------------------------------------
    ! PRINT_TABLE - Display the Romberg table
    !---------------------------------------------------------------------------
    ! PURPOSE:
    !   Prints the Romberg table of approximations in a formatted manner.
    !   Shows the triangular structure of the J(i,j) values.
    !
    ! INPUTS:
    !   J - Romberg table of approximations
    !   n - Number of rows to print (must be ≤ MAXITER)
    !
    ! OUTPUT:
    !   Formatted table printed to standard output
    !---------------------------------------------------------------------------
    subroutine print_table(J, n)
        real(8), intent(in) :: J(MAXITER, MAXITER)
        integer, intent(in) :: n
        integer :: i, l
        
        ! Validate input
        if (n <= 0 .or. n > MAXITER) then
            write(*,'(A)') 'Error: Invalid table size for printing'
            return
        endif
        
        write(*,*)
        write(*,'(A)') 'Romberg Table (J values):'
        write(*,'(A)') '-------------------------'
        
        do i = 1, n
            write(*,'(A,I0,A)', advance='no') '  Row ', i, ': '
            do l = 1, i
                write(*,'(F12.8,2X)', advance='no') J(i, l)
            end do
            write(*,*)
        end do
        write(*,*)
        
    end subroutine print_table

end module romberg_utils