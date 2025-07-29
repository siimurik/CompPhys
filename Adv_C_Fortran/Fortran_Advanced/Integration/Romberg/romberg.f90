!===============================================================================
! ROMBERG INTEGRATION METHOD - DOCUMENTATION
!===============================================================================
!
! THEORY:
! -------
! The Romberg integration method is an extrapolation technique that uses the
! trapezoidal rule with successive halving of the step size to achieve high
! accuracy. It is based on Richardson extrapolation, which eliminates 
! successive terms in the error expansion of the trapezoidal rule.
!
! The trapezoidal rule has an error expansion of the form:
! I - T(h) = c₁h² + c₂h⁴ + c₃h⁶ + ...
! where I is the exact integral and T(h) is the trapezoidal approximation.
!
! By computing T(h) and T(h/2) and using Richardson extrapolation, we can
! eliminate the h² term to get a more accurate approximation.
!
! ALGORITHM STRUCTURE:
! -------------------
! The method creates a triangular table of approximations J(i,j) where:
! - J(i,1) = trapezoidal rule with 2^(i-1) intervals
! - J(i,j) = Richardson extrapolation of J(i,j-1) and J(i-1,j-1)
!
! The table looks like:
!     J₁₁
!     J₂₁  J₂₂
!     J₃₁  J₃₂  J₃₃
!     J₄₁  J₄₂  J₄₃  J₄₄
!     ...
!
! RICHARDSON EXTRAPOLATION FORMULA:
! ---------------------------------
! J(i,j) = J(i,j-1) + [J(i,j-1) - J(i-1,j-1)] / (4^(j-1) - 1)
!
! This formula eliminates the next term in the error expansion, giving:
! - J(i,2): eliminates h² term (denominator = 4¹ - 1 = 3)
! - J(i,3): eliminates h⁴ term (denominator = 4² - 1 = 15)
! - J(i,4): eliminates h⁶ term (denominator = 4³ - 1 = 63)
!
! ERROR ESTIMATION:
! ----------------
! The error estimate for J(i,i) is given by:
! ε(i,i) = [J(i,i) - J(i-1,i-1)] / (4^(i-1) - 1)
!
! This estimates the error in the diagonal element and is used for
! convergence testing.
!
! STEP-BY-STEP PROCEDURE:
! -----------------------
! 1. Calculate J₁₁ using trapezoidal rule with h = (b-a)
! 2. Calculate J₂₁ using trapezoidal rule with h = (b-a)/2
! 3. Calculate error estimate: ε₂₁ = (J₂₁ - J₁₁)/(4-1)
! 4. Check if |ε₂₁| ≤ TOL. If yes, stop with result J₂₂ = J₂₁ + ε₂₁
! 5. If not converged, continue with J₃₁ using h = (b-a)/4
! 6. Calculate J₃₂ using Richardson extrapolation
! 7. Calculate error estimate: ε₃₂ = (J₃₂ - J₂₂)/(16-1)
! 8. Check convergence. If yes, stop with J₃₃ = J₃₂ + ε₃₂
! 9. Continue until convergence or maximum iterations reached
!
! CONVERGENCE CRITERIA:
! --------------------
! The method stops when |ε(i,i)| ≤ TOL, where TOL is the specified tolerance.
! The final result is J(i,i) after applying the error correction.
!
! ADVANTAGES:
! -----------
! - Rapid convergence for smooth functions
! - Automatic error estimation
! - Efficient use of previously computed values
! - Higher order accuracy than simple trapezoidal rule
!
! LIMITATIONS:
! -----------
! - Requires smooth functions (continuous derivatives)
! - May not work well for functions with singularities
! - Convergence may be slow for oscillatory functions
!
! IMPLEMENTATION NOTES:
! --------------------
! - The trapezoidal rule is implemented for exactly n intervals
! - J(i,1) uses 2^(i-1) intervals, so J₁₁ uses 1 interval, J₂₁ uses 2, etc.
! - The table is filled row by row, with Richardson extrapolation applied
! - Error estimates are calculated for diagonal elements only
! - The method terminates when error estimate is within tolerance
!
! MATHEMATICAL BACKGROUND:
! -----------------------
! The denominators in Richardson extrapolation (3, 15, 63, ...) come from
! the formula 4^k - 1 where k is the extrapolation level:
! - Level 1: 4¹ - 1 = 3   (eliminates h² term)
! - Level 2: 4² - 1 = 15  (eliminates h⁴ term)
! - Level 3: 4³ - 1 = 63  (eliminates h⁶ term)
!
! This pattern continues, with each level eliminating the next even power
! of h in the error expansion.
!
! EXAMPLE APPLICATION:
! -------------------
! Problem 1: ∫₀² e^(-x) dx with TOL = 10⁻³
! - Exact value: 1 - e⁻² ≈ 0.864664716
! - Expected to converge quickly due to smooth function
!
! Problem 2: ∫₀² (π/4)x⁴cos(πx/4) dx with TOL = 10⁻⁴
! - More challenging due to oscillatory behavior
! - May require more iterations to achieve convergence
!
! Refernece: Advanced Engineering Mathematics 10th Edition by Erwin Kreyszig
!            PAGE 840, PROBLEM SET 19.5, EX. 26
!===============================================================================

program romberg_integration
    implicit none
    
    ! Parameters
    integer, parameter :: MAXITER = 10
    real(8), parameter :: TOL1 = 1.0d-3    ! Tolerance for first problem
    real(8), parameter :: TOL2 = 1.0d-4    ! Tolerance for second problem
    real(8), parameter :: PI = 4.D0 * atan(1.D0)
    
    ! Variables
    real(8) :: J(MAXITER, MAXITER)  ! Romberg table
    real(8) :: a, b, tol
    integer :: n
    logical :: converged
    
    ! Problem 1: Integral of e^(-x) from 0 to 2
    write(*,*) '=============================================='
    write(*,*) 'Problem 1: Integral of e^(-x) from 0 to 2'
    write(*,*) 'Tolerance = ', TOL1
    write(*,*) '=============================================='
    write(*,*)
    
    a = 0.0d0
    b = 2.0d0
    tol = TOL1
    
    call romberg_method(f1, a, b, tol, J, n, converged)
    
    write(*,*) 'Exact value: ', 1.0d0 - exp(-2.0d0)
    write(*,*) 'Converged: ', converged
    if (converged) then
        write(*,*) 'Final result: ', J(n, n)
    endif
    
    write(*,*)
    write(*,*) '=============================================='
    write(*,*) 'Problem 2: Integral of (π/4)x^4 cos(πx/4) from 0 to 2'
    write(*,*) 'Tolerance = ', TOL2
    write(*,*) '=============================================='
    write(*,*)
    
    a = 0.0d0
    b = 2.0d0
    tol = TOL2
    
    call romberg_method(f2, a, b, tol, J, n, converged)
    
    write(*,*) 'Converged: ', converged
    if (converged) then
        write(*,*) 'Final result: ', J(n, n)
    endif

    !---------------- end of main code ----------------

contains

    ! Function 1: f(x) = e^(-x)
    real(8) function f1(x)
        real(8), intent(in) :: x
        f1 = exp(-x)
    end function f1
    
    ! Function 2: f(x) = (π/4)x^4 cos(πx/4)
    real(8) function f2(x)
        real(8), intent(in) :: x
        f2 = (PI/4.0d0) * x**4 * cos(PI * x / 4.0d0)
    end function f2
    
    ! Trapezoidal rule
    real(8) function trapezoidal_rule(f, a, b, n)
        interface
            real(8) function f(x)
                real(8), intent(in) :: x
            end function f
        end interface
        real(8), intent(in) :: a, b
        integer, intent(in) :: n
        
        real(8) :: h, sum_val, x
        integer :: i
        
        h = (b - a) / real(n, 8)
        sum_val = 0.5d0 * (f(a) + f(b))
        
        do i = 1, n - 1
            x = a + real(i, 8) * h
            sum_val = sum_val + f(x)
        end do
        
        trapezoidal_rule = h * sum_val
    end function trapezoidal_rule
    
    ! Main Romberg integration routine
    subroutine romberg_method(f, a, b, tol, J, final_n, converged)
        interface
            real(8) function f(x)
                real(8), intent(in) :: x
            end function f
        end interface
        real(8), intent(in) :: a, b, tol
        real(8), intent(out) :: J(MAXITER, MAXITER)
        integer, intent(out) :: final_n
        logical, intent(out) :: converged
        
        integer :: i, l, n_points
        real(8) :: epsilon_val, power_of_4
        
        ! Initialize
        J = 0.0d0
        converged = .false.
        
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
        
        ! Check convergence
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
            write(*,'(A, i1, i1, A, i1, i1, A, i1, i1, A, F16.14)') 'J_', i, i, ' = J_', i, i-1, ' + ε_', i, i, ' = ', J(i, i)
            
            ! Check convergence
            if (abs(epsilon_val) <= tol) then
                converged = .true.
                final_n = i
                write(*,'(A, i1, i1, A, E19.12, A, E10.4)') 'Converged! |ε_', i, i, '| = ', abs(epsilon_val), ' ≤ ', tol
                exit
            else
                write(*,'(A, i1, i1, A, E19.12, A, E10.4)') 'Not converged: |ε_', i, i, '| = ', abs(epsilon_val), ' > ', tol
            endif
        end do
        
        if (.not. converged) then
            write(*,'(a)') 'Did not converge within maximum iterations'
            final_n = MAXITER
        endif
        
        ! Print the Romberg table
        call print_table(J, final_n)
        
    end subroutine romberg_method
    
    ! Print the Romberg table
    subroutine print_table(J, n)
        real(8), intent(in) :: J(MAXITER, MAXITER)
        integer, intent(in) :: n
        integer :: i, l
        
        write(*,*)
        write(*,*) 'Romberg Table (J values):'
        write(*,*) '-------------------------'
        do i = 1, n
            write(*,'(A,I0,A)', advance='no') '  Row ', i, ': '
            do l = 1, i
                write(*,'(F12.8,2X)', advance='no') J(i, l)
            end do
            write(*,*)
        end do
        write(*,*)
        
    end subroutine print_table
    
end program romberg_integration