!=========================================================
! Compile and rund with:
!   gfortran nls.f90 minpack.f90 -o solver && ./solver 
!=========================================================

program nonlinear_solver_example
    use minpack_module, only: hybrd
    implicit none
    external :: fcn
    
    ! Problem parameters
    integer, parameter :: n = 2                  ! Number of equations/variables
    real(8), parameter :: tol = 1.0d-4           ! Solution tolerance
    real(8), parameter :: initial_guess(n) = &   ! Starting point
                          [10.0d0, 10.0d0]       ! Bad guess example: [15.0d0, -2.0d0]
    
    ! Solver variables
    integer :: info, func_evals
    real(8) :: solution(n), residuals(n)
    real(8) :: jacobian(n,n), workspace_r(n*(n+1)/2)
    real(8) :: diag(n), qtf(n), workspace_1(n), workspace_2(n), &
               workspace_3(n), workspace_4(n)
    
    ! Solve the system
    call solve_system(solution, residuals, info, func_evals)
    
    ! Display results
    call print_results(solution, residuals, info, func_evals)
    
contains

    subroutine solve_system(x, fvec, status, eval_count)
        real(8), intent(out) :: x(n), fvec(n)
        integer, intent(out) :: status, eval_count
        
        x = initial_guess
        
        call hybrd(fcn, n, x, fvec, tol, 1000, & 
                  n-1, n-1, 0.0d0, diag, 1, 100.0d0, 1, &
                  status, eval_count, jacobian, n, workspace_r, &
                  size(workspace_r), qtf, workspace_1, workspace_2, &
                  workspace_3, workspace_4)
    end subroutine solve_system


    subroutine print_results(x, fvec, status, eval_count)
        real(8), intent(in) :: x(n), fvec(n)
        integer, intent(in) :: status, eval_count
        character(len=30) :: status_message
        
        select case(status)
            case(0); status_message = "Invalid input"
            case(1); status_message = "Solution converged"
            case(2); status_message = "Max evaluations reached"
            case(3); status_message = "Tolerance too small"
            case(4); status_message = "Poor progress (Jacobian)"
            case(5); status_message = "Poor progress (iterations)"
            case default; status_message = "Unknown status"
        end select
        
        print *, "----------------------------------------"
        print *, " NONLINEAR SYSTEM SOLVER RESULTS"
        print *, "----------------------------------------"
        print '(A,2F12.6)', " Solution:       ", x
        print '(A,2ES12.4)', " Residuals:      ", fvec
        print '(A,I0,A)',    " Function evals: ", eval_count, " calls"
        print '(A,I0,A,A)',  " Solver status:  ", status, " - ", trim(status_message)
        print *, "----------------------------------------"
        
        if (maxval(abs(fvec)) < tol) then
            print *, " Verification: Solution satisfies tolerance!"
        else
            print *, " Warning: Residuals above tolerance level"
        end if
    end subroutine print_results

end program nonlinear_solver_example

subroutine fcn(n, x, fvec)
    integer, intent(in) :: n
    real(8), intent(in) :: x(n)
    real(8), intent(out) :: fvec(n)
    
    fvec(1) = -13.0d0 + x(1) + ((-x(2) + 5.0d0)*x(2) - 2.0d0)*x(2)
    fvec(2) = -29.0d0 + x(1) + ((x(2) + 1.0d0)*x(2) - 14.0d0)*x(2)
end subroutine fcn