!=========================================================
! Compile and rund with:
!   gfortran fdjac.f90 fmin.f90 lnsrch.f90 lubksb.f90 ludcmp.f90 main_newt.f90 newt.f90 -o newt
!=========================================================

PROGRAM test_newton
    IMPLICIT NONE
    INTEGER, PARAMETER :: n = 2
    INTEGER :: eval_count
    REAL    :: x(n)
    LOGICAL :: check

    ! Initial guess
    x = [20.0, 0.05]

    ! Globally convergent Newton's method
    CALL newt(x, n, check, eval_count)

    ! Display results
    CALL print_results(x, n, check, eval_count)
END PROGRAM

SUBROUTINE funcv(n, x, fvec)
    INTEGER n
    REAL x(n), fvec(n)

    ! Example: Solve the system
    ! f1(x) = x1^2 + x2^2 - 4
    ! f2(x) = x1*x2 - 1
    fvec(1) = x(1)**2 + x(2)**2 - 4.0
    fvec(2) = x(1)*x(2) - 1.0
END SUBROUTINE

SUBROUTINE print_results(x, n, check, eval_count)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, eval_count
    REAL,    INTENT(IN) :: x(n) 
    LOGICAL, INTENT(IN) :: check

    IF (check) THEN
        PRINT *, "Converged to a local minimum, try a different initial guess."
    ELSE
        PRINT *, "----------------------------------------"
        PRINT *, " GLOBALLY CONVERGENT NEWTON'S METHOD"
        PRINT *, "----------------------------------------"
        PRINT '(A, *(F12.6))', " Solution:       ",  x
        print '(A,I5,A)',    " Function evals: ", eval_count, " calls"
        PRINT *, "----------------------------------------"
    END IF


END SUBROUTINE print_results

