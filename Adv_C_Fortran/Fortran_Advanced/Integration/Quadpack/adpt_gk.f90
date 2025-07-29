!===================================================================
! Compile and execute with:
!       gfortran -Wall -Wextra -o int adpt_gk.f90 quadpack.f90
!===================================================================

program integrate_example
    implicit none
    double precision :: a, b, epsabs, epsrel, result, abserr
    integer :: key, limit, neval, ier
    double precision, dimension(100) :: alist, blist, rlist, elist
    integer, dimension(100) :: iord
    integer :: last
    character(len=80) :: rule_description
    double precision :: analytical_result

    ! Print title screen and documentation
    write(*,'(A)') '================================================================================'
    write(*,'(A)') '              ADAPTIVE GAUSS-KRONROD QUADRATURE INTEGRATION'
    write(*,'(A)') '================================================================================'
    write(*,'(A)') 'Method: QUADPACK DQAGE Subroutine'
    write(*,'(A)') 'Algorithm: Adaptive integration using Gauss-Kronrod quadrature rules'
    write(*,*)
    write(*,'(A)') 'METHODOLOGY:'
    write(*,'(A)') '• Combines Gauss quadrature (optimal for smooth functions)'
    write(*,'(A)') '• With Kronrod extension (adds points for error estimation)'
    write(*,'(A)') '• Adaptive refinement subdivides intervals where error is large'
    write(*,'(A)') '• Continues until convergence or maximum subdivisions reached'
    write(*,*)

    ! Define the integration limits
    a = 0.0d0
    b = 1.0d0

    ! Set the desired absolute and relative accuracy
    epsabs = 1.0d-6
    epsrel = 1.0d-6

    ! Choose the integration rule (1 to 6)
    key = 2  ! Using 10-21 point Gauss-Kronrod pair

    ! Set the limit for the number of subintervals
    limit = 100

    ! Set rule description based on key value
    select case(key)
        case(1)
            rule_description = '7-15 Point Gauss-Kronrod Pair'
        case(2)
            rule_description = '10-21 Point Gauss-Kronrod Pair'
        case(3)
            rule_description = '15-31 Point Gauss-Kronrod Pair'
        case(4)
            rule_description = '20-41 Point Gauss-Kronrod Pair'
        case(5)
            rule_description = '25-51 Point Gauss-Kronrod Pair'
        case(6)
            rule_description = '30-61 Point Gauss-Kronrod Pair'
        case default
            rule_description = 'Unknown Rule'
    end select

    write(*,'(A)') 'INTEGRATION PARAMETERS:'
    write(*,'(A,A,A,I0,A)') 'Rule:     ', trim(rule_description), ' (KEY=', key, ')'
    write(*,'(A)') 'Function: f(x) = exp(-x²)'
    write(*,'(A,F4.1,A,F4.1,A)') 'Domain:   [', a, ',', b, ']'
    write(*,'(A,ES9.2)') 'Absolute tolerance: ', epsabs
    write(*,'(A,ES9.2)') 'Relative tolerance: ', epsrel
    write(*,'(A,I0)') 'Max subintervals:   ', limit
    write(*,*)
    write(*,'(A)') '================================================================================'
    write(*,'(A)') '                              INTEGRATION RESULTS'
    write(*,'(A)') '================================================================================'

    ! Call the dqage subroutine
    call dqage(f, a, b, epsabs, epsrel, key, limit, result, abserr, &
                neval, ier, alist, blist, rlist, elist, iord, last)

    ! Analytical result for comparison (integral of exp(-x²) from 0 to 1)
    analytical_result = 0.746824132812427d0

    ! Check for errors and print results
    if (ier == 0) then
        write(*,*)
        write(*,'(A)') 'SUCCESS: Integration completed successfully'
        write(*,*)
        write(*,'(A)') 'NUMERICAL RESULTS:'
        write(*,'(A,ES20.12)') 'Integral result:           ', result
        write(*,'(A,ES12.4)') 'Estimated absolute error:  ', abserr
        write(*,'(A,I0)') 'Number of function evals:  ', neval
        write(*,'(A,I0)') 'Number of subintervals:    ', last
        write(*,*)
        write(*,'(A)') 'ACCURACY ASSESSMENT:'
        write(*,'(A,ES20.12)') 'Analytical value:          ', analytical_result
        write(*,'(A,ES12.4)') 'Actual absolute error:     ', abs(result - analytical_result)
        write(*,'(A,ES12.4)') 'Error estimate accuracy:   ', abs(abserr - abs(result - analytical_result))
        write(*,*)
        if (abs(result - analytical_result) <= abserr) then
            write(*,'(A)') 'STATUS: Error estimate is reliable (actual error ≤ estimated error)'
        else
            write(*,'(A)') 'STATUS: Error estimate may be conservative'
        end if
    else
        write(*,*)
        write(*,'(A)') 'ERROR IN INTEGRATION:'
        write(*,'(A,I0)') 'Error code (ier):          ', ier
        select case(ier)
            case(1)
                write(*,'(A)') 'Error description:         Maximum number of subdivisions reached'
                write(*,'(A)') 'Suggestion:                Increase limit or relax tolerance'
            case(2)
                write(*,'(A)') 'Error description:         Roundoff error detected'
                write(*,'(A)') 'Suggestion:                Function may be discontinuous or tolerance too strict'
            case(3)
                write(*,'(A)') 'Error description:         Bad integrand behavior'
                write(*,'(A)') 'Suggestion:                Function may have singularities or discontinuities'
            case(6)
                write(*,'(A)') 'Error description:         Invalid input parameters'
                write(*,'(A)') 'Suggestion:                Check integration limits and tolerance values'
            case default
                write(*,'(A)') 'Error description:         Unknown error occurred'
        end select
        write(*,*)
        write(*,'(A)') 'PARTIAL RESULTS:'
        write(*,'(A,ES20.12)') 'Best approximation:        ', result
        write(*,'(A,ES12.4)') 'Error estimate:            ', abserr
        write(*,'(A,I0)') 'Function evaluations used: ', neval
        write(*,'(A,I0)') 'Subintervals processed:    ', last
    end if

    write(*,*)
    write(*,'(A)') '================================================================================'
    write(*,'(A)') 'QUADRATURE RULE INFORMATION:'
    write(*,'(A)') 'The Gauss-Kronrod method uses nested quadrature rules where:'
    write(*,'(A)') '• Gauss points provide the basic integral approximation'
    write(*,'(A)') '• Kronrod points (including Gauss points) give error estimation'
    write(*,'(A)') '• Higher key values = more points = better accuracy but more evaluations'
    write(*,'(A)') '================================================================================'

contains

    ! Define the integrand function
    double precision function f(x)
        double precision :: x
        f = exp(-x**2)  ! Example: f(x) = e^(-x^2)
    end function f

end program integrate_example
