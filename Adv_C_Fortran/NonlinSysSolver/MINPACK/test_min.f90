!===========================================================
! Compile and rund with:
!   gfortran test_min.f90 minpack.f90 -o solver && ./solver 
!===========================================================

program minpack_example
    use minpack_module
    implicit none
    integer, parameter :: n = 2
    integer :: info, nfev
    real(8) :: x(n), fvec(n), xtol
    real(8) :: diag(n), fjac(n,n), r(n*(n+1)/2), qtf(n), wa1(n), wa2(n), wa3(n), wa4(n)
    external :: fcn
    
    ! Initial guess
    x = [15.0d0, 10.0d0]
    
    ! Tolerance
    xtol = 1.0d-6
    
    ! Call HYBRD solver
    call hybrd(fcn, n, x, fvec, xtol, 1000, n-1, n-1, 0.0d0, diag, &
               1, 100.0d0, 1, info, nfev, fjac, n, r, n*(n+1)/2, qtf, &
               wa1, wa2, wa3, wa4)
    
    ! Output results
    print *, 'Solution:'
    print *, 'x(1) =', x(1)
    print *, 'x(2) =', x(2)
    print *, 'Function values at solution:'
    print *, 'f(1) =', fvec(1)
    print *, 'f(2) =', fvec(2)
    print *, 'Number of function evaluations:', nfev
    print *, 'Termination flag (info):', info
    
end program minpack_example

! User-defined function (same as your CALFUN)
subroutine fcn(n, x, fvec)
    integer, intent(in) :: n
    real(8), intent(in) :: x(n)
    real(8), intent(out) :: fvec(n)
    
    fvec(1) = -13.0d0 + x(1) + ((-x(2) + 5.0d0)*x(2) - 2.0d0)*x(2)
    fvec(2) = -29.0d0 + x(1) + ((x(2) + 1.0d0)*x(2) - 14.0d0)*x(2)
end subroutine fcn