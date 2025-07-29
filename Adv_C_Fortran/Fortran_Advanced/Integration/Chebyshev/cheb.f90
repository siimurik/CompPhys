PROGRAM ChebyshevIntegral
    INTEGER, PARAMETER         :: n = 50
    DOUBLE PRECISION           :: a, b, c(n), cint(n)
    DOUBLE PRECISION, EXTERNAL :: func

    ! Set the interval [a, b] and initialize arrays c and cint
    a = -2.0D0
    b =  2.0D0
    c =  0.0D0
    cint = 0.0D0

    ! Call the chebft subroutine to compute Chebyshev coefficients for the function
    call chebft(a, b, c, n, func)
    
    ! Call the chint subroutine to compute Chebyshev coefficients for the integral
    call chint(a, b, c, cint, n)

    ! Print the result (integral) obtained from Chebyshev coefficients
    print*, "Result: ", cint(1)

END PROGRAM ChebyshevIntegral

! Function to be approximated using Chebyshev polynomials
FUNCTION func(x)
    DOUBLE PRECISION :: func, x
    func = EXP(-x**2)
END FUNCTION func

! Subroutine to compute Chebyshev coefficients for a given function
SUBROUTINE chebft(a, b, c, n, func)
    !Chebyshev fit: Given a function func, lower and upper limits of the interval [a,b], and
    !a maximum degree n, this routine computes the n coefficients ck such that func(x) ≈ Pn
    !k=1 ckTk−1(y)] − c1/2, where y and x are related by (5.8.10). This routine is to be
    !used with moderately large n (e.g., 30 or 50), the array of c’s subsequently to be truncated
    !at the smaller value m such that cm+1 and subsequent elements are negligible.
    !Parameters: Maximum expected value of n, and π.
    INTEGER             :: n 
    INTEGER, PARAMETER  :: NMAX = 50
    DOUBLE PRECISION    :: a, b, c(n)
    DOUBLE PRECISION, EXTERNAL  :: func
    DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793d0
    INTEGER             :: j, k
    DOUBLE PRECISION    :: bma, bpa, fac, y, f(NMAX)
    DOUBLE PRECISION    :: sum

    ! Check if n exceeds the maximum allowed value
    if (n > NMAX) then
        write(*, *) 'Warning! Array size larger than 50! Results may contain anomalous answers.'
    end if

    ! Calculate parameters for Chebyshev fit
    bma = 0.5D0*(b-a)
    bpa = 0.5D0*(b+a)

    ! Evaluate the function at Chebyshev points and store in array f
    do k = 1, n
        y = cos(PI*(k-0.5D0)/n)
        f(k) = func(y*bma+bpa)
    end do

    ! Compute Chebyshev coefficients using the coefficients formula
    fac = 2.0D0/n
    do j = 1, n
        sum = 0.0D0
        do k = 1, n
            sum = sum + f(k)*cos(PI*(j-1)*((k-0.5D0)/n))
        end do
        c(j) = fac*sum
    end do

    return
END SUBROUTINE chebft

! Subroutine to compute Chebyshev coefficients for the integral of a given function
SUBROUTINE chint(a, b, c, cint, n)
    INTEGER          :: n
    DOUBLE PRECISION :: a, b, c(n), cint(n)
    INTEGER          :: j
    DOUBLE PRECISION :: con, fac, sum

    ! Calculate normalization factor and initialize variables
    con = 0.25D0*(b-a)
    sum = 0.0D0
    fac = 1.0D0

    ! Compute Chebyshev coefficients for the integral
    do j = 2, n-1
        cint(j) = con*(c(j-1)-c(j+1))/(j-1)
        sum     = sum + fac*cint(j)
        fac     = -fac
    end do

    ! Special case for the last coefficient
    cint(n) = con*c(n-1)/(n-1)
    sum     = sum + fac*cint(n)

    ! Set the constant of integration
    cint(1) = 2.0D0*sum

    return
END SUBROUTINE chint
