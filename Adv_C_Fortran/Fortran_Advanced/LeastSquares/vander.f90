SUBROUTINE vander(x, w, q, n)
    IMPLICIT NONE
    INTEGER, INTENT(IN)           :: n
    DOUBLE PRECISION, INTENT(IN)  :: x(n),q(n)
    DOUBLE PRECISION, INTENT(OUT) :: w(n)
    ! Solves the Vandermonde linear system ∑N i=1 xk−1
    ! i wi = qk (k = 1, ..., N). Input consists
    ! of the vectors x(1:n) and q(1:n); the vector w(1:n) is output.
    ! Parameters: NMAX is the maximum expected value of n.
    INTEGER i,j,k, NMAX
    PARAMETER (NMAX=100)
    DOUBLE PRECISION b,s,t,xx,c(NMAX)
    if (n .eq. 1) then
        w(1) = q(1)
    else
        do i = 1, n !! Initialize array.
            c(i) = 0.d0
        end do
        c(n) = -x(1)    !! Coefficients of the master polynomial are found by recur-
        do i = 2, n     !! sion.
            xx = -x(i)
            do j = n+1-i, n-1
                c(j) = c(j) + xx*c(j+1)
            end do
            c(n) = c(n) + xx
        end do
        do i = 1, n !! Each subfactor in turn
            xx = x(i)
            t = 1.d0
            b = 1.d0
            s = q(n)
            do k = n, 2, -1 !! is synthetically divided,
                b = c(k) + xx*b
                s = s + q(k-1)*b !! matrix-multiplied by the right-hand side,
                t = xx*t + b
            end do
            w(i) = s/t !! and supplied with a denominator.
        end do
    end if
    return
END SUBROUTINE 