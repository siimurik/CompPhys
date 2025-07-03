SUBROUTINE sprsax(sa, ija, x, b, n)
    INTEGER n,ija(*)
    REAL b(n),sa(*),x(n)
    ! Multiply a matrix in row-index sparse storage arrays sa and ija by a vector x(1:n), giving
    ! a vector b(1:n).
    INTEGER i,k
    if (ija(1).ne.n+2) stop 'mismatched vector and matrix in sprsax'
    do i = 1, n
        b(i) = sa(i)*x(i) !! Start with diagonal term.
        do k = ija(i), ija(i+1)-1 !! Loop over off-diagonal terms.
            b(i) = b(i) + sa(k)*x(ija(k))
        end do
    end do
    return
END SUBROUTINE sprsax