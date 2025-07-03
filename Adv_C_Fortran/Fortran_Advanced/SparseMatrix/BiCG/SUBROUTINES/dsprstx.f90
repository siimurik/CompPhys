SUBROUTINE dsprstx(sa,ija,x,b,n)
    INTEGER n,ija(*)
    DOUBLE PRECISION b(n),sa(*),x(n)
    ! Multiply the transpose of a matrix in row-index sparse storage arrays sa and ija by a
    ! vector x(1:n), giving a vector b(1:n).
    INTEGER i,j,k
    if (ija(1).ne.n+2) stop 'mismatched vector and matrix in sprstx'
    do i = 1, n    !! Start with diagonal terms.
        b(i) = sa(i)*x(i)
    end do
    do i = 1, n    !! Loop over off-diagonal terms.
        do k = ija(i), ija(i+1)-1
            j = ija(k)
            b(j) = b(j) + sa(k)*x(i)
        end do
    end do
    return
END SUBROUTINE dsprstx