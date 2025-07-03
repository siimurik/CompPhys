SUBROUTINE asolve(n, b, x, itrnsp)
    INTEGER :: n,itrnsp,ija,NMAX,i
    DOUBLE PRECISION :: x(n), b(n), sa
    PARAMETER (NMAX=1000)
    COMMON /mat/ sa(NMAX),ija(NMAX) !! The matrix is stored somewhere.
    do i = 1, n
        x(i) = b(i)/sa(i)
    end do
    return
END SUBROUTINE asolve