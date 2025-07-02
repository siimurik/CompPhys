SUBROUTINE rsolv(a,n,np,d,b)
    INTEGER n, np
    REAL a(np,np), b(n), d(n)
    ! Solves the set of n linear equations R·x = b, where R is an upper triangular matrix stored
    ! in a and d. a and d are input as the output of the routine qrdcmp and are not modified.
    ! b(1:n) is input as the right-hand side vector, and is overwritten with the solution vector
    ! on output.
    INTEGER i,j
    REAL sum
    b(n) = b(n)/d(n)
    do i = n-1, 1, -1
        sum = 0.
        do j = i+1, n
            sum = sum + a(i,j)*b(j)
        enddo
        b(i) = (b(i)-sum)/d(i)
    enddo
    return
END SUBROUTINE