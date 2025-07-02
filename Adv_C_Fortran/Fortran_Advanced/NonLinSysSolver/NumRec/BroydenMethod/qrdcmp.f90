SUBROUTINE qrdcmp(a,n,np,c,d,sing)
    INTEGER n,np
    REAL a(np,np),c(n),d(n)
    LOGICAL sing
    ! Constructs the QR decomposition of a(1:n,1:n), with physical dimension np. The upper
    ! triangular matrix R is returned in the upper triangle of a, except for the diagonal elements
    ! of R which are returned in d(1:n). The orthogonal matrix Q is represented as a product of
    ! n − 1 Householder matrices Q1 ... Qn-1, where Qj = 1 - uj ⊗ uj /cj . The ith component
    ! of uj is zero for i = 1, ..., j - 1 while the nonzero components are returned in a(i,j) for
    ! i = j, ..., n. sing returns as true if singularity is encountered during the decomposition,
    ! but the decomposition is still completed in this case.
    INTEGER i,j,k
    REAL scale,sigma,sum,tau
    sing = .false.
    do k = 1, n-1
        scale = 0.
        do i = k, n
            scale = max(scale, abs(a(i,k)))
        enddo
        if (scale .eq. 0.) then !! Singular case.
            sing = .true.
            c(k) = 0.
            d(k) = 0.
            else !! Form Qk and Qk·A.
            do i = k, n
                a(i,k) = a(i,k)/scale
            enddo
            sum = 0.
            do i = k, n
                sum = sum + a(i,k)**2
            enddo
            sigma  = sign(sqrt(sum),a(k,k))
            a(k,k) = a(k,k)+sigma
            c(k)   = sigma*a(k,k)
            d(k)   = -scale*sigma
            do j = k+1, n
                sum = 0.
                do i = k, n
                    sum = sum + a(i,k)*a(i,j)
                enddo
                tau = sum/c(k)
                do i = k, n
                    a(i,j)=a(i,j)-tau*a(i,k)
                enddo
            enddo
        endif
    enddo
    d(n) = a(n,n)
    if(d(n) .eq. 0.) sing = .true.
    return
END SUBROUTINE