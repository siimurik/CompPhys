PROGRAM main
    IMPLICIT NONE
    INTEGER, PARAMETER :: n = 9, M = 4  ! Degree of polynomial + 1
    DOUBLE PRECISION :: x(n), y(n), k(n), V(n, M), x_sol(M)

    ! Input data
    x = [1.2D0, 1.3D0, 1.4D0, 1.5D0, 1.6D0, 1.7D0, 1.8D0, 1.9D0, 2.0D0]
    y = [9.08D0, 10.43D0, 11.9D0, 13.48D0, 15.19D0, 17.03D0, 19.01D0, 21.13D0, 23.39D0]
    k = [1.0D0, 1.0D0, 2.0D0, 5.0D0, 1.0D0, 4.0D0, 2.0D0, 2.0D0, 1.0D0]

    ! Construct weighted Vandermonde matrix
    CALL construct_vander(x, k, n, M, V)

    ! Solve using QR decomposition
    CALL solve_qr(V, y, k, n, M, x_sol)

    ! Print solution (flipped to match MATLAB's convention)
    PRINT *, "Solution (flipped) X ="
    PRINT *, x_sol
END PROGRAM main

SUBROUTINE solve_qr(V, y, k, n, M, x_sol)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, M
    DOUBLE PRECISION, INTENT(IN) :: V(n, M), y(n), k(n)
    DOUBLE PRECISION, INTENT(OUT) :: x_sol(M)
    DOUBLE PRECISION :: a(n, M), c(M), d(M), qt(n, n), u(n), vv(n), work(n), r(M, M)
    LOGICAL :: sing
    INTEGER :: i, j

    ! Copy V into a (qrdcmp will modify it)
    a = V

    ! Perform QR decomposition (modified for rectangular matrices)
    CALL qrdcmp_rect(a, n, M, c, d, sing)
    IF (sing) THEN
        PRINT *, "Warning: Matrix is singular or nearly singular."
    END IF

    ! Extract R (upper triangular part of a)
    r = 0.0D0
    DO j = 1, M
        DO i = 1, MIN(j, n)
            r(i, j) = a(i, j)
        END DO
    END DO

    ! Construct Q^T * (k .* y)
    work = k * y
    DO i = 1, M
        x_sol(i) = SUM(a(:, i) * work(:))
    END DO

    ! Solve R * x = Q^T * (k .* y) using back substitution
    DO i = M, 1, -1
        x_sol(i) = x_sol(i) / r(i, i)
        DO j = 1, i-1
            x_sol(j) = x_sol(j) - r(j, i) * x_sol(i)
        END DO
    END DO

    ! Flip the solution to match MATLAB's convention
    x_sol = x_sol(M:1:-1)
END SUBROUTINE solve_qr

SUBROUTINE qrdcmp_rect(a, m, n, c, d, sing)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m, n
    DOUBLE PRECISION, INTENT(INOUT) :: a(m, n)
    DOUBLE PRECISION, INTENT(OUT) :: c(n), d(n)
    LOGICAL, INTENT(OUT) :: sing
    INTEGER :: i, j, k
    DOUBLE PRECISION :: scale, sigma, sum, tau

    sing = .FALSE.
    DO k = 1, n
        scale = 0.0D0
        DO i = k, m
            scale = MAX(scale, ABS(a(i, k)))
        END DO
        IF (scale == 0.0D0) THEN
            sing = .TRUE.
            c(k) = 0.0D0
            d(k) = 0.0D0
        ELSE
            DO i = k, m
                a(i, k) = a(i, k) / scale
            END DO
            sum = 0.0D0
            DO i = k, m
                sum = sum + a(i, k)**2
            END DO
            sigma = SIGN(SQRT(sum), a(k, k))
            a(k, k) = a(k, k) + sigma
            c(k) = sigma * a(k, k)
            d(k) = -scale * sigma
            DO j = k+1, n
                sum = 0.0D0
                DO i = k, m
                    sum = sum + a(i, k) * a(i, j)
                END DO
                tau = sum / c(k)
                DO i = k, m
                    a(i, j) = a(i, j) - tau * a(i, k)
                END DO
            END DO
        END IF
    END DO
END SUBROUTINE qrdcmp_rect

SUBROUTINE construct_vander(x, k, n, M, V)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, M
    DOUBLE PRECISION, INTENT(IN) :: x(n), k(n)
    DOUBLE PRECISION, INTENT(OUT) :: V(n, M)
    INTEGER :: i, j

    DO j = M, 1, -1
        DO i = 1, n
            V(i, j) = (x(i)**(j-1)) * k(i)
        END DO
    END DO
END SUBROUTINE construct_vander

SUBROUTINE qrdcmp(a,n,np,c,d,sing)
    INTEGER n,np
    DOUBLE PRECISION a(np,np),c(n),d(n)
    LOGICAL sing
    ! Constructs the QR decomposition of a(1:n,1:n), with physical dimension np. The upper
    ! triangular matrix R is returned in the upper triangle of a, except for the diagonal elements
    ! of R which are returned in d(1:n). The orthogonal matrix Q is represented as a product of
    ! n − 1 Householder matrices Q1 ... Qn-1, where Qj = 1 - uj ⊗ uj /cj . The ith component
    ! of uj is zero for i = 1, ..., j - 1 while the nonzero components are returned in a(i,j) for
    ! i = j, ..., n. sing returns as true if singularity is encountered during the decomposition,
    ! but the decomposition is still completed in this case.
    INTEGER i,j,k
    DOUBLE PRECISION scale,sigma,sum,tau
    sing = .false.
    do k = 1, n-1
        scale = 0.d0
        do i = k, n
            scale = max(scale, abs(a(i,k)))
        enddo
        if (scale .eq. 0.) then !! Singular case.
            sing = .true.
            c(k) = 0.d0
            d(k) = 0.d0
            else !! Form Qk and Qk·A.
            do i = k, n
                a(i,k) = a(i,k)/scale
            enddo
            sum = 0.d0
            do i = k, n
                sum = sum + a(i,k)**2
            enddo
            sigma  = sign(sqrt(sum),a(k,k))
            a(k,k) = a(k,k)+sigma
            c(k)   = sigma*a(k,k)
            d(k)   = -scale*sigma
            do j = k+1, n
                sum = 0.d0
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

SUBROUTINE qrupdt(r,qt,n,np,u,v)
    INTEGER n,np
    DOUBLE PRECISION r(np,np),qt(np,np),u(np),v(np)
    ! USES rotate
    ! Given the QR decomposition of some n × n matrix, calculates the QR decomposition of
    ! the matrix Q·(R + u⊗v). The matrices r and qt have physical dimension np. Note that
    ! QT is input and returned in qt.
    INTEGER i,j,k
    do k = n, 1, -1 !! Find largest k such that u(k) .ne. 0.
        if (u(k).ne.0.) goto 1
    enddo
    k=1
1   do i=k-1,1,-1       !! Transform R + u⊗v to upper Hes-
        call rotate(r,qt,n,np,i,u(i),-u(i+1)) !! senberg.
        if (u(i).eq.0.) then
            u(i) = abs(u(i+1))
        else if (abs(u(i)).gt.abs(u(i+1))) then
            u(i)=abs(u(i))*sqrt(1.d0+(u(i+1)/u(i))**2)
        else
            u(i)=abs(u(i+1))*sqrt(1.d0+(u(i)/u(i+1))**2)
        endif
    enddo
    do j = 1, n
        r(1,j) = r(1,j) + u(1)*v(j)
    enddo
    do i = 1, k-1                   !! Transform upper Hessenberg matrix
        call rotate(r,qt,n,np,i,r(i,i),-r(i+1,i)) !! to upper triangular.
    enddo
    return
END SUBROUTINE

SUBROUTINE rotate(r,qt,n,np,i,a,b)
    INTEGER n,np,i
    DOUBLE PRECISION a,b,r(np,np),qt(np,np)
    ! Given n×n matrices r and qt of physical dimension np, carry out a Jacobi rotation on rows i
    ! and i + 1 of each matrix. a and b are the parameters of the rotation: cos θ = a/√a2 + b2,
    ! sin θ = b/√a2 + b2.
    INTEGER j
    DOUBLE PRECISION c,fact,s,w,y
    if(a.eq.0.)then !! Avoid unnecessary overflow or underflow.
        c = 0.
        s = sign(1.d0,b)
    else if (abs(a).gt.abs(b)) then
        fact = b/a
        c = sign(1.d0/sqrt(1.d0+fact**2),a)
        s = fact*c
    else
        fact = a/b
        s    = sign(1.d0/sqrt(1.d0+fact**2),b)
        c    = fact*s
    endif
    do j = i, n ! Premultiply r by Jacobi rotation.
        y = r(i,j)
        w = r(i+1,j)
        r(i,j)   = c*y - s*w
        r(i+1,j) = s*y + c*w
    enddo
    do j = 1, n !! Premultiply qt by Jacobi rotation.
        y=qt(i,j)
        w=qt(i+1,j)
        qt(i,j)=c*y-s*w
        qt(i+1,j)=s*y+c*w
    enddo
    return
END