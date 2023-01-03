!  =============================================================================
!   $ gfortran -o dgesv dgesv_solve.f90  -llapack
!   $ ./dgsev
! In GNU Plot:
!   $ plot 'out.csv' with lines
!  =============================================================================
!
!  DGESV Example.
!  ==============
!
!  The program computes the solution to the system of linear
!  equations with a square matrix A and multiple
!  right-hand sides B, where A is the coefficient matrix:
!
!    6.80  -6.05  -0.45   8.32  -9.67
!   -2.11  -3.30   2.58   2.71  -5.14
!    5.66   5.36  -2.70   4.35  -7.26
!    5.97  -4.44   0.27  -7.17   6.08
!    8.23   1.08   9.04   2.14  -6.87
!
!  and B is the right-hand side matrix:
!
!    4.02  -1.56   9.81
!    6.19   4.00  -4.09
!   -8.22  -8.67  -4.57
!   -7.57   1.75  -8.61
!   -3.03   2.86   8.99
!
!  Description.
!  ============
!
!  The routine solves for X the system of linear equations A*X = B,
!  where A is an n-by-n matrix, the columns of matrix B are individual
!  right-hand sides, and the columns of X are the corresponding
!  solutions.
!
!  The LU decomposition with partial pivoting and row interchanges is
!  used to factor A as A = P*L*U, where P is a permutation matrix, L
!  is unit lower triangular, and U is upper triangular. The factored
!  form of A is then used to solve the system of equations A*X = B.
!
!  Example Program Results.
!  ========================
!
! DGESV Example Program Results
!
! Solution
!  -0.80  -0.39   0.96
!  -0.70  -0.55   0.22
!   0.59   0.84   1.90
!   1.32  -0.10   5.36
!   0.57   0.11   4.04
!
! Details of LU factorization
!   8.23   1.08   9.04   2.14  -6.87
!   0.83  -6.94  -7.92   6.55  -3.99
!   0.69  -0.67 -14.18   7.24  -5.19
!   0.73   0.75   0.02 -13.82  14.19
!  -0.26   0.44  -0.59  -0.34  -3.43
!
! Pivot indices
!      5      5      3      4      5
!  =============================================================================
!
!     .. Parameters ..
program dgesv_solver
    IMPLICIT NONE
    INTEGER          N, NRHS
    PARAMETER        ( N = 201, NRHS = 1 )
    INTEGER          LDA, LDB, i
    PARAMETER        ( LDA = N, LDB = N )
    INTEGER          INFO
    INTEGER          IPIV( N )
    DOUBLE PRECISION A( LDA, N ), B( LDB, NRHS )
    DOUBLE PRECISION h
    DOUBLE PRECISION X(N)
    DOUBLE PRECISION start_time, end_time, elapsed_time

    h = 2.D0/float(N-1)
    do i = 1, N
        X(i) = (i-1)*h
    end do
    A = 0.D0
    B = 0.D0

    A(1,1) =  1.D0/h + h/2.D0*(4.D0-X(1))
    A(1,2) = -1.D0/h
    B(1,1)   = h/2.D0*(X(1)+5.D0) - 1.D0

    do i = 2, N-1
        A(i,i-1) = -1.D0/h
        A(i,i)   =  2.D0/h + h*(4.D0-X(i))
        A(i,i+1) = -1.D0/h
        B(i,1)     = h*(X(i)+5.D0)
    end do
    A(N,N-1)= -1.D0/h
    A(N,N)  =  1.D0/h + h/2.D0*(4.D0-X(N))
    B(N,1)  = h/2.D0*(X(N)+5.D0) - 1.D0

    ! record the starting time
    call cpu_time(start_time)

!     Solve the equations A*X = B.
    CALL DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )

    ! record the ending time
    call cpu_time(end_time)
!     Check for the exact singularity.

    IF( INFO.GT.0 ) THEN
        WRITE(*,*)'The diagonal element of the triangular factor of A,'
        WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
        WRITE(*,*)'A is singular; the solution could not be computed.'
        STOP
    END IF

    WRITE(*,*)'DGESV Program Results'
    do i = 1, 3
        write(*,'(6f12.5)') X(i), B(i,1)
    end do
    WRITE(*,*) '            ...'
    do i = N-2, N
        write(*,'(6f12.5)') X(i), B(i,1)
    end do

!     Print solution.
    !CALL PRINT_MATRIX( 'Solution', N, NRHS, B, LDB )
    open(unit = 11, file="out.csv")
    do i = 1, N
        write(11,*) X(i), B(i,1)
    end do
    close(11)

    ! compute the elapsed time
    elapsed_time = end_time - start_time

    ! print the elapsed time
    write(*,9) elapsed_time
9   format ( 'Elapsed time: ', F8.6, ' seconds.')

!     Print details of LU factorization.
    !CALL PRINT_MATRIX( 'Details of LU factorization', N, N, A, LDA )

!     Print pivot indices.
    !CALL PRINT_INT_VECTOR( 'Pivot indices', N, IPIV )
end program
!     End of DGESV Example.

!  =============================================================================
!     Auxiliary routine: printing a matrix.

    SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
    CHARACTER*(*)    DESC
    INTEGER          M, N, LDA
    DOUBLE PRECISION A( LDA, * )
    INTEGER          I, J

    WRITE(*,*)
    WRITE(*,*) DESC
    DO I = 1, M
        WRITE(*,9998) ( A( I, J ), J = 1, N )
    END DO

9998 FORMAT( 11(:,1X,F6.2) )
    RETURN
    END

!     Auxiliary routine: printing a vector of integers.

    SUBROUTINE PRINT_INT_VECTOR( DESC, N, A )
    CHARACTER*(*)    DESC
    INTEGER          N
    INTEGER          A( N )
    INTEGER          I

    WRITE(*,*)
    WRITE(*,*) DESC
    WRITE(*,9999) ( A( I ), I = 1, N )
9999 FORMAT( 11(:,1X,I6) )
    RETURN
    END
