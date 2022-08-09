!* Fortran source code is found in dgemm_example.f
! Compile and execute with:
!   $ gfortran -o matmul_lapack matmul_lapack.f90 -L/usr/local/lib -llapack -lblas
!   $ ./matmul_lapack
!
PROGRAM   MAIN

    IMPLICIT NONE

    REAL(8)                 :: ALPHA, BETA, T1, T2, TIME
    INTEGER(4)              :: I, J
    INTEGER(4), PARAMETER   :: M=2000, K=200, N=1000
    REAL(8), DIMENSION(M,K) :: A
    REAL(8), DIMENSION(K,N) :: B
    REAL(8), DIMENSION(M,N) :: C

    PRINT *, "This example computes real matrix C=alpha*A*B+beta*C"
    PRINT *, "using LAPACK package function dgemm, where A, B, and C"
    PRINT *, "are matrices and alpha and beta are double precision "
    PRINT *, "scalars"
    PRINT *, ""

    PRINT *, "Initializing data for matrix multiplication C=A*B for "
    PRINT 10, " matrix A(",M," x",K, ") and matrix B(", K," x", N, ")"
10      FORMAT(a,I5,a,I5,a,I5,a,I5,a)
    PRINT *, ""
    ALPHA = 1.0 
    BETA = 0.0

    PRINT *, "Intializing matrix data"
    PRINT *, ""
    DO I = 1, M
    DO J = 1, K
        A(I,J) = (I-1) * K + J
    END DO
    END DO

    DO I = 1, K
    DO J = 1, N
        B(I,J) = -((I-1) * N + J)
    END DO
    END DO

    DO I = 1, M
    DO J = 1, N
        C(I,J) = 0.0
    END DO
    END DO

    PRINT *, "Computing matrix product using LAPACK DGEMM "
    PRINT *, "subroutine"
    CALL CPU_TIME(T1)
    CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M)
    CALL CPU_TIME(T2)
    TIME = T2 - T1
    PRINT *, "Computations completed."
    WRITE (*,15) TIME
    PRINT *, ""
    
15      FORMAT(/'Calculation time is ', e9.4, ' seconds.')

    PRINT *, "Top left corner of matrix A:"
    PRINT 20, ((A(I,J), J = 1,MIN(K,6)), I = 1,MIN(M,6))
    PRINT *, ""

    PRINT *, "Top left corner of matrix B:"
    PRINT 20, ((B(I,J),J = 1,MIN(N,6)), I = 1,MIN(K,6))
    PRINT *, ""

20      FORMAT(6(F12.0,1x))

    PRINT *, "Top left corner of matrix C:"
    PRINT 30, ((C(I,J), J = 1,MIN(N,6)), I = 1,MIN(M,6))
    PRINT *, ""

30      FORMAT(6(ES12.4,1x))

    PRINT *, "Example completed."
    STOP 

END PROGRAM