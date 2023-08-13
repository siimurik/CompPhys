!---------------------------------------------------------------------------------
! Compile and execute with:
!   $ gfortran -o matmul_lapack matmul_lapack.f90 -L/usr/local/lib -llapack -lblas
! Or
!   $ ifx -o matmul_lapack matmul_lapack.f90 -llapack -lblas
!   $ ./matmul_lapack
!---------------------------------------------------------------------------------
! For multithreading uncomment the line
!   CALL OPENBLAS_SET_NUM_THREADS(8)
! and compile with
!   $ gfortran matmul_lapack.f90 -o matmul_lapack -lopenblas -lpthread
!---------------------------------------------------------------------------------
PROGRAM   MAIN

    IMPLICIT NONE

    DOUBLE PRECISION    :: ALPHA, BETA
    INTEGER             :: I, J
    !INTEGER(4), PARAMETER   :: M=2000, K=200, N=1000
    INTEGER, PARAMETER  :: M=5000, K=5000, N=5000
    DOUBLE PRECISION, DIMENSION(M,K) :: A
    DOUBLE PRECISION, DIMENSION(K,N) :: B
    DOUBLE PRECISION, DIMENSION(M,N) :: C
    INTEGER :: START_TIME, END_TIME, ELAPSED_TIME, RATE
    REAL :: ELAPSED_SECONDS

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
    !DO I = 1, M
    !    DO J = 1, K
    !        A(I,J) = (I-1) * K + J
    !    END DO
    !END DO

    CALL RANDOM_NUMBER(A)

    !DO I = 1, K
    !    DO J = 1, N
    !        B(I,J) = -((I-1) * N + J)
    !    END DO
    !END DO

    CALL RANDOM_NUMBER(B)

    C = 0.0
    
    ! Specify the number of cores to be used
    CALL OPENBLAS_SET_NUM_THREADS(8)

    PRINT *, "Computing matrix product using LAPACK DGEMM "
    PRINT *, "subroutine"

    CALL SYSTEM_CLOCK(count=START_TIME, count_rate=RATE)

    CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M)
    
    CALL SYSTEM_CLOCK(count=END_TIME)

    ELAPSED_TIME = END_TIME - START_TIME
    ELAPSED_SECONDS = REAL(ELAPSED_TIME) / REAL(RATE)

    PRINT *, "Computations completed."
    WRITE (*,15) ELAPSED_SECONDS
    PRINT *, ""
    
15      FORMAT(/'Calculation time is', F6.3, ' seconds.')
    !PRINT *, "Elapsed time:", ELAPSED_SECONDS, "seconds"

    PRINT *, "Top left corner of matrix A:"
    PRINT 20, ((A(I,J), J = 1,MIN(K,6)), I = 1,MIN(M,6))
    PRINT *, ""

    PRINT *, "Top left corner of matrix B:"
    PRINT 20, ((B(I,J),J = 1,MIN(N,6)), I = 1,MIN(K,6))
    PRINT *, ""

20      FORMAT(6(F12.6,1x))

    PRINT *, "Top left corner of matrix C:"
    PRINT 30, ((C(I,J), J = 1,MIN(N,6)), I = 1,MIN(M,6))
    PRINT *, ""

30      FORMAT(6(ES12.4,1x))

    PRINT *, "Example completed."
    STOP 

END PROGRAM
