/*
 Compile and execute with:
    $ gcc dgesv_Csolve.c -o dg -llapacke
    $ ./dg
*/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <lapacke.h>

#define N 5001
#define NRHS 1
#define LDA N
#define LDB N
int main(void)
{
    int n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
    int ipiv[N];
    //int N = N-1;
    double h = 2.0/(double)(N-1);
    //static double A[N*N];
    //static double X[N];
    //static double y[N];

    double *A = malloc(lda * n * sizeof(double));
    double *X = malloc(ldb * nrhs * sizeof(double));
    double *y = malloc(ldb * nrhs * sizeof(double));


    for (int i = 0; i < N; i++)
    {
        X[i] = i*h;
    }
    //A[0][0] = 1.0/h + h/2*(4-X[0]);
    //A[0][1] = -1.0/h;
    A[0] = 1.0/h + h/2*(4-X[0]);
    A[1] = -1.0/h;
    y[0] = h/2*(X[0]+5) - 1;

    for (int i = 1; i < N; i++)
    {
        //A[i][i-1] = -1.0/h;
        //A[i][i]   =  2.0/h + h*(4-X[i]);
        //A[i][i+1] = -1.0/h;
        A[i*N + i-1] = -1.0/h;
        A[i*N + i]   =  2.0/h + h*(4-X[i]);
        A[i*N + i+1] = -1.0/h;
        y[i]     = h*(X[i]+5);
    }

    //A[N][N-1]   = -1.0/h;
    //A[N][N] = 1.0/h + h/2*(4-X[N]);
    A[(N-1)*N + (N-2)]   = -1.0/h;
    A[(N-1)*N + (N-1)] = 1.0/h + h/2*(4-X[N-1]);
    y[N-1]     = h/2*(X[N-1]+5) - 1;

    // Print out the values of A, y, and u
    printf("\nLower right elements of matrix A:\n");
    for (int i = N-5; i < N; i++)
    {
        for (int j = N-5; j < N; j++)
        {
            //printf("%lf ", A[i][j]);
            printf("%lf ", A[i*N + j]);

        }
        printf("\n");
    }
    printf("\nX:\t\ty:\n");
    for (int i = N-5; i < N; i++)
    {
        printf("%lf\t%lf\n ", X[i], y[i]);
    }


    // Record the start time
    clock_t start = clock();

    // Call LAPACKE_dgesv
    info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, A, lda, ipiv, y, ldb );

    // Record the end time
    clock_t end = clock();

    /* Check for the exact singularity */
    if( info > 0 ) {
        printf( "The diagonal element of the triangular factor of A,\n" );
        printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
        printf( "the solution could not be computed.\n" );
        exit( 1 );
    }

    // Calculate the elapsed time in seconds
    double elapsed_time = (end - start) / (double)CLOCKS_PER_SEC;

    // printing out the first and last 3 elements
    printf("\nDGESV in C Program Results\n");
    for (int i = 0; i < 3; i++) {
        printf("%f %f\n", X[i], y[i]);
    }
    printf("\t...\n");

    for (int i = N-3; i < N; i++) {
        printf("%f %f\n", X[i], y[i]);
    }

    printf("\nElapsed time: %e s\n", elapsed_time);
    printf("\n");

    // Free the memory that was allocated for M and V

    free(A);
    free(X);
    free(y);


    return 0;
}
