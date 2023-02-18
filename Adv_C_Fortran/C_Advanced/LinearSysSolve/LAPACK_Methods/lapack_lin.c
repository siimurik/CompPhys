/*
======================================================================================================
 Compile and execute with:
    $ gcc lapack_lin.c -o lap -llapacke
    $ ./lap
---    
 On Mac:   
    $ gcc laplin.c -o lin -I/opt/homebrew/opt/lapack/include -L/opt/homebrew/opt/lapack/lib -llapacke
======================================================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>
#include <time.h>

#define N 100
#define NRHS 1
#define LDA N
#define LDB NRHS

int main (){
    int n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
    int ipiv[N];

    // Dynamically allocate memory for M and V
    double *M = malloc(LDA * N * sizeof(double));
    double *V = malloc(LDB * N * sizeof(double));

    // Generate random elements for M and V
    srand(time(NULL));
    for (int i = 0; i < LDA * N; i++) {
        M[i] = (double)rand() / RAND_MAX;
    }
    for (int i = 0; i < LDB * N; i++) {
        V[i] = (double)rand() / RAND_MAX;
    }

    // Record the start time
    clock_t start = clock();

    // Call LAPACKE_dgesv
    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, M, lda, ipiv, V, ldb);

    // Record the end time
    clock_t end = clock();

    // Calculate the elapsed time in seconds
    double elapsed_time = (end - start) / (double)CLOCKS_PER_SEC;

    printf("Solution to the system of linear equations:");
    for (int i = 0; i < N; i++) {
        printf("\nx[%d] = %f", i, V[i]);
    }
    printf("\nElapsed time: %e s\n", elapsed_time);
    printf("\n");

    // Free the memory that was allocated for M and V
    free(M);
    free(V);

    return 0;
}
