/*
 Compile and execute with:
    $ gcc laplin.c -o lin -llapacke
    $ ./lin
*/

#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>
#include <time.h>

#define N 3
#define NRHS 1
#define LDA N
#define LDB NRHS

int main (){
    int n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
    int ipiv[N];

    // Define a system of linear equations
    // 3x + 2y - z = 2
    // 2x - 2y + 4z = -3
    // -x + y/2 - z = z/2
    double M[LDA * N] = {
        3.0,  2.0, -1.0, 
        2.0, -2.0,  4.0, 
        -1.0,  0.5, -1.0
    };

    double V[LDB * N] = {
        2.0, 
        -3.0, 
        0.5
    };

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
    return 0;
}
