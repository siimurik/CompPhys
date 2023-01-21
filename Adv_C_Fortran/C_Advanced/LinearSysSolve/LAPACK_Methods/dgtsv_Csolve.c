/*
========================================================
 Compile and execute with:
    $ gcc dgtsv_Csolve.c -o diag -llapacke
    $ ./diag
========================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <lapacke.h>

#define N 5001
#define NRHS 1
#define LDA N
#define LDB N

int main(void) {
    int n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
    double h;
    static double dl[N-1], d[N], du[N-1], a[LDA*N], b[LDB*NRHS];
    double start_time, end_time, elapsed_time;
    static double x[N];
    int i;

    // Initializing
    h = 2.0/((double)N-1);
    for (i = 0; i < N; i++) {
        x[i] = i*h;
    }
    for (i = 0; i < N-1; i++) {
        dl[i] = 0.0;
        d[i]  = 0.0;
        du[i] = 0.0;
    }
    for (i = 0; i < LDB*NRHS; i++) {
        b[i] = 0.0;
    }

    // Boundary conditions
    d[0] =  1.0/h + h/2.0*(4.0-x[0]);
    b[0] = h/2.0*(x[0]+5.0) - 1.0;

    // Elements of the main diagonal D and
    // righthand side vector B
    for (i = 1; i < N-1; i++) {
        d[i] =  2.0/h + h*(4.0-x[i]);
        b[i] = h*(x[i]+5.0);
    }

    // Elements of off-diagonal elements
    for (i = 0; i < N-1; i++) {
        dl[i] = -1.0/h;
        du[i] = -1.0/h;
    }

    // Final values
    d[N-1] = 1.0/h + h/2.0*(4.0-x[N-1]);
    b[N-1] = h/2.0*(x[N-1]+5.0) - 1.0;

    // record the starting time
    start_time = clock();

    // Solve the equations A*X = B.
    info = LAPACKE_dgtsv(LAPACK_COL_MAJOR, n, nrhs, dl, d, du, b, ldb);

    // record the ending time
    end_time = clock();

    // checking for problems
    if (info > 0) {
        printf("The diagonal element of the triangular factor of A,\n");
        printf("U(%d,%d) is zero, so that A is singular; the solution could not be computed.\n", info, info);
        return 1;
    }

    // printing out the first and last 3 elements
    printf("DGTSV in C Program Results\n");
    for (i = 0; i < 3; i++) {
        printf("%f %f\n", x[i], b[i]);
    }
    printf("            ...\n");

    for (i = N-2; i < N; i++) {
    printf("%f %f\n", x[i], b[i]);
    }

    // Print the solution into a separate file
    FILE *file = fopen("out_2.csv", "w");
    for (i = 0; i < N; i++) {
    fprintf(file, "%f %f\n", x[i], b[i]);
    }
    fclose(file);

    // compute the elapsed time
    elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // print the elapsed time
    printf("Elapsed time: %lf seconds.\n", elapsed_time);

    return 0;
}

