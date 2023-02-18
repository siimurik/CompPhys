/*
===========================================================
 Compile and execute with:
    $ gcc -o sym sym_tridiag.c -lgsl
    $ ./sym
===========================================================
 On Mac, use these commands to find what are the '-I' flags
 and '-L' libraries, use the commands:
    $ gsl-config --cflags
    $ gsl-config --libs
 Using this, you can compile by using:
    $ gcc sym_tridiag.c -o sym -I/opt/homebrew/Cellar/gsl/2.7.1/include -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas
 ===========================================================
 This function solves the general N-by-N system A x = b
 where A is symmetric tridiagonal (N \geq 2). The off-
 diagonal vector e must be one element shorter than the
 diagonal vector diag. The form of A for the 4-by-4 case
 is shown below
    A = [   d0  e0  0   0
            e0  d1  e1  0
            0   e1  d2  e2
            0   0   e2  d3  ]
===========================================================
*/
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_linalg.h>

#define N 201
int main(void)
{
    // Create the matrix
    gsl_vector *d = gsl_vector_alloc (N);
    gsl_vector *e = gsl_vector_alloc (N-1);
    gsl_vector *x = gsl_vector_alloc (N);
    gsl_vector *b = gsl_vector_alloc (N);
    static double X[N];
    int i, j;
    struct timespec start, stop;

    printf("Dimension of tridiagonal system: %d x %d.\n", N-1, N-1);

    // Define step size
    double h = 2.0/(double)(N-1);

    // A vector of 200 elements spanning
    // from 0.0 to 2.0
    for (int i = 0; i < N; i++)
    {
        X[i] = i*h;
    }

    for (i = 0; i < N-1; i++) {
        gsl_vector_set(e, i, 0.0);
    }

    for (i = 0; i < N; i++) {
        gsl_vector_set(b, i, 0.0);
        gsl_vector_set(d, i, 0.0);
    }

    // Boundary conditions
    gsl_vector_set(d, 0, 1.0/h + h/2.0*(4.0-X[0]));
    gsl_vector_set(b, 0, h/2.0*(X[0]+5.0) - 1.0);

    // Elements of the main diagonal D and
    // righthand side vector B
    for (i = 1; i < N-1; i++) {
        gsl_vector_set(d, i, 2.0/h + h*(4.0-X[i]));
        gsl_vector_set(b, i, h*(X[i]+5.0));
    }

    // Elements of off-diagonal elements
    for (i = 0; i < N-1; i++) {
        gsl_vector_set(e, i, -1.0/h);
    }

    // Final values
    gsl_vector_set(d, N-1, 1.0/h + h/2.0*(4.0-X[N-1]));
    gsl_vector_set(b, N-1, h/2.0*(X[N-1]+5.0) - 1.0);

    // Print out vectors
    printf("\nLower and upper diagonal elements: \ni:\te:\n");
    for (int i = N-6; i < N-1; i++)
    {
        printf("%d\t%lf\n ", i, gsl_vector_get(e,i));
    }

    printf("\nMain diagonal elements: \ni:\td:\n");
    for (int i = N-5; i < N; i++)
    {
        printf("%d\t%lf\n ", i, gsl_vector_get(d,i));
    }

    // Get starting time
    clock_gettime(CLOCK_MONOTONIC, &start);

    gsl_linalg_solve_symm_tridiag(d, e, b, x);

    // Record the ending time
    clock_gettime(CLOCK_MONOTONIC, &stop);
    double time_taken = (stop.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (stop.tv_nsec - start.tv_nsec)) * 1e-9;

    // Write solutions into separate file
    FILE *fp;
    fp = fopen("output.dat", "w");
    for (int i = 0; i < N; i++)
    {
        fprintf(fp, "%g \t%g\n", X[i], gsl_vector_get(x, i));
    }
    fclose(fp);

    // Print out first five elements of the solution
    printf ("\nSolution: \ni \t u \n");
    for (int i = N-5; i < N; i++)
        printf ("%d \t %g\n", i, gsl_vector_get (x, i));

    // Print out elapsed time
    printf("\nElapsed time: %lf seconds.\n", time_taken);

    // Free allocated memory
    gsl_vector_free (x);
    gsl_vector_free (b);
    gsl_vector_free (d);
    gsl_vector_free (e);

    return 0;
}
