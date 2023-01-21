/*
========================================================
 Compile and execute with:
    $ gcc -o lus lu_solver.c  -lgsl -lgslcblas
    $ ./lus
========================================================
*/
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_linalg.h>
#define n 201
int main(void)
{
    // Matrix A has dimensions MxN
    int N = n-1;

    // Declare arrays to hold the matrix A and vector B
    double h = 2.0/(double)(N);
    static double A[n*n];
    static double X[n];
    static double y[n];

    for (int i = 0; i < n; i++)
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
        A[i*n + i-1] = -1.0/h;
        A[i*n + i]   =  2.0/h + h*(4-X[i]);
        A[i*n + i+1] = -1.0/h;
        y[i]     = h*(X[i]+5);
    }

    //A[N][N-1]   = -1.0/h;
    //A[N][N] = 1.0/h + h/2*(4-X[N]);
    A[N*n + N-1]   = -1.0/h;
    A[N*n + N] = 1.0/h + h/2*(4-X[N]);
    y[N]     = h/2*(X[N]+5) - 1;

/*
    // Set the values of A and B
    // by writing matrix values into vectors
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            AA[i*n + j] = A[i][j];
        }
        B[i] = y[i];
    }
*/

    // Print out the values of A
    printf("\nLower right elements of matrix A:\n");
    for (int i = n-5; i < n; i++)
    {
        for (int j = n-5; j < n; j++)
        {
            //printf("%lf ", A[i][j]);
            printf("%lf ", A[i*n + j]);

        }
        printf("\n");
    }
    printf("\nX:\t\ty:\n");
    for (int i = n-5; i < n; i++)
    {
        printf("%lf\t%lf\n ", X[i], y[i]);
    }


    gsl_matrix_view m = gsl_matrix_view_array (A, n, n);
    gsl_vector_view b = gsl_vector_view_array (y, n);
    gsl_vector *x = gsl_vector_alloc (n);

    int s;

    gsl_permutation * p = gsl_permutation_alloc (n);
    gsl_linalg_LU_decomp (&m.matrix, p, &s);

    // Get starting time
    clock_t start = clock();

    gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

    // Get end time
    clock_t stop = clock();

    // Write solutions into separate file
    FILE *fp;
    fp = fopen("output.dat", "w");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%g \t%g\n", X[i], gsl_vector_get(x, i));
    }
    fclose(fp);

    // Print out first five elements of the solution
    printf ("\nSolution: \ni \t u \n");
    //gsl_vector_fprintf (stdout, x, "%g");
    for (int i = 0; i < 5; i++)
        //printf ("element %d = %g\n", i, gsl_vector_get (x, i));
        printf ("%d \t %g\n", i, gsl_vector_get (x, i));

    // Calculate the elapsed time in seconds
    double elapsed_time = (double)(stop - start) / CLOCKS_PER_SEC;
    printf("\nElapsed time: %g sec.\n", elapsed_time);

    gsl_permutation_free (p);
    gsl_vector_free (x);

    return 0;
}
