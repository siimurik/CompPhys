/*
========================================================
 Compile and execute with:
    $ gcc -o lus lu_solver.c -lgsl -lgslcblas
 or
    $ gcc lu_solver.c -o lus -I/opt/homebrew/Cellar/gsl/2.7.1/include -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas
    $ ./lus
========================================================
 This is a linear system solver code that uses GSL's 
 LU-decomposition and then solves the system of linear
 equations using the decomposed matrix. This code works
 most effectively for small matrices.
========================================================
*/
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_linalg.h>

#define n 201
int main(void)
{    
    // Create the matrix
    gsl_matrix *A = gsl_matrix_alloc (n+1, n+1);
    gsl_vector *x = gsl_vector_alloc (n+1);
    gsl_vector *y = gsl_vector_alloc (n+1);
    int s;
    double h = 2.0/(double)(n-1);
    static double X[n+1];

    // A vector of 200 elements spanning
    // from 0.0 to 2.0
    for (int i = 0; i < n; i++)
    {
        X[i] = i*h;
    }

    // Setting left boundary conditions
    gsl_matrix_set(A, 0, 0, 1.0/h + h/2*(4-X[0]));
    gsl_matrix_set(A, 0, 1, -1.0/h);
    gsl_vector_set(y, 0, h/2*(X[0]+5) - 1) ;

    // Writing the diagonal elements
    for (int i = 1; i < n; i++)
    {
        gsl_matrix_set(A, i, i-1, -1.0/h);
        gsl_matrix_set(A, i, i, 2.0/h + h*(4-X[i]));
        gsl_matrix_set(A, i, i+1, -1.0/h);
        gsl_vector_set(y, i, h*(X[i]+5));
    }

    // Setting right boundary conditions
    gsl_matrix_set(A, n, n-1, -1.0/h);
    gsl_matrix_set(A, n, n, 1.0/h + h/2*(4-X[n]));
    gsl_vector_set(y, n, h/2*(X[n]+5) - 1);

    // Print the initial matrix
    printf("Initial Matrix\n");
    for (int i = n-5; i < n; i++)
        {
        for (int j = n-5; j < n; j++)
            {
                printf("%16.4f ",gsl_matrix_get(A, i, j));
            }
        printf("\n");
        }
    printf("\n");

    // Print out vectors
    printf("\nX:\t\ty:\n");
    for (int i = n-5; i < n; i++)
    {
        printf("%lf\t%lf\n ", X[i], gsl_vector_get(y,i));
    }

    // Perform LU-decomposition
    gsl_permutation * p = gsl_permutation_alloc (n+1);
    gsl_linalg_LU_decomp (A, p, &s);

    // Get starting time
    clock_t start = clock();

    gsl_linalg_LU_solve (A, p, y, x);

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
    for (int i = n-5; i < n; i++)
        //printf ("element %d = %g\n", i, gsl_vector_get (x, i));
        printf ("%d \t %g\n", i, gsl_vector_get (x, i));

    // Calculate the elapsed time in seconds
    double elapsed_time = (double)(stop - start) / CLOCKS_PER_SEC;
    printf("\nElapsed time: %g sec.\n", elapsed_time);

    // Free allocated memory
    gsl_permutation_free (p);
    gsl_vector_free (x);
    gsl_vector_free (y);
    gsl_matrix_free (A);

    return 0;
}
