/*
==============================================================================
Best possible command for fast results:
    $ gcc -o gsl_openblas gsl_matmul.c -lgsl -lopenblas -lpthread -lm
Help with compiling commands
    $ gcc -I/usr/include/gsl -lgsl -lgslcblas -lm gsl_matmul.c -o gsl_matmul
    $ ./a.out
    // Result of 5000x5000 multiplication
      time = 669.75406 sec. 
==============================================================================                            
*/

#include <time.h>
#include <stdio.h>
#include <gsl/gsl_blas.h>

int
main (void)
{   
    struct timespec start, stop;
    int n = 5000;

    printf ("Dimension of matrices: %d x %d.\n", n, n);
    int size = n*n;
    double *aq = malloc (sizeof (double) * size);
    double *bq = malloc (sizeof (double) * size);
    double *cq = malloc (sizeof (double) * size);

    int i,j;

    for ( i=0;i<n;i++){
        for ( j=0;j<n;j++){
            aq[i+j*n] = 4.0 *(i+9+j);
            bq[i+j*n] = 4.0 *(11-i-2*j);
        }
    }

    gsl_matrix_view A = gsl_matrix_view_array(aq, n, n);
    gsl_matrix_view B = gsl_matrix_view_array(bq, n, n);
    gsl_matrix_view C = gsl_matrix_view_array(cq, n, n);

    // Set the number of threads for OpenBLAS
    // Actually not even necessary to explicidly write
    //openblas_set_num_threads(8); // Set the desired number of threads

    /* Compute C = A B */
    
    // Get starting time
    clock_gettime(CLOCK_MONOTONIC, &start);

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &A.matrix, &B.matrix,
                  0.0, &C.matrix);
    
    // Get end time
    clock_gettime(CLOCK_MONOTONIC, &stop);

    // Calculate the elapsed time in seconds
    double time_taken = (stop.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (stop.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("\nMatrix multiplication took %lf seconds.\n", time_taken);

    return 0;
}
