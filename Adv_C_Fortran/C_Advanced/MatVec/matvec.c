/* Compile with:
 	$ gcc matvec.c -o matvec -lcblas -lm
 	$ ./matvec
====================================================
 Matrix-Vector multiplication using CBLAS and 
 dynamically allocated arrays. So this baby is supa
 fast and optimised.
----------------------------------------------------
 Solves the following system:
    {x} * [A] = {y}
----------------------------------------------------
 Input:
 m - Size of Column ( the number of rows ) 
 n -  Size of Row ( the number of columns ) 
 lda - Leading dimension of n * m matrix is n 
 A - matrix, with dimensions n*m 
 x - vector x with n elements
----------------------------------------------------
 Output:
 Resulting vector y with dimension n.
====================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cblas.h"
#include <time.h>

int main ( )
{
    CBLAS_LAYOUT Layout;
    CBLAS_TRANSPOSE transa;

    double *a, *x, *y;
    double alpha, beta;
    int m, n, lda, incx, incy, i;

    Layout = CblasColMajor;
    transa = CblasNoTrans;

    m = 3; /* Size of Column ( the number of rows ) */
    n = 3; /* Size of Row ( the number of columns ) */
    lda = 3; /* Leading dimension of 3 * 3 matrix is 3 */
    incx = 1;
    incy = 1;
    alpha = 1;
    beta = 0;

    a = (double *)malloc(sizeof(double)*m*n);
    x = (double *)malloc(sizeof(double)*n);
    y = (double *)malloc(sizeof(double)*n);

    static double values[] = {  0.1, 0.2, 0.7, 
                                0.4, 0.5, 0.1, 
                                0.0, 0.1, 0.9   }; // a 3x3 matrix

    //Copy values to dynamically allocated matrix
    memcpy(a, values, sizeof(double) * m * n);

    /* The elements of x*/
    x[0] = 0.7;
    x[1] = 0.1;
    x[2] = 0.2;

    /*Result vector, which is initially empty*/
    y[0] = 0;
    y[1] = 0;
    y[2] = 0;

    // Start the clock
    clock_t start = clock();

    // Perform matrix-vector multiplication using the BLAS library
    cblas_dgemv( Layout, transa, m, n, alpha, a, lda, x, 
                    incx, beta, y, incy );
    
    // Stop the clock
    clock_t end = clock();

    // Calculate the elapsed time in seconds
    double elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;

    /* Print y */
    for( i = 0; i < n; i++ )
        printf(" y%d = %f\n", i, y[i]);

    // Print the elapsed time
    printf("Calculations took %5.4e seconds.\n", elapsed_time);
    
    /*Freeing allocated memory*/
    free(a);
    free(x);
    free(y);
    return 0;
}
