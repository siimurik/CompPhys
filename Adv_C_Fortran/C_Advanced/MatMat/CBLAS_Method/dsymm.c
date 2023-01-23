/*
=========================================================
 This code needs the BLAS library to be installed:
    $ sudo apt install libblas-dev
=========================================================
 Compile and execute
    $ gcc dsymm.c -o sym -lblas
    $ ./sym
    or
    $ gcc dsymm.c -O2 -lblas
    $ ./a.out
=========================================================
 This currently the fastest possible matrix multiplcation 
 alogrithm I have written, that is used for symmetric 
 matrices.
=========================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <time.h>

#define dim 5000
int main()
{
    // Define the dimensions of the matrices
    int m = dim, n = dim, k = dim, lda = dim, ldb = dim, ldc = dim;
    double alpha = 1.0, beta = 0.0;
    printf("Dimension of matrix matrix multiplicaton: %d x %d.\n", m,n);

    // Allocate memory for the matrices
    //double *A = malloc(m * k * sizeof(double)); // 1000x1000 matrix
    //double *B = malloc(k * n * sizeof(double)); // 1000x1000 matrix
    //double *C = malloc(m * n * sizeof(double)); // 1000x1000 matrix
    
    // Static approach
    static double A[dim][dim];
    static double B[dim][dim];
    static double C[dim][dim];

    // fill the matrices with some values
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < k; j++)
        {
            A[i][j] = ((double)rand()/(double)(RAND_MAX));
            B[i][j] = ((double)rand()/(double)(RAND_MAX));
	        //A[i * m + j] = ((double)rand()/(double)(RAND_MAX));
	        //B[i * n + j] = ((double)rand()/(double)(RAND_MAX));
	    }
    }

    // print the result (only the first few elements)
    printf("\nPrinting out the first 5x5 elements of matrix A:\n");
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            printf("%lf ", A[i][j]);
            //printf("%lf ", A[i*dim + j]);
        }
        printf("\n");
    }
    printf("\nPrinting out the first 5x5 elements of matrix B:\n");
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            printf("%lf ", B[i][j]);
            //printf("%lf ", B[i*dim + j]);
        }
        printf("\n");
    }

    // Get starting time
    clock_t tic = clock();

    // perform matrix multiplication using dsymm()
    // For dynamic allocation
    //cblas_dsymm(CblasRowMajor,CblasLeft ,CblasUpper,m,n,alpha, A, lda, B, ldb, beta, C, ldc);
    
    // For static allocation
    cblas_dsymm(CblasRowMajor,CblasLeft ,CblasUpper,m,n,alpha, &A[0][0], lda, &B[0][0], ldb, beta, &C[0][0], ldc);
    // Get end time
    clock_t toc = clock();

    // print the result (only the first few elements)
    printf("\nPrinting out the first 5x5 elements of the result matrix:\n");
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            printf("%lf ", C[i][j]);
            //printf("%lf ", C[i*dim + j]);
        }
        printf("\n");
    }
    
    // Calculate the elapsed time in seconds
    double elapsed_time = (double)(toc - tic) / CLOCKS_PER_SEC;
    printf("\nMatrix multiplication took %lf seconds.\n", elapsed_time);

    // free the allocated memory if 
    // dynamically allocated
    //free(A);
    //free(B);
    //free(C);

    return 0;
}
