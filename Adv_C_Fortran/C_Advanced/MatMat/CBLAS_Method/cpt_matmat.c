/*
To perform fast matrix multiplication in C programming, 
you can use the BLAS (Basic Linear Algebra Subprograms) library. 
This library provides optimized routines for performing common 
linear algebra operations, such as matrix multiplication, 
on large matrices.

To use the BLAS library, you will first need to install it. 
Here is an example of how to do this on a Linux system:
    > sudo dnf install libblas-dev
---------------------------------------------    
 Compile and execute with:
    > gcc -o cpt cpt_matmat.c -lcblas -lblas
    > ./cpt
*/
#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <time.h>

int main()
{   
    // Define two matrices
    double A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9}; // 3x3 matrix
    double B[9] = {9, 8, 7, 6, 5, 4, 3, 2, 1}; // 3x3 matrix

    // Define a result matrix
    double C[9]; // 3x3 matrix

    // Start the clock
     clock_t start = clock();
    // Perform matrix multiplication using the BLAS library
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1.0, A, 3, B, 3, 0.0, C, 3);
    // Stop the clock
    clock_t end = clock();

    // Calculate the elapsed time in seconds
    double elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;

    // Print the result matrix
    printf("Result matrix:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%f ", C[i*3 + j]);
        }
      printf("\n");
    }

    // Print the elapsed time
    printf("Calculation took %5.4e seconds.\n", elapsed_time);

    return 0;
}
