/*
 Compile and execute with:
    > gcc cblas_matmul_vec.c -o vecform -lblas
    > ./vecform
*/
#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <time.h>

int main() {
    // Define the dimensions of the matrices
    int m = 1000, n = 1000, k = 1000;

    // Allocate memory for the matrices
    double *A = malloc(m * k * sizeof(double)); // 1000x1000 matrix
    double *B = malloc(k * n * sizeof(double)); // 1000x1000 matrix
    double *C = malloc(m * n * sizeof(double)); // 1000x1000 matrix

    // Initialize the matrices with random values
    for (int i = 0; i < m * k; i++) {
        A[i] = (double)rand() / RAND_MAX;
    }
    for (int i = 0; i < k * n; i++) {
        B[i] = (double)rand() / RAND_MAX;
    }

    // Start the clock
     clock_t start = clock();
    // Perform the matrix multiplication
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, A, k, B, n, 0.0, C, n);
    // Stop the clock
    clock_t end = clock();

    // Calculate the elapsed time in seconds
    double elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;

    // Print the first 10 elements of the result matrix
    for (int i = 0; i < 10; i++) {
        printf("C[%d] = %f\n", i, C[i]);
    }
    
    // Print the elapsed time
    printf("Calculation took %f seconds.\n", elapsed_time);

    return 0;
}
