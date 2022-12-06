/*
 Compile and execute with
    $ gcc lapack_matmul.c -o lapmul -llapacke -lcblas
    $ ./lapmul 
*/
#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <lapacke.h>

int main() {
    // Matrix dimensions
    int m = 2;
    int k = 3;
    int n = 4;

    // Matrix data A[m*k]
    double A[6] = {
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0
    };
    // Matrix data B[k*n]
    double B[12] = {
        7.0, 8.0, 9.0, 10.0,
        11.0, 12.0, 13.0, 14.0,
        15.0, 16.0, 17.0, 18.0
    };

    // Result matrix
    double C[m*n];

    // Perform matrix multiplication
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, 1.0, A, k, B, n, 0.0, C, n);

    // Print result matrix
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
        printf("%f ", C[i*n + j]);
        }
        printf("\n");
    }

    return 0;
}
