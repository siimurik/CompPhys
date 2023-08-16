/*
 Compile and excute with:
    $ gcc det.c -o det -llapacke
    $ ./det 
*/

#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>
#define dim 5

int main() {
    int n = dim; // matrix is 3x3
    int lda = n;
    //double A[9] = {3, -1, 2, -3, 3, -1, 6, 0, 4};
    double *A = (double *)malloc(n * n * sizeof(double)); // Input matrix (column-major order)
    double initValues[dim][dim] = {
        {0.154930, 0.682167, 0.270557, 0.330188, 0.440037},
        {0.590628, 0.945112, 0.548013, 0.103298, 0.519449},
        {0.312663, 0.116302, 0.847965, 0.853670, 0.358802},
        {0.363340, 0.902817, 0.776792, 0.069543, 0.617297},
        {0.848290, 0.762891, 0.093645, 0.752094, 0.088619}
    };

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i * n + j] = initValues[i][j];
        }
    }
    //int ipiv[n];
    int *ipiv = (int *)malloc(n * sizeof(int)); // Pivot array
    int info;

    // Print matrix A
    for (int i = 0; i < n; i++) {
        printf("  [");
        for (int j = 0; j < n; j++) {
            printf("%f", A[i * n + j]);
            if (j < n - 1) printf(", ");
        }
        printf("]");
        if (i < n - 1) printf(",\n");
    }

    // Perform LU-decomposition to get ipiv values
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A, lda, ipiv);
    if (info < 0) {
        printf("\nError: Argument %d had an illegal value\n", -info);
        return -1;
    } else if (info > 0) {
        printf("\nError: Matrix is singular\n");
        return -1;
    }

    // Calculate the determinant
    double determinant = 1.0;
    for (int i = 0; i < n; i++) {
        determinant *= A[i * n + i];
    }
    int sign = 1;
    for (int i = 0; i < n; i++) {
        if (ipiv[i] != i + 1) {
            sign = -sign;
        }
    }
    determinant *= sign;

    // Print the result
    printf("\nDeterminant: %lf\n", determinant);

    // Free memory
    free(A);
    free(ipiv);

    return 0;
}
