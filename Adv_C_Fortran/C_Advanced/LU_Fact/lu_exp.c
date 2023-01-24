#include <stdio.h>
#include <lapacke.h>

int main() {
    int N = 3; // matrix is 3x3
    int LDA = N;
    double A[9] = {3, -1, 2, -3, 3, -1, 6, 0, 4};
    int IPIV[N];
    int INFO;


    INFO = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N, N, A, LDA, IPIV);
    if (INFO < 0) {
        printf("Error: Argument %d had an illegal value\n", -INFO);
        return -1;
    } else if (INFO > 0) {
        printf("Error: Matrix is singular\n");
        return -1;
    }

    printf("The LU factorization of A is: \n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%.4f ", A[i * N + j]);
        }
        printf("\n");
    }

    printf("The L matrix is: \n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i > j) {
                printf("%.4f ", A[i * N + j]);
            } else if (i == j) {
                printf("1.0000 ");
            } else {
                printf("0.0000 ");
            }
        }
        printf("\n");
    }
    //printing the U matrix
    printf("The U matrix is: \n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i <= j) {
                printf("%.4f ", A[i * N + j]);
            } else {
                printf("0.0000 ");
            }
        }
        printf("\n");
    }


    return 0;
}
