#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double rand_range(double min, double max);
double **random_matrix(int rows, int cols);
void free_matrix(double **matrix, int rows);
void print_matrix(double **matrix, int rows, int cols);
void iterativeMatmulLoopReorderSubMatrix(
        double A[][K],
        double B[][N],
        double C[][N],
        const int M,
        const int K,
        const int loM,
        const int hiM,
        const int loN,
        const int hiN,
        const int loK,
        const int hiK); 
void mtMatmul(double A[][K], double B[][N], double C[][N], int M, int N, int K);

int main() {
    int rows = 3;
    int cols = 3;
    int M = 8;

    // Generate a random matrix with the specified dimensions.
    //double **random_mat = random_matrix(rows, cols);
    double **matA = random_matrix(rows, cols);
    double **matB = random_matrix(rows, cols);
    double **matC = random_matrix(rows, cols);

    if (random_mat == NULL) {
        return 1;
    }

    printf("Random matrix:\n");
    // Print the elements of the random matrix.
    //print_matrix(random_mat, rows, cols);
    print_matrix(matC, rows, cols);

    // Free the memory allocated for the random matrix.
    //free_matrix(random_mat, rows);
    free_matrix(matA, rows);
    free_matrix(matB, rows);
    free_matrix(matC, rows);
    return 0;
}

// Generate a random floating-point number within the given range.
double rand_range(double min, double max) {
    return min + (max-min) * ((double)rand()/RAND_MAX);
}

// Generate a random matrix with the specified number of rows and columns.
double **random_matrix(int rows, int cols) {
    // Allocate memory for an array of pointers to double arrays (rows).
    double **matrix = (double **)malloc(rows * sizeof(double *));

    if (matrix == NULL) {
        perror("Memory allocation failed");
        return NULL;
    }

    // Allocate memory for each row of the matrix.
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double *)malloc(cols * sizeof(double));

        if (matrix[i] == NULL) { // Check matrix[i], not matrix
            perror("Memory allocation failed");
            return NULL;
        }
    }
    
    // Seed the random number generator based on the current time.
    srand((unsigned int)time(NULL));

    // Fill the matrix with random numbers.
    for (int i = 0; i < rows; i++) { // Use i, not j
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = rand_range(0.0, 1.0);
        }
    }

    return matrix;
}

// Free the memory allocated for the matrix.
void free_matrix(double **matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);    // Free memory for each row.
    }
    free(matrix);           // Free memory for the array of pointers 
}

// Function for printing a 2D matrix
void print_matrix(double **matrix, int rows, int cols) {
    int max_print_size = 6;

    printf("[\n");

    if (rows <= max_print_size && cols <= max_print_size) {
        for (int i = 0; i < rows; i++) {
            printf("  [");
            for (int j = 0; j < cols; j++) {
                printf("%f", matrix[i][j]);
                if (j < cols - 1) printf(", ");
            }
            printf("]");
            if (i < rows - 1) printf(",\n");
        }
    } else {
        for (int i = 0; i < max_print_size; i++) {
            printf("  [");
            for (int j = 0; j < max_print_size; j++) {
                printf("%f", matrix[i][j]);
                if (j < max_print_size - 1) printf(", ");
            }
            printf(", ...");
            printf("]");
            if (i < max_print_size - 1) printf(",\n");
        }
        printf(",\n  ...\n");
        printf("  ...\n");
        printf("  ...");
    }

    printf("\n]\n");
}

void iterativeMatmulLoopReorderSubMatrix(
        double A[][K],
        double B[][N],
        double C[][N],
        const int M,
        const int K,
        const int loM,
        const int hiM,
        const int loN,
        const int hiN,
        const int loK,
        const int hiK) {

    for (int m = loM; m < hiM; m++) {
        for (int k = loK; k < hiK; k++) {
            for (int n = loN; n < hiN; n++) {
                C[m][n] += A[m][k] * B[k][n];
            }
        }
    }
}

void mtMatmul(double A[][K], double B[][N], double C[][N], int M, int N, int K) {
    const int m2 = M / 2;
    const int n2 = N / 2;
    const int k2 = K / 2;

    // First part of C11: A_11*B_11
    iterativeMatmulLoopReorderSubMatrix(A, B, C, M, K, 0, m2, 0, n2, 0, k2);

    // First part of C12: A_11*B_12
    iterativeMatmulLoopReorderSubMatrix(A, B, C, M, K, 0, m2, n2 + 1, N, 0, k2);

    // First part of C21: A_21*B_11
    iterativeMatmulLoopReorderSubMatrix(A, B, C, M, K, m2 + 1, M, 0, n2, 0, k2);

    // First part of C22: A_21*B_12
    iterativeMatmulLoopReorderSubMatrix(A, B, C, M, K, m2 + 1, M, n2 + 1, N, 0, k2);

    // Second part of C11: A_12*B_21
    iterativeMatmulLoopReorderSubMatrix(A, B, C, M, K, 0, m2, 0, n2, k2 + 1, K);

    // Second part of C12: A_12*B_22
    iterativeMatmulLoopReorderSubMatrix(A, B, C, M, K, 0, m2, n2 + 1, N, k2 + 1, K);

    // Second part of C21: A_22*B_21
    iterativeMatmulLoopReorderSubMatrix(A, B, C, M, K, m2 + 1, M, 0, n2, k2 + 1, K);

    // Second part of C22: A_22*B_22
    iterativeMatmulLoopReorderSubMatrix(A, B, C, M, K, m2 + 1, M, n2 + 1, N, k2 + 1, K);
}
