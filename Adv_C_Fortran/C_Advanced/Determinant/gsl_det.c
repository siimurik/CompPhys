/* Compile and execute with:
    $ gcc gsl_det.c -o gsl_det -lgsl
    $ ./gsl_det <rows> <cols>

Example:
    $./gsl_det 3 3 
Enter the elements of the matrix (3x3):
5 8 3
6 1 0 
5 9 2
Determinant: 61
*/
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

typedef gsl_matrix Matrix;

void read_matrix(Matrix *matrix, int rows, int cols);
double compute_determinant(Matrix *matrix);

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <rows> <cols>\n", argv[0]);
        return EXIT_FAILURE;
    }

    int rows = atoi(argv[1]);
    int cols = atoi(argv[2]);

    if (rows != cols) {
        fprintf(stderr, "Error: Matrix must be square to compute determinant.\n");
        return EXIT_FAILURE;
    }

    Matrix *matrix = gsl_matrix_alloc(rows, cols);

    read_matrix(matrix, rows, cols);

    double determinant = compute_determinant(matrix);

    printf("Determinant: %g\n", determinant);

    gsl_matrix_free(matrix);

    return EXIT_SUCCESS;
}

// Function to read matrix from the input
void read_matrix(Matrix *matrix, int rows, int cols) {
    printf("Enter the elements of the matrix (%dx%d):\n", rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double element;
            if (scanf("%lf", &element) != 1) {
                fprintf(stderr, "Error reading element (%d, %d)\n", i, j);
                exit(EXIT_FAILURE);
            }
            gsl_matrix_set(matrix, i, j, element);
        }
    }
}

// Function to compute the determinant of a square matrix
double compute_determinant(Matrix *matrix) {
    int signum;
    gsl_permutation *p = gsl_permutation_alloc(matrix->size1);

    Matrix *lu = gsl_matrix_alloc(matrix->size1, matrix->size2);
    gsl_matrix_memcpy(lu, matrix);

    gsl_linalg_LU_decomp(lu, p, &signum);
    double det = gsl_linalg_LU_det(lu, signum);

    gsl_permutation_free(p);
    gsl_matrix_free(lu);

    return det;
}