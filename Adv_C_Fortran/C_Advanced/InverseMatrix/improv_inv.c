#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

// Typedefs for better readability
typedef gsl_matrix Matrix;
typedef gsl_permutation Permutation;

// Function declarations
Matrix* create_matrix(int size);
void fill_matrix_random(Matrix* matrix, int size);
void print_matrix(const char* title, Matrix* matrix, int size);
int decompose_matrix(Matrix* matrix, Permutation* p, int* signum);
int invert_matrix(Matrix* matrix, Permutation* p, Matrix* inverse);
void cleanup(Matrix* A, Matrix* Ainverse, Permutation* p);

int main() {
    int N = 4;
    int signum;
    int status;

    Matrix* A = create_matrix(N);
    Matrix* Ainverse = create_matrix(N);
    Permutation* p = gsl_permutation_alloc(N);

    fill_matrix_random(A, N);

    print_matrix("Initial Matrix", A, N);

    status = decompose_matrix(A, p, &signum);
    printf("Status of decomposition: %d\n", status);

    status = invert_matrix(A, p, Ainverse);
    printf("Status of inversion: %d\n", status);

    print_matrix("Inverse Matrix", Ainverse, N);

    cleanup(A, Ainverse, p);

    return 0;
}

// Function to create a matrix
Matrix* create_matrix(int size) {
    return gsl_matrix_alloc(size, size);
}

// Function to fill the matrix with random values
void fill_matrix_random(Matrix* matrix, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            gsl_matrix_set(matrix, i, j, drand48());
        }
    }
}

// Function to print a matrix
void print_matrix(const char* title, Matrix* matrix, int size) {
    printf("%s\n", title);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            printf("%16.6f ", gsl_matrix_get(matrix, i, j));
        }
        printf("\n");
    }
    printf("\n");
}

// Function to perform LU decomposition on the matrix
int decompose_matrix(Matrix* matrix, Permutation* p, int* signum) {
    return gsl_linalg_LU_decomp(matrix, p, signum);
}

// Function to invert the matrix using the LU decomposition
int invert_matrix(Matrix* matrix, Permutation* p, Matrix* inverse) {
    return gsl_linalg_LU_invert(matrix, p, inverse);
}

// Function to free the allocated memory
void cleanup(Matrix* A, Matrix* Ainverse, Permutation* p) {
    gsl_matrix_free(A);
    gsl_matrix_free(Ainverse);
    gsl_permutation_free(p);
}
