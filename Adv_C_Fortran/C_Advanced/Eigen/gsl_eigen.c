/* Compile and execute with:
    $ gcc -o gsl_eigen gsl_eigen.c -lgsl -lgslcblas -lm
    $ ./gsl_eigen <rows> <cols>

Example:
    $ ./gsl_eigen 5 5
Enter the elements of the matrix (row-wise):
    -1.01   0.86  -4.60   3.31  -4.81
     3.98   0.53  -7.04   5.29   3.55
     3.30   8.26  -3.89   8.20  -1.51
     4.43   4.96  -7.66  -7.33   6.18
     7.31  -6.43  -6.16   2.47   5.58
Eigenvalues:
(   2.86,  10.76) (   2.86, -10.76) ( -10.46,   0.00) (  -0.69,  -4.70) (  -0.69,   4.70) 

Right Eigenvectors:
(   0.02,   0.20) (   0.02,  -0.20) (   0.46,   0.00) (   0.73,   0.06) (   0.73,  -0.06) 
(   0.48,  -0.05) (   0.48,   0.05) (   0.34,   0.00) (  -0.03,   0.01) (  -0.03,  -0.01) 
(   0.31,  -0.41) (   0.31,   0.41) (   0.31,   0.00) (   0.17,   0.31) (   0.17,  -0.31) 
(   0.40,   0.09) (   0.40,  -0.09) (  -0.74,   0.00) (  -0.08,   0.07) (  -0.08,  -0.07) 
(   0.48,   0.24) (   0.48,  -0.24) (   0.16,   0.00) (  -0.33,   0.47) (  -0.33,  -0.47) 

Left Eigenvectors:
(   0.11,   0.43) (   0.11,  -0.43) (  -0.02,   0.00) (  -0.01,  -0.01) (  -0.01,   0.01) 
(   0.69,  -0.18) (   0.69,   0.18) (   0.03,   0.00) (   0.03,  -0.00) (   0.03,   0.00) 
(   0.03,   0.08) (   0.03,  -0.08) (   0.02,   0.00) (  -0.14,  -0.36) (  -0.14,   0.36) 
(  -0.13,  -0.49) (  -0.13,   0.49) (  -0.01,   0.00) (  -0.75,   0.34) (  -0.75,  -0.34) 
(   0.16,   0.08) (   0.16,  -0.08) (   1.00,   0.00) (   0.40,   0.09) (   0.40,  -0.09) 

NOTE:
For symmetric (or Hermitian) matrices, the eigenvectors computed by 
GSL are orthogonal when associated with distinct eigenvalues. However, 
for non-symmetric matrices, the eigenvectors are generally not 
orthogonal. GSL does not enforce orthogonality but instead computes 
the eigenvectors based on the matrix properties. For non-symmetric 
matrices, left and right eigenvectors satisfy a biorthogonality 
condition, meaning that the left eigenvector corresponding to one 
eigenvalue is orthogonal to the right eigenvector corresponding to 
any other eigenvalue.
*/
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>

typedef struct {
    size_t rows;
    size_t cols;
} MatrixDimensions;

typedef struct {
    gsl_matrix_complex *eigenvectors;
    gsl_vector_complex *eigenvalues;
} EigenResults;

// Function declarations
MatrixDimensions parse_dimensions(int argc, char *argv[]);
gsl_matrix* create_matrix(MatrixDimensions dims);
void fill_matrix(gsl_matrix *matrix);
EigenResults compute_eigen(gsl_matrix *matrix);
gsl_matrix* transpose_matrix(const gsl_matrix *matrix);
void print_eigenvalues(const gsl_vector_complex *eigenvalues, size_t n);
void print_eigenvectors(const gsl_matrix_complex *eigenvectors, size_t n);
void print_eigen_results(const EigenResults results, const char *label);
void free_eigen_results(EigenResults results);

int main(int argc, char *argv[]) {
    // Parse matrix dimensions from command line
    MatrixDimensions dims = parse_dimensions(argc, argv);

    // Create and fill the matrix
    gsl_matrix *matrix = create_matrix(dims);
    fill_matrix(matrix);

    // Compute eigenvalues and right eigenvectors
    EigenResults right_results = compute_eigen(matrix);

    // Compute left eigenvectors by transposing the matrix
    gsl_matrix *transposed_matrix = transpose_matrix(matrix);
    EigenResults left_results = compute_eigen(transposed_matrix);

    // Print eigenvalues
    printf("Eigenvalues:\n");
    print_eigenvalues(right_results.eigenvalues, dims.rows);

    // Print right eigenvectors
    print_eigen_results(right_results, "Right Eigenvectors");

    // Print left eigenvectors
    print_eigen_results(left_results, "Left Eigenvectors");

    // Free memory
    free_eigen_results(right_results);
    free_eigen_results(left_results);
    gsl_matrix_free(matrix);
    gsl_matrix_free(transposed_matrix);

    return 0;
}

// Function to parse matrix dimensions from command line
MatrixDimensions parse_dimensions(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <rows> <cols>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    MatrixDimensions dims;
    dims.rows = (size_t) atoi(argv[1]);
    dims.cols = (size_t) atoi(argv[2]);

    return dims;
}

// Function to create a matrix of given dimensions
gsl_matrix* create_matrix(MatrixDimensions dims) {
    gsl_matrix *matrix = gsl_matrix_alloc(dims.rows, dims.cols);
    if (matrix == NULL) {
        fprintf(stderr, "Failed to allocate matrix.\n");
        exit(EXIT_FAILURE);
    }
    return matrix;
}

// Function to fill the matrix with user input
void fill_matrix(gsl_matrix *matrix) {
    printf("Enter the elements of the matrix (row-wise):\n");
    for (size_t i = 0; i < matrix->size1; ++i) {
        for (size_t j = 0; j < matrix->size2; ++j) {
            double value;
            scanf("%lf", &value);
            gsl_matrix_set(matrix, i, j, value);
        }
    }
}

// Function to transpose the matrix
gsl_matrix* transpose_matrix(const gsl_matrix *matrix) {
    gsl_matrix *transposed_matrix = gsl_matrix_alloc(matrix->size2, matrix->size1);
    gsl_matrix_transpose_memcpy(transposed_matrix, matrix);
    return transposed_matrix;
}

// Function to compute the eigenvalues and eigenvectors of the matrix
EigenResults compute_eigen(gsl_matrix *matrix) {
    EigenResults results;
    size_t n = matrix->size1;

    // Allocate workspace and results
    gsl_eigen_nonsymmv_workspace *workspace = gsl_eigen_nonsymmv_alloc(n);

    results.eigenvalues = gsl_vector_complex_alloc(n);
    results.eigenvectors = gsl_matrix_complex_alloc(n, n);

    // Compute eigenvalues and eigenvectors
    gsl_eigen_nonsymmv(matrix, results.eigenvalues, results.eigenvectors, workspace);

    // Sort the eigenvalues and eigenvectors (optional)
    gsl_eigen_nonsymmv_sort(results.eigenvalues, results.eigenvectors, GSL_EIGEN_SORT_ABS_DESC);

    gsl_eigen_nonsymmv_free(workspace);

    return results;
}

// Function to print the eigenvalues
void print_eigenvalues(const gsl_vector_complex *eigenvalues, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        gsl_complex eval = gsl_vector_complex_get(eigenvalues, i);
        printf("( %6.2f, %6.2f) ", GSL_REAL(eval), GSL_IMAG(eval));
    }
    printf("\n");
}

// Function to print the eigenvectors
void print_eigenvectors(const gsl_matrix_complex *eigenvectors, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            gsl_complex evec = gsl_matrix_complex_get(eigenvectors, i, j);
            printf("( %6.2f, %6.2f) ", GSL_REAL(evec), GSL_IMAG(evec));
        }
        printf("\n");
    }
}

// Function to print the results with a label
void print_eigen_results(const EigenResults results, const char *label) {
    printf("\n%s:\n", label);
    print_eigenvectors(results.eigenvectors, results.eigenvectors->size1);
}

// Function to free the memory used by eigenvalues and eigenvectors
void free_eigen_results(EigenResults results) {
    gsl_vector_complex_free(results.eigenvalues);
    gsl_matrix_complex_free(results.eigenvectors);
}
