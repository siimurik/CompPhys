/*
 Optimized matrix multiplication with BLAS and OpenMP
 Compile with:
    gcc -O3 -march=native -fopenmp optimized_matmul.c -o matmul -lopenblas -lm
 On Mac:
    clang -O3 -march=native -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp optimized_matmul.c -o matmul -lopenblas -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <cblas.h>

#define DIM 5000
#define DISPLAY_SIZE 5
#define CLOCK_MONOTONIC 1

typedef struct {
    double *data;  // Contiguous memory block
    int rows;
    int cols;
} Matrix;

// Function prototypes
Matrix matrix_create(int rows, int cols);
Matrix matrix_random(int rows, int cols);
void matrix_free(Matrix mat);
void matrix_print(Matrix mat, const char *title);
void matrix_multiply_blas(const Matrix *A, const Matrix *B, Matrix *C);
void parallel_initialize(Matrix *mat);

int main() {
    struct timespec start, end;
    double elapsed;
    
    // Create matrices
    Matrix A = matrix_create(DIM, DIM);
    Matrix B = matrix_create(DIM, DIM);
    Matrix C = matrix_create(DIM, DIM);
    
    // Parallel initialization
    printf("Initializing matrices...\n");
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    #pragma omp parallel sections
    {
        #pragma omp section
        parallel_initialize(&A);
        
        #pragma omp section
        parallel_initialize(&B);
    }
    
    // Initialize C to zero in parallel
    #pragma omp parallel for
    for (int i = 0; i < DIM*DIM; i++) {
        C.data[i] = 0.0;
    }
    
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Initialization completed in %.3f seconds\n\n", elapsed);
    
    // Print sample of matrices
    matrix_print(A, "Matrix A");
    matrix_print(B, "Matrix B");
    
    // Matrix multiplication
    printf("\nMultiplying matrices (%dx%d) × (%dx%d)...\n", DIM, DIM, DIM, DIM);
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    matrix_multiply_blas(&A, &B, &C);
    
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("\nMultiplication completed in %.3f seconds\n", elapsed);
    
    // Print sample of result
    matrix_print(C, "\nResult matrix C");
    
    // Free memory
    matrix_free(A);
    matrix_free(B);
    matrix_free(C);
    
    return 0;
}

// Create a matrix with contiguous memory
Matrix matrix_create(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (double *)aligned_alloc(64, rows * cols * sizeof(double));
    if (mat.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    return mat;
}

// Parallel matrix initialization
void parallel_initialize(Matrix *mat) {
    unsigned int seed = time(NULL) + omp_get_thread_num();
    
    #pragma omp parallel for
    for (int i = 0; i < mat->rows * mat->cols; i++) {
        mat->data[i] = (double)rand_r(&seed) / RAND_MAX;
    }
}

// Free matrix memory
void matrix_free(Matrix mat) {
    free(mat.data);
}

// Print a sample of the matrix
void matrix_print(Matrix mat, const char *title) {
    printf("%s:\n", title);
    
    int print_size = (mat.rows > DISPLAY_SIZE || mat.cols > DISPLAY_SIZE) ? 
                     DISPLAY_SIZE : mat.rows;
    
    for (int i = 0; i < print_size; i++) {
        for (int j = 0; j < print_size; j++) {
            printf("%8.4f ", mat.data[i * mat.cols + j]);
        }
        if (mat.cols > print_size) printf("...");
        printf("\n");
    }
    if (mat.rows > print_size) printf(" ...\n");
    printf("\n");
}

// Optimized matrix multiplication using BLAS
void matrix_multiply_blas(const Matrix *A, const Matrix *B, Matrix *C) {
    // Use column-major order for maximum BLAS performance
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                A->rows, B->cols, A->cols,
                1.0, A->data, A->rows,
                B->data, B->rows,
                0.0, C->data, C->rows);
}