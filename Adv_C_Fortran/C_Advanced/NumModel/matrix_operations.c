/* Compile and execute with
    $ gcc -o matrix_operations matrix_operations.c -lgsl -lm
    $ ./matrix_operations
*/
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <math.h>

typedef gsl_matrix Matrix;

void print_matrix(Matrix *m, const char *name) {
    size_t i, j;
    printf("%s =\n", name);
    for (i = 0; i < m->size1; i++) {
        for (j = 0; j < m->size2; j++) {
            printf("%7.3f ", gsl_matrix_get(m, i, j));
        }
        printf("\n");
    }
    printf("\n");
}

Matrix* matrix_alloc(size_t rows, size_t cols) {
    return gsl_matrix_alloc(rows, cols);
}

void matrix_free(Matrix *m) {
    gsl_matrix_free(m);
}

void matrix_set_from_array(Matrix *m, const double *data) {
    size_t i, j, idx = 0;
    for (i = 0; i < m->size1; i++) {
        for (j = 0; j < m->size2; j++) {
            gsl_matrix_set(m, i, j, data[idx++]);
        }
    }
}

void matrix_add(Matrix *A, Matrix *B, Matrix *result) {
    gsl_matrix_add(A, B);
    gsl_matrix_memcpy(result, A);
}

void matrix_sub(Matrix *A, Matrix *B, Matrix *result) {
    gsl_matrix_sub(A, B);
    gsl_matrix_memcpy(result, A);
}

void matrix_elementwise_mult(Matrix *A, Matrix *B, Matrix *result) {
    size_t i, j;
    for (i = 0; i < A->size1; i++) {
        for (j = 0; j < A->size2; j++) {
            double val = gsl_matrix_get(A, i, j) * gsl_matrix_get(B, i, j);
            gsl_matrix_set(result, i, j, val);
        }
    }
}

void matrix_elementwise_div(Matrix *A, Matrix *B, Matrix *result) {
    size_t i, j;
    for (i = 0; i < A->size1; i++) {
        for (j = 0; j < A->size2; j++) {
            double val = gsl_matrix_get(A, i, j) / gsl_matrix_get(B, i, j);
            gsl_matrix_set(result, i, j, val);
        }
    }
}

void matrix_elementwise_pow(Matrix *A, Matrix *B, Matrix *result) {
    size_t i, j;
    for (i = 0; i < A->size1; i++) {
        for (j = 0; j < A->size2; j++) {
            double val = pow(gsl_matrix_get(A, i, j), gsl_matrix_get(B, i, j));
            gsl_matrix_set(result, i, j, val);
        }
    }
}

Matrix* matrix_inverse(Matrix *A) {
    int s;
    Matrix *inverse = matrix_alloc(A->size1, A->size2);
    gsl_matrix *LU = gsl_matrix_alloc(A->size1, A->size2);
    gsl_matrix_memcpy(LU, A);
    
    gsl_permutation *p = gsl_permutation_alloc(A->size1);
    gsl_linalg_LU_decomp(LU, p, &s);
    gsl_linalg_LU_invert(LU, p, inverse);

    gsl_permutation_free(p);
    gsl_matrix_free(LU);
    
    return inverse;
}

Matrix* matrix_transpose(Matrix *A) {
    Matrix *transpose = matrix_alloc(A->size2, A->size1);
    gsl_matrix_transpose_memcpy(transpose, A);
    return transpose;
}

double matrix_determinant(Matrix *A) {
    int s;
    gsl_matrix *LU = gsl_matrix_alloc(A->size1, A->size2);
    gsl_matrix_memcpy(LU, A);
    gsl_permutation *p = gsl_permutation_alloc(A->size1);

    gsl_linalg_LU_decomp(LU, p, &s);
    double det = gsl_linalg_LU_det(LU, s);

    gsl_permutation_free(p);
    gsl_matrix_free(LU);
    return det;
}

void matrix_multiply(Matrix *A, Matrix *B, Matrix *result) {
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, B, 0.0, result);
}

int main() {
    printf("Ülesanne 14.\n");

    Matrix *A = matrix_alloc(2, 4);
    Matrix *B = matrix_alloc(2, 4);

    // Initializing matrix A
    // Initialize matrix A using a static array
    const double A_data[] = {
        1.0, -1.0, 3.0, 0.0,
        2.0,  0.0, -1.0, 2.0
    };
    matrix_set_from_array(A, A_data);

    // Initializing matrix B
    const double B_data[] = {
        4.0, 2.0, -1.0, 2.0,
        3.0, 1.0,  5.0, 0.0
    };
    matrix_set_from_array(B, B_data);

    print_matrix(A, "Matrix A");
    print_matrix(B, "Matrix B");

    // Matrix inverse
    if (A->size1 == A->size2) {
        Matrix *A_inv = matrix_inverse(A);
        print_matrix(A_inv, "Inverse of A");
        matrix_free(A_inv);
    }

    // Transpose
    Matrix *A_trans = matrix_transpose(A);
    print_matrix(A_trans, "Transpose of A");
    matrix_free(A_trans);

    // Determinant
    if (A->size1 == A->size2) {
        double det = matrix_determinant(A);
        printf("Determinant of A = %f\n\n", det);
    }

    // Matrix addition
    Matrix *A_plus_B = matrix_alloc(A->size1, A->size2);
    matrix_add(A, B, A_plus_B);
    print_matrix(A_plus_B, "A + B");

    // Matrix subtraction
    Matrix *A_minus_B = matrix_alloc(A->size1, A->size2);
    matrix_sub(A, B, A_minus_B);
    print_matrix(A_minus_B, "A - B");

    // Element-wise multiplication
    Matrix *A_mult_B = matrix_alloc(A->size1, A->size2);
    matrix_elementwise_mult(A, B, A_mult_B);
    print_matrix(A_mult_B, "A .* B");

    // Matrix multiplication (Requires A to have compatible dimensions)
    if (A->size2 == B->size1) {
        Matrix *A_mat_B = matrix_alloc(A->size1, B->size2);
        matrix_multiply(A, B, A_mat_B);
        print_matrix(A_mat_B, "A * B");
        matrix_free(A_mat_B);
    }

    // Element-wise division
    Matrix *A_div_B = matrix_alloc(A->size1, A->size2);
    matrix_elementwise_div(A, B, A_div_B);
    print_matrix(A_div_B, "A ./ B");

    // Element-wise power
    Matrix *A_pow_B = matrix_alloc(A->size1, A->size2);
    matrix_elementwise_pow(A, B, A_pow_B);
    print_matrix(A_pow_B, "A .^ B");

    // Clean up
    matrix_free(A);
    matrix_free(B);
    matrix_free(A_plus_B);
    matrix_free(A_minus_B);
    matrix_free(A_mult_B);
    matrix_free(A_div_B);
    matrix_free(A_pow_B);

    return 0;
}
