/* Compile and execute with
    $ gcc -o linear_solver linear_solver.c -lgsl -lm
    $ ./linear_solver
*/
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

typedef gsl_matrix Matrix;
typedef gsl_vector Vector;

void print_vector(Vector *v, const char *name) {
    size_t i;
    printf("%s =\n", name);
    for (i = 0; i < v->size; i++) {
        printf("%7.4f\n", gsl_vector_get(v, i));
    }
    printf("\n");
}

Matrix* matrix_alloc(size_t rows, size_t cols) {
    return gsl_matrix_alloc(rows, cols);
}

Vector* vector_alloc(size_t size) {
    return gsl_vector_alloc(size);
}

void matrix_free(Matrix *m) {
    gsl_matrix_free(m);
}

void vector_free(Vector *v) {
    gsl_vector_free(v);
}

void matrix_set_from_array(Matrix *m, const double *data) {
    size_t i, j, idx = 0;
    for (i = 0; i < m->size1; i++) {
        for (j = 0; j < m->size2; j++) {
            gsl_matrix_set(m, i, j, data[idx++]);
        }
    }
}

void vector_set_from_array(Vector *v, const double *data) {
    for (size_t i = 0; i < v->size; i++) {
        gsl_vector_set(v, i, data[i]);
    }
}

Vector* solve_linear_system(Matrix *A, Vector *b) {
    int s;
    Matrix *LU = matrix_alloc(A->size1, A->size2);
    gsl_matrix_memcpy(LU, A);  // Copy A to preserve the original matrix

    gsl_permutation *p = gsl_permutation_alloc(A->size1);
    gsl_linalg_LU_decomp(LU, p, &s);

    Vector *x = vector_alloc(b->size);
    gsl_linalg_LU_solve(LU, p, b, x);  // Solve Ax = b

    gsl_permutation_free(p);
    matrix_free(LU);

    return x;
}

int main() {
    printf("Ülesanne 10.\n");

    // Define matrix A
    Matrix *A = matrix_alloc(3, 3);
    const double A_data[] = {
        67, 12, 16,
        -54, 71, 43,
        -82, 35, 88
    };
    matrix_set_from_array(A, A_data);

    // Define vector b
    Vector *b = vector_alloc(3);
    const double b_data[] = {31, 15, -14};
    vector_set_from_array(b, b_data);

    // Solve the system A * x = b
    Vector *x = solve_linear_system(A, b);

    // Print solution vector x
    print_vector(x, "Solution x");

    // Clean up
    matrix_free(A);
    vector_free(b);
    vector_free(x);

    return 0;
}
