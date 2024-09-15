/* Compile and execute with
    $ $ gcc -o bicgstab_solver bicgstab_solver.c -lgsl -lm
    $ ./bicgstab_solver
*/
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <math.h>

typedef gsl_matrix Matrix;
typedef gsl_vector Vector;

void print_vector(Vector *v, const char *name) {
    printf("%s =\n", name);
    for (size_t i = 0; i < v->size; i++) {
        printf("%7.4f\n", gsl_vector_get(v, i));
    }
    printf("\n");
}

// Function to compute the matrix-vector product
void matrix_vector_product(const Matrix *A, const Vector *v, Vector *result) {
    gsl_blas_dgemv(CblasNoTrans, 1.0, A, v, 0.0, result);
}

// Custom BiCGSTAB solver
int bicgstab(const Matrix *A, const Vector *b, Vector *x, double tol, int max_iters, double *error, int *total_iters) {
    size_t n = b->size;
    Vector *r = gsl_vector_alloc(n);
    Vector *r_hat0 = gsl_vector_alloc(n);
    Vector *p = gsl_vector_alloc(n);
    Vector *v = gsl_vector_alloc(n);
    Vector *s = gsl_vector_alloc(n);
    Vector *t = gsl_vector_alloc(n);
    Vector *temp = gsl_vector_alloc(n);

    double rho_old = 1.0, alpha = 1.0, omega = 1.0;
    double rho_new, beta, tau;

    // r = b - A * x
    matrix_vector_product(A, x, temp);  // temp = A * x
    gsl_vector_memcpy(r, b);            // r = b
    gsl_vector_sub(r, temp);            // r = b - A * x
    gsl_vector_memcpy(r_hat0, r);       // r_hat0 = r

    double norm_b = gsl_blas_dnrm2(b);  // norm(b)
    double zeta = gsl_blas_dnrm2(r);    // norm(r)
    double errtol = tol * norm_b;
    *error = zeta;

    int k = 0;
    while (zeta > errtol && k < max_iters) {
        gsl_blas_ddot(r_hat0, r, &rho_new);  // rho_new = r_hat0' * r

        if (rho_new == 0) {
            printf("Bi-CGSTAB breakdown: rho_new = 0\n");
            break;
        }

        if (k > 0) {
            beta = (rho_new / rho_old) * (alpha / omega);
            gsl_vector_memcpy(temp, v);
            gsl_vector_scale(temp, omega);
            gsl_vector_sub(p, temp);
            gsl_vector_scale(p, beta);
            gsl_vector_add(p, r);  // p = r + beta * (p - omega * v)
        } else {
            gsl_vector_memcpy(p, r);  // p = r
        }

        matrix_vector_product(A, p, v);  // v = A * p
        gsl_blas_ddot(r_hat0, v, &tau);  // tau = r_hat0' * v

        if (tau == 0) {
            printf("Bi-CGSTAB breakdown: tau = 0\n");
            break;
        }

        alpha = rho_new / tau;
        gsl_vector_memcpy(s, r);       // s = r
        gsl_vector_scale(v, alpha);    // v = alpha * v
        gsl_vector_sub(s, v);          // s = r - alpha * v

        matrix_vector_product(A, s, t);  // t = A * s
        gsl_blas_ddot(t, t, &tau);       // tau = t' * t

        if (tau == 0) {
            printf("Bi-CGSTAB breakdown: t = 0\n");
            break;
        }

        gsl_blas_ddot(t, s, &omega);   // omega = t' * s / (t' * t)
        omega /= tau;

        if (omega == 0) {
            printf("Bi-CGSTAB breakdown: omega = 0\n");
            break;
        }

        gsl_vector_memcpy(temp, p);
        gsl_vector_scale(temp, alpha);     // temp = alpha * p
        gsl_vector_add(x, temp);           // x = x + alpha * p

        gsl_vector_memcpy(temp, s);
        gsl_vector_scale(temp, omega);     // temp = omega * s
        gsl_vector_add(x, temp);           // x = x + omega * s

        gsl_vector_memcpy(r, s);
        gsl_vector_scale(t, omega);
        gsl_vector_sub(r, t);              // r = s - omega * t

        zeta = gsl_blas_dnrm2(r);          // zeta = norm(r)
        rho_old = rho_new;                 // update rho_old

        k++;
    }

    *total_iters = k;

    gsl_vector_free(r);
    gsl_vector_free(r_hat0);
    gsl_vector_free(p);
    gsl_vector_free(v);
    gsl_vector_free(s);
    gsl_vector_free(t);
    gsl_vector_free(temp);

    return zeta <= errtol ? 0 : -1;
}

int main() {
    printf("Ülesanne 10 (BiCGSTAB).\n");

    // Define matrix A
    Matrix *A = gsl_matrix_alloc(3, 3);
    const double A_data[] = {
        67, 12, 16,
        -54, 71, 43,
        -82, 35, 88
    };
    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 3; j++)
            gsl_matrix_set(A, i, j, A_data[i * 3 + j]);

    // Define vector b
    Vector *b = gsl_vector_alloc(3);
    const double b_data[] = {31, 15, -14};
    for (size_t i = 0; i < 3; i++)
        gsl_vector_set(b, i, b_data[i]);

    // Initial guess x0 (set to zero vector)
    Vector *x = gsl_vector_alloc(3);
    gsl_vector_set_zero(x);

    // Parameters for BiCGSTAB
    double tol = 1e-12;
    int max_iters = 100;
    double error;
    int total_iters;

    // Solve using BiCGSTAB
    int status = bicgstab(A, b, x, tol, max_iters, &error, &total_iters);

    // Print results
    if (status == 0) {
        printf("BiCGSTAB converged in %d iterations.\n", total_iters);
        print_vector(x, "Solution x");
    } else {
        printf("BiCGSTAB did not converge within %d iterations.\n", max_iters);
    }

    // Clean up
    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(x);

    return 0;
}
