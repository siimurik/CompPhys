#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

typedef gsl_matrix Matrix;
typedef gsl_vector Vector;

void print_vector(Vector *v, const char *name) {
    printf("%s =\n", name);
    for (size_t i = 0; i < v->size; i++) {
        printf("%7.3f\n", gsl_vector_get(v, i));
    }
    printf("\n");
}

// Matrix-vector product function for Ax
void matrix_vector_product(const Matrix *A, const Vector *x, Vector *result) {
    gsl_blas_dgemv(CblasNoTrans, 1.0, A, x, 0.0, result);
}

// TFQMR solver function
int tfqmr(const Matrix *A, const Vector *b, Vector *x, double tol, int max_iters, double *error, int *total_iters) {
    size_t n = b->size;

    // Allocate memory for necessary vectors
    Vector *r = gsl_vector_alloc(n);
    Vector *w = gsl_vector_alloc(n);
    Vector *y1 = gsl_vector_alloc(n);
    Vector *y2 = gsl_vector_alloc(n);
    Vector *v = gsl_vector_alloc(n);
    Vector *u1 = gsl_vector_alloc(n);
    Vector *u2 = gsl_vector_alloc(n);
    Vector *d = gsl_vector_alloc(n);
    Vector *temp = gsl_vector_alloc(n);

    double tau, rho, rho_new, sigma, alpha, theta = 0, eta = 0, c = 0;
    double norm_r0, tau_sq;
    double norm_b = gsl_blas_dnrm2(b);  // norm(b)
    double errtol = tol * norm_b;       // error tolerance

    // Initializations
    gsl_vector_set_zero(x);  // Set initial guess x to 0
    matrix_vector_product(A, x, temp);  // temp = A * x
    gsl_vector_memcpy(r, b);            // r = b - A * x (since x is 0, r = b)
    gsl_vector_sub(r, temp);            // r = b - A * x
    gsl_vector_memcpy(w, r);            // w = r
    gsl_vector_memcpy(y1, r);           // y1 = r

    matrix_vector_product(A, y1, v);    // v = A * y1
    gsl_vector_memcpy(u1, v);           // u1 = v

    norm_r0 = gsl_blas_dnrm2(r);        // norm of initial residual
    tau = norm_r0;                      // initial tau
    rho = tau * tau;                    // rho = tau^2
    tau_sq = tau * tau;
    gsl_vector_set_zero(d);             // d = 0

    *error = tau;
    int k = 0;

    // Main iteration loop
    while (k < max_iters) {
        k++;

        // Step 1: Compute sigma = r' * v
        gsl_blas_ddot(r, v, &sigma);

        if (sigma == 0) {
            printf("TFQMR breakdown: sigma = 0\n");
            break;
        }

        alpha = rho / sigma;

        for (int j = 0; j < 2; j++) {
            int m = 2 * (k - 1) + j + 1;

            if (j == 1) {
                // y2 = y1 - alpha * v
                gsl_vector_memcpy(y2, y1);
                gsl_vector_scale(v, alpha);
                gsl_vector_sub(y2, v);

                // u2 = A * y2
                matrix_vector_product(A, y2, u2);
            }

            Vector *u = (j == 0) ? u1 : u2;
            Vector *y = (j == 0) ? y1 : y2;

            // Update w = w - alpha * u
            gsl_vector_memcpy(temp, u);
            gsl_vector_scale(temp, alpha);
            gsl_vector_sub(w, temp);

            // Update d = y + (theta^2 * eta / alpha) * d
            gsl_vector_memcpy(temp, y);
            gsl_vector_scale(d, (theta * theta * eta / alpha));
            gsl_vector_add(d, temp);

            // theta = norm(w) / tau
            double norm_w = gsl_blas_dnrm2(w);
            theta = norm_w / tau;

            // c = 1 / sqrt(1 + theta * theta)
            c = 1.0 / sqrt(1.0 + theta * theta);

            // tau = tau * theta * c
            tau = tau * theta * c;

            // eta = c * c * alpha
            eta = c * c * alpha;

            // x = x + eta * d
            gsl_vector_memcpy(temp, d);
            gsl_vector_scale(temp, eta);
            gsl_vector_add(x, temp);

            // Check for convergence
            if (tau * sqrt(m + 1) <= errtol) {
                *total_iters = k;
                gsl_vector_free(r);
                gsl_vector_free(w);
                gsl_vector_free(y1);
                gsl_vector_free(y2);
                gsl_vector_free(v);
                gsl_vector_free(u1);
                gsl_vector_free(u2);
                gsl_vector_free(d);
                gsl_vector_free(temp);
                return 0;  // Converged
            }
        }

        // rho_new = r' * w
        gsl_blas_ddot(r, w, &rho_new);

        if (rho_new == 0) {
            printf("TFQMR breakdown: rho_new = 0\n");
            break;
        }

        // beta = rho_new / rho
        double beta = rho_new / rho;
        rho = rho_new;

        // y1 = w + beta * y2
        gsl_vector_memcpy(y1, w);
        gsl_vector_memcpy(temp, y2);
        gsl_vector_scale(temp, beta);
        gsl_vector_add(y1, temp);

        // u1 = A * y1
        matrix_vector_product(A, y1, u1);

        // v = u1 + beta * (u2 + beta * v)
        gsl_vector_memcpy(v, u2);
        gsl_vector_scale(v, beta);
        gsl_vector_add(v, u1);
        gsl_vector_scale(u2, beta);
        gsl_vector_add(v, u2);

        *error = tau;
    }

    *total_iters = k;

    // Free allocated memory
    gsl_vector_free(r);
    gsl_vector_free(w);
    gsl_vector_free(y1);
    gsl_vector_free(y2);
    gsl_vector_free(v);
    gsl_vector_free(u1);
    gsl_vector_free(u2);
    gsl_vector_free(d);
    gsl_vector_free(temp);

    return -1;  // Did not converge within max_iters
}

int main() {
    printf("TFQMR Solver for Linear Systems.\n");

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

    // Parameters for TFQMR
    double tol = 1e-6;
    int max_iters = 100;
    double error;
    int total_iters;

    // Solve using TFQMR
    int status = tfqmr(A, b, x, tol, max_iters, &error, &total_iters);

    // Print results
    if (status == 0) {
        printf("TFQMR converged in %d iterations.\n", total_iters);
        print_vector(x, "Solution x");
    } else {
        printf("TFQMR did not converge within %d iterations.\n", max_iters);
    }

    // Clean up
    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(x);

    return 0;
}
