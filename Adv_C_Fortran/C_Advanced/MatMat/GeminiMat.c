#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

// Define the matrix dimensions (N x N) and the block size.
// N=512 is a good test size. A BLOCK_SIZE of 32 or 64 often optimizes L1/L2 cache hits.
#define N 5000
#define BLOCK_SIZE 64

/**
 * @brief Initializes an NxN matrix with random double values.
 * @param M The matrix to initialize.
 * @param n The size of the matrix (N).
 */
void init_matrix(double *M, int n) {
    // Seed the random number generator only once
    static int seeded = 0;
    if (!seeded) {
        srand(time(NULL));
        seeded = 1;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // Assign a small random value for predictable testing
            M[i * n + j] = (double)rand() / (RAND_MAX / 10.0);
        }
    }
}

/**
 * @brief Performs highly optimized blocked matrix multiplication C = A * B.
 * * This function uses three nested loops over the block indices (ii, jj, kk)
 * to maximize data locality (cache reuse) for matrices A, B, and C.
 * It uses OpenMP to parallelize the block iterations.
 *
 * @param A Left-hand side matrix (N x N).
 * @param B Right-hand side matrix (N x N).
 * @param C Result matrix (N x N).
 * @param n Dimension of the matrices (N).
 * @param bs The block size (BLOCK_SIZE).
 */
void multiply_matrix_blocked(const double *A, const double *B, double *C, int n, int bs) {
    // Initialize C to zero. This is crucial as the block multiplication is an accumulation.
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i * n + j] = 0.0;
        }
    }

    // Blocked Matrix Multiplication: C(ii, jj) += A(ii, kk) * B(kk, jj)
    // Parallelize over the block row (ii) and block column (jj) indices.
    // The 'collapse(2)' clause tells OpenMP to merge the first two loops for load balancing.
    #pragma omp parallel for collapse(2) schedule(static)
    for (int ii = 0; ii < n; ii += bs) {
        for (int jj = 0; jj < n; jj += bs) {
            for (int kk = 0; kk < n; kk += bs) {
                
                // Perform multiplication for the current block C(i, j) = A(i, k) * B(k, j)
                // The innermost loops are not parallelized as they are designed to execute fast
                // within a CPU core, keeping the block data localized in the cache.
                for (int i = ii; i < ii + bs && i < n; i++) {
                    for (int k = kk; k < kk + bs && k < n; k++) {
                        // This optimization focuses on reusing A[i * n + k] across the inner j-loop.
                        // This is a crucial step for achieving high performance (register reuse).
                        register double temp_A_ik = A[i * n + k];
                        for (int j = jj; j < jj + bs && j < n; j++) {
                            C[i * n + j] += temp_A_ik * B[k * n + j];
                        }
                    }
                }
            }
        }
    }
}


/**
 * @brief Main function to run the simulation, measure time, and verify result.
 */
int main() {
    printf("--- OpenMP Blocked Matrix Multiplication ---\n");
    printf("Matrix Size: N=%d x %d\n", N, N);
    printf("Block Size: BS=%d\n", BLOCK_SIZE);
    printf("Number of Threads: %d\n", omp_get_max_threads());

    // 1. Memory Allocation (using aligned memory would offer more speed but requires platform-specific functions)
    double *A = (double *)malloc(N * N * sizeof(double));
    double *B = (double *)malloc(N * N * sizeof(double));
    double *C = (double *)malloc(N * N * sizeof(double));

    if (!A || !B || !C) {
        perror("Failed to allocate memory");
        return 1;
    }

    // 2. Initialization
    init_matrix(A, N);
    init_matrix(B, N);

    // 3. Execution and Timing
    double start_time = omp_get_wtime();

    multiply_matrix_blocked(A, B, C, N, BLOCK_SIZE);

    double end_time = omp_get_wtime();

    // 4. Output Results
    double time_taken = end_time - start_time;
    double ops = 2.0 * N * N * N; // 2 operations (1 mult + 1 add) per element
    double gflops = (ops / time_taken) / 1e9;

    printf("\nExecution Time: %.4f seconds\n", time_taken);
    printf("Performance: %.2f GFLOP/s\n", gflops);
    printf("--------------------------------------------\n");
    // Example output of one element for simple verification
    // printf("C[0][0] = %.4f\n", C[0]);

    // 5. Cleanup
    free(A);
    free(B);
    free(C);

    return 0;
}