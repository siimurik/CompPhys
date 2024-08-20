/* Compile and execute with:
    $ gcc lu_solver2.c -o lu2 -lgsl -lm
    $ ./lu2 <rows> <cols>

Example:
    $ ./lu2 3 3
Insert values of matrix (3x3):
A[1][1]: 3
A[1][2]: -1
A[1][3]: 14
A[2][1]: 2
A[2][2]: 2
A[2][3]: 3
A[3][1]: 1
A[3][2]: -12
A[3][3]: -18
Insert values of vector (size 3):
b[1]: 7
b[2]: 0
b[3]: 33
Solution vector x:
x[1] = 2.53846
x[2] = -2.23325
x[3] = -0.203474

NOTE: You can also input all values at once and it will also work.
*/
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>

// Define a type for a matrix system
typedef struct {
    size_t rows;
    size_t cols;
    gsl_matrix *matrix;
    gsl_vector *vector;
    gsl_vector *solution;
} LinearSystem;

// Function prototypes
void read_matrix(LinearSystem *system);
void solve_system(LinearSystem *system);
void print_solution(const LinearSystem *system);
void free_system(LinearSystem *system);

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <rows> <cols>\n", argv[0]);
        return EXIT_FAILURE;
    }

    // Parse rows and columns from command line
    size_t rows = atoi(argv[1]);
    size_t cols = atoi(argv[2]);

    if (rows != cols) {
        fprintf(stderr, "Error: Number of rows must equal number of columns for a square matrix.\n");
        return EXIT_FAILURE;
    }

    // Initialize the linear system
    LinearSystem system;
    system.rows = rows;
    system.cols = cols;
    system.matrix = gsl_matrix_alloc(rows, cols);
    system.vector = gsl_vector_alloc(rows);
    system.solution = gsl_vector_alloc(rows);

    // Read matrix and vector from user input
    read_matrix(&system);

    // Solve the linear system
    solve_system(&system);

    // Print the solution
    print_solution(&system);

    // Free allocated memory
    free_system(&system);

    return EXIT_SUCCESS;
}

// Function to read matrix and vector elements from the user
void read_matrix(LinearSystem *system) {
    printf("Insert values of matrix (%zux%zu):\n", system->rows, system->cols);
    for (size_t i = 0; i < system->rows; i++) {
        for (size_t j = 0; j < system->cols; j++) {
            double value;
            printf("A[%zu][%zu]: ", i + 1, j + 1);
            scanf("%lf", &value);
            gsl_matrix_set(system->matrix, i, j, value);
        }
    }

    printf("Insert values of vector (size %zu):\n", system->rows);
    for (size_t i = 0; i < system->rows; i++) {
        double value;
        printf("b[%zu]: ", i + 1);
        scanf("%lf", &value);
        gsl_vector_set(system->vector, i, value);
    }
}

// Function to solve the linear system using LU decomposition
void solve_system(LinearSystem *system) {
    gsl_permutation *p = gsl_permutation_alloc(system->rows);
    int signum;

    gsl_linalg_LU_decomp(system->matrix, p, &signum);
    gsl_linalg_LU_solve(system->matrix, p, system->vector, system->solution);

    gsl_permutation_free(p);
}

// Function to print the solution vector
void print_solution(const LinearSystem *system) {
    printf("Solution vector x:\n");
    for (size_t i = 0; i < system->rows; i++) {
        printf("x[%zu] = %g\n", i + 1, gsl_vector_get(system->solution, i));
    }
}

// Function to free the allocated memory for the linear system
void free_system(LinearSystem *system) {
    gsl_matrix_free(system->matrix);
    gsl_vector_free(system->vector);
    gsl_vector_free(system->solution);
}
