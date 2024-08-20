/*=============================================================================
 Compile and execute:
    $ gcc indepSchEp.c -o indepSchEp -llapacke -lm
    $ ./indepSchEp
-------------------------------------------------------------------------------
Description:
This code numerically solves the Time Independent Schrödinger Equation (TISE) 
for a confined quantum particle, which is an example of the regular Sturm-
Liouville eigenvalue problem. 

The numerical method employed involves constructing and solving a tridiagonal
matrix that approximates the differential operator in the TISE. The eigenvalues 
and eigenvectors of this matrix correspond to the energy levels and quantum 
states of the system, respectively. The code leverages LAPACK routines for 
efficient computation of these eigenvalues and eigenvectors. Additionally, 
the user can explore specific eigenvectors and ranges within them, facilitating 
detailed analysis of the quantum states.
===============================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

// Typedefs for readability
typedef double (*PotentialFunc)(double);

typedef struct {
    int size;
    double *d;  // The diagonal elements of the array.
    double *e;  // The off-diagonal elements of the array.
} TridiagonalMatrix;

// Declare functions
double mL2V(double y);
TridiagonalMatrix* init_tridiagonal_matrix(int N, double (*potential)(double));
void solve_eigenproblem(TridiagonalMatrix *matrix, double *eigenvalues, double *eigenvectors);
void print_matrix(const TridiagonalMatrix *matrix);
void print_eigenvalues(int N, double *eigenvalues);
void print_eigenvector_range(int N, double *eigenvectors, int vec_index, int start_index, int end_index);

// Main function
int main() {
    int N = 100;
    double *eigenvalues = malloc(N * sizeof(double));
    double *eigenvectors = malloc(N * N * sizeof(double));

    TridiagonalMatrix *matrix = init_tridiagonal_matrix(N, mL2V);
    printf("The initial tridiagonal matrix:\n");
    print_matrix(matrix);

    // Solve the eigenvalue problem using LAPACKE_dstev subroutine
    // NB! Values of 'matrix' are overwritten!
    solve_eigenproblem(matrix, eigenvalues, eigenvectors);

    // Print the eigenvalues
    print_eigenvalues(25, eigenvalues);

    // User specifies which eigenvector and range of elements to print
    int vec_index   =  0;   // For example, the first eigenvector
    int start_index = 75;   // Start at index 75 (second-to-last element)
    int end_index   = 99;   // End at index 99 (last element)
    print_eigenvector_range(N, eigenvectors, vec_index, start_index, end_index);

    // Free allocated memory
    free(matrix->d);
    free(matrix->e);
    free(matrix);
    free(eigenvalues);
    free(eigenvectors);

    return 0;
}

// Function to calculate the potential mL2V
double mL2V(double y) {
    return 1000 * sin(20 * y) * pow(y, 4);
    // return 1000 * pow((y - 0.5), 2); // 2nd example
}

// Function to initialize the tridiagonal matrix
TridiagonalMatrix* init_tridiagonal_matrix(int N, double (*potential)(double)) {
    double dy = 1.0 / N;
    TridiagonalMatrix *matrix = malloc(sizeof(TridiagonalMatrix));
    matrix->size = N;
    matrix->d = malloc(N * sizeof(double));
    matrix->e = malloc((N - 1) * sizeof(double));

    for (int i = 0; i < N; i++) {
        double y = i * dy;
        matrix->d[i] = 1.0 / (dy * dy) + potential(y);
        if (i < N - 1) {
            matrix->e[i] = -1.0 / (2.0 * dy * dy);
        }
    }
    return matrix;
}

// Function to solve the eigenvalue problem
void solve_eigenproblem(TridiagonalMatrix *matrix, double *eigenvalues, double *eigenvectors) {
    LAPACKE_dstev(LAPACK_ROW_MAJOR, 'V', matrix->size, matrix->d, matrix->e, eigenvectors, matrix->size);
    for (int i = 0; i < matrix->size; i++) {
        eigenvalues[i] = matrix->d[i];
    }
}

// Function to print the matrix (for debugging purposes)
void print_matrix(const TridiagonalMatrix *matrix) {
    int max_print_size = 6;

    for (int i = 0; i < max_print_size; i++) {
        for (int j = 0; j <max_print_size; j++) {
            if (i == j) {
                printf("%12.5e ", matrix->d[i]);
            } else if (abs(i - j) == 1) {
                printf("%12.5e ", matrix->e[i < j ? i : j]);
            } else {
                printf("%12.5e ", 0.0);
            }
        }
        printf("\n");
    }
}

void print_eigenvalues(int N, double *eigenvalues) {
    printf("\nEigenvalues:\n");
    for (int i = 0; i < N; i++) {
        printf("%12.5e ", eigenvalues[i]);
        if ((i + 1) % 5 == 0) { // Print 5 elements per line for better readability
            printf("\n");
        }
    }
    printf("\n");
}

// Enhanced function to print a specific range of elements from an eigenvector
void print_eigenvector_range(int N, double *eigenvectors, int vec_index, int start_index, int end_index) {
    if (vec_index < 0 || vec_index >= N) {
        printf("Error: Eigenvector index out of bounds. Please choose a value between 0 and %d.\n", N-1);
        return;
    }
    if (start_index < 0 || start_index >= N || end_index < 0 || end_index >= N || start_index > end_index) {
        printf("Error: Invalid range. Please choose start_index and end_index between 0 and %d, with start_index <= end_index.\n", N-1);
        return;
    }

    printf("Eigenvector %d (elements %d to %d):\n", vec_index + 1, start_index + 1, end_index + 1);
    for (int i = start_index; i <= end_index; i++) {
        printf("%12.5e ", eigenvectors[i * N + vec_index]);
        if ((i - start_index + 1) % 5 == 0) { // Print 5 elements per line for better readability
            printf("\n");
        }
    }
    printf("\n");
}
