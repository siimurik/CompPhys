#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

// Define the matrix structure
typedef struct {
    int rows;
    int cols;
    double **data;
} Matrix;

// Function to create a zero matrix
Matrix matZeros(int rows, int cols) {
    Matrix matrix;
    matrix.rows = rows;
    matrix.cols = cols;
    matrix.data = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++) {
        matrix.data[i] = (double *)calloc(cols, sizeof(double));
    }
    return matrix;
}

// Function to free a matrix
void matFree(Matrix matrix) {
    for (int i = 0; i < matrix.rows; i++) {
        free(matrix.data[i]);
    }
    free(matrix.data);
}

// The p_diffuse function
Matrix p_diffuse(Matrix *matrix, bool in_x_direction) {
    // Initialize new matrices with zeros
    Matrix new_matrix = matZeros(matrix->rows, matrix->cols);
    Matrix c_init = matZeros(matrix->rows, matrix->cols);

    // Set the diffusion coefficient
    double a = 200.0;
    const double pi = 3.141592653589793;

    // Loop through each element in the matrix
    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < matrix->cols; j++) {
            double w1, w2;
            c_init.data[i][j] = sin(j * pi * 2 / (double)matrix->cols); // sine on i-axis
            // Calculate weights based on the direction of diffusion
            if (in_x_direction) {
                // If diffusion is in the x-direction
                w1 = a * cos(j * 2.0 * pi / matrix->cols) - floor(a * cos(j * 2.0 * pi / matrix->cols));
                w2 = 1 - w1;

                // Handle the special case when j is 0
                if (j == 0) {
                    // Wrap around to the last column for the diffusion calculation
                    new_matrix.data[i][j] = matrix->data[i][matrix->cols - 1] * w1 + matrix->data[i][j] * w2 + c_init.data[i][j];
                } else {
                    // Standard diffusion calculation
                    new_matrix.data[i][j] = matrix->data[i][j - 1] * w1 + matrix->data[i][j] * w2 + c_init.data[i][j];
                }
            } else {
                // If diffusion is in the y-direction
                w1 = a * cos(i * 2.0 * pi / matrix->rows) - floor(a * cos(i * 2.0 * pi / matrix->rows));
                w2 = 1 - w1;

                // Handle the special case when i is 0
                if (i == 0) {
                    // Wrap around to the last row for the diffusion calculation
                    new_matrix.data[i][j] = matrix->data[matrix->rows - 1][j] * w1 + matrix->data[i][j] * w2 + c_init.data[i][j];
                } else {
                    // Standard diffusion calculation
                    new_matrix.data[i][j] = matrix->data[i - 1][j] * w1 + matrix->data[i][j] * w2 + c_init.data[i][j];
                }
            }
        }
    }

    // Free c_init matrix
    matFree(c_init);

    // Return the resulting matrix
    return new_matrix;
}

int main() {
    int rows = 5, cols = 5;
    Matrix matrix = matZeros(rows, cols);

    // Initialize the input matrix with some values
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix.data[i][j] = i + j;
        }
    }

    // Print the resulting matrix
    printf("Inital matrix:\n");
    for (int i = 0; i < matrix.rows; i++) {
        for (int j = 0; j < matrix.cols; j++) {
            printf("%f ", matrix.data[i][j]);
        }
        printf("\n");
    }

    int k, i, j;
    double a = 200;
    Matrix c = matZeros(rows, cols);
    int N = 5;
    const double pi = 3.141592653589793;
    for (k = 0; k < 2; k++) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                if (k % 2 == 0) {
                    c.data[i][j] = matrix.data[i][(j + (int)(a * cos(i * pi * 2 / N))) & N - 1];
                    printf("%d\n", (j + (int)(a * cos(i * pi * 2 / N))) & N - 1);
                } else {
                    c.data[i][j] = matrix.data[(i + (int)(a * cos(j * pi * 2 / N))) & N - 1][j];
                }
            }
        }
        c = p_diffuse(&c, k % 2 == 0);
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                matrix.data[i][j] = c.data[i][j];
            }
        }
    }

    // Print the resulting matrix
    printf("Final matrix:\n");
    for (int i = 0; i < matrix.rows; i++) {
        for (int j = 0; j < matrix.cols; j++) {
            printf("%f ", matrix.data[i][j]);
        }
        printf("\n");
    }

    // Call the diffusion function
    //Matrix result = p_diffuse(&matrix, true);

    // Free matrices
    matFree(matrix);
    //matFree(result);

    return 0;
}
