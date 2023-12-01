/*
===============================================================================
 This code needs the FFTW libray to run the code:
    $ sudo apt install libfftw3-dev
===============================================================================
 Complile and execute with:
    $ gcc mixFast.c -o mixF -lm 
    $ ./mixF
===============================================================================
 How to set up and compile on a Mac:
    $ brew install fftw
    $ brew --prefix fftw
    > /opt/homebrew/opt/fftw # The path where fftw is installed using homebrew
    $ gcc -o mix mixing.c -L/opt/homebrew/opt/fftw/lib -lfftw3
    $ ./mix
 ---
===============================================================================
Documentation:
 The fastest possible version to perform turbulent mixing using bitwise 
 manipulation.
===============================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
//#include <complex.h>
//#include <fftw3.h>

#define N 512
#define pi 3.14159265358979323846

typedef struct {
    double *data;
    int rows;
    int cols;
} MatrixVector;

MatrixVector    matvecZeros(int rows, int cols);
void            matvecFree(MatrixVector matvec);
MatrixVector    matMulVector(const MatrixVector *A, const MatrixVector *B);
void            printMatrixVector(MatrixVector matvec);
MatrixVector    p_diffuse(MatrixVector *matrix, bool in_x_direction);
MatrixVector    iterate(int no_of_times, MatrixVector *matrix);
void            matrix2csv(char *filename, MatrixVector *matrix);

int main() {
    struct timespec start, stop;
    //Declaraion of initial variables
    double a = 200.0;
    MatrixVector c = matvecZeros(N, N);
    int i, j, k;

    for (i = 0; i < c.rows; i++) { //same code that was in lecture: initial sine
        for (j = 0; j < c.cols; j++) {
            c.data[i*c.rows + j] = sin(j * pi * 2 / (double) N); //sine on i-axis
        }
    }
    
    // Ask the user an input for mixing cycles
    int numCycles;
    printf("Input the number of mixing cycles:\n");
    scanf("%d", &numCycles);

    // Get starting time
    clock_gettime(CLOCK_MONOTONIC, &start);

    // Iterate the diffusion process 'numCycles' of times
    c = iterate(numCycles, &c);
    
    // Get end time
    clock_gettime(CLOCK_MONOTONIC, &stop);

    // Print matrix
    printf("c = ");
    printMatrixVector(c);

    // Calculate the elapsed time in seconds
    double time_taken = (stop.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (stop.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("\nDiffusion iterations took %.3e seconds.\n", time_taken);

    // Write contour plot data into a CSV file
    char str[12] = "coloursF";
    matrix2csv(str, &c);

    // Free allocated memory
    matvecFree(c);
    return 0;
}

// Function to create a new MatrixVector and initialize with zeros
MatrixVector matvecZeros(int rows, int cols) {
    MatrixVector matvec;
    matvec.rows = rows;
    matvec.cols = cols;
    matvec.data = (double *)calloc(rows * cols, sizeof(double));

    if (matvec.data == NULL) {
        fprintf(stderr, "Memory allocation failed for MatrixVector.\n");
        matvec.rows = 0; // Reset rows in case of failure
    }

    return matvec;
}

// Function to free the memory allocated for MatrixVector
void matvecFree(MatrixVector matvec) {
    free(matvec.data);
}

// Function to print a MatrixVector
void printMatrixVector(MatrixVector matvec) {
    int max_print_size = 5;

    printf("[\n");

    if (matvec.rows <= max_print_size) {
        for (int i = 0; i < matvec.rows; i++) {
            printf("  [");
            for (int j = 0; j < matvec.cols; j++) {
                printf("%f", matvec.data[i * matvec.cols + j]);
                if (j < matvec.cols - 1) printf(", ");
            }
            printf("]");
            if (i < matvec.rows - 1) printf(",\n");
        }
    } else {
        for (int i = 0; i < max_print_size; i++) {
            printf("  [");
            for (int j = 0; j < max_print_size; j++) {
                printf("%f", matvec.data[i * matvec.cols + j]);
                if (j < matvec.cols - 1) printf(", ");
            }
            printf(" ...");
            printf("]");
            if (i < max_print_size - 1) printf(",\n");
        }
        printf(",\n\t\t\t  ...\n");
    }
    printf("]\n");
}

MatrixVector p_diffuse(MatrixVector *matrix, bool in_x_direction) {
    /*=============================================================
    Documentation:
    ---------------------------------------------------------------
    The p_diffuse function is designed to perform a diffusion 
    operation on a given matrix in either the x-direction or the 
    y-direction, based on the value of the boolean parameter 
    in_x_direction. This diffusion process is a crucial step in 
    simulating turbulent mixing. 
    =============================================================*/
    // Initialize new matrices with zeros
    MatrixVector c_new  = matvecZeros(matrix->rows, matrix->cols);
    MatrixVector c_init = matvecZeros(matrix->rows, matrix->cols);
    
    int rows = matrix->rows; 
    int cols = matrix->cols; 

    // Set the diffusion coefficient
    double a = 200.0;

    // Loop through each element in the matrix
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double w1, w2;
            c_init.data[i*rows + j] = sin(j * pi * 2 / (double) N); //sine on i-axis
            // Calculate weights based on the direction of diffusion
            if (in_x_direction) {
                // If diffusion is in the x-direction
                w1 = a * cos(j * 2.0 * pi / cols) - floor(a * cos(j * 2.0 * pi / cols));
                w2 = 1 - w1;

                // Handle the special case when j is 0
                if (j == 0) {
                    // Wrap around to the last column for the diffusion calculation
                    c_new.data[i*rows + j] = matrix->data[i*rows + (cols-1)] * w1 + matrix->data[i*rows + j] * w2 + c_init.data[i*rows + j];
                } else {
                    // Standard diffusion calculation
                    c_new.data[i*rows + j] = matrix->data[i*rows + (j-1)] * w1 + matrix->data[i*rows + j] * w2 + c_init.data[i*rows + j];
                }
            } else {
                // If diffusion is in the y-direction
                w1 = a * cos(i * 2.0 * pi / rows) - floor(a * cos(i * 2.0 * pi / rows));
                w2 = 1 - w1;

                // Handle the special case when i is 0
                if (i == 0) {
                    // Wrap around to the last row for the diffusion calculation
                    c_new.data[i*rows + j] = matrix->data[(rows-1)*rows + j] * w1 + matrix->data[i*rows + j] * w2 + c_init.data[i*rows + j];
                } else {
                    // Standard diffusion calculation
                    c_new.data[i*rows + j] = matrix->data[(i-1)*rows + j] * w1 + matrix->data[i*rows + j] * w2 + c_init.data[i*rows + j];
                }
            }
        }
    }
    matvecFree(c_init);

    // Return the resulting matrix
    return c_new;
    
    //matvecFree(c_new);
}

MatrixVector iterate(int no_of_times, MatrixVector *matrix) {
    /*=======================================================
    Documentation:
    ---------------------------------------------------------
    The iterate function is responsible for performing a 
    series of mixing and diffusion operations on a given 
    matrix over a specified number of iterations. This 
    function plays a key role in simulating turbulent mixing.
    =======================================================*/
    int k, i, j;
    double a = 200;
    MatrixVector c = matvecZeros(matrix->rows, matrix->rows);

    int rows = matrix->rows; 
    int cols = matrix->cols; 

    for (k = 0; k < no_of_times; k++) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                if (k % 2 == 0) { //even => mixes one way, odd => another way
                    c.data[i*rows + j] = matrix->data[i*rows + ((j + (int)(a * cos(i * pi * 2 / N))) & N - 1)];
                } else {
                    c.data[i*rows + j] = matrix->data[((i + (int)(a * cos(j * pi * 2 / N))) & N - 1)*rows + j];
                }
            }
        }
        c = p_diffuse(&c, k % 2 == 0);
        for (i = 0; i < N; i++) { //put elements of new matrix into the old one, continue changing matrix c.
            for (j = 0; j < N; j++) {
                matrix->data[i*rows + j] = c.data[i*rows + j];
            }
        }
    }
    for (i = 0; i < N; i++) { //put elements of new matrix into the old one, continue changing matrix c.
        for (j = 0; j < N; j++) {
            c.data[i*rows + j] = matrix->data[i*rows + j];
        }
    }
    return c;
    //matvecFree(c);
}

void matrix2csv(char *filename, MatrixVector *matrix) {
    printf("\nCreating %s.csv file for matrix.", filename);
    FILE *fp;
    int i, j;
    filename = strcat(filename, ".csv");
    fp = fopen(filename, "w+");

    if (fp == NULL) {
        fprintf(stderr, "Error opening file %s for writing.\n", filename);
        return;
    }

    for (i = 0; i < matrix->rows; i++) {
        for (j = 0; j < matrix->cols; j++) {
            fprintf(fp, "%8.6e, ", matrix->data[i*matrix->rows + j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("\n%s file created.\n", filename);
}