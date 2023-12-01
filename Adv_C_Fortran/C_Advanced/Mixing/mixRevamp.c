/*
===============================================================================
 This code needs the FFTW libray to run the code:
    $ sudo apt install libfftw3-dev
===============================================================================
 Complile and execute with:
    $ gcc mixRevamp.c -o mixR -lm -lfftw3
    $ ./mix2
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
 This is the Revamped version of the 'mixing.c' code. First of all, it is 
 finally working correctly. In the orginal version 'c_init' served no purpose
 and was basically a dummy argument adding a matrix of zeroes. Instead, what it
 was supposed to do was add the original sine wave back into the mix. Now, this
 function works and 'p_diffuse' together with 'iterate' give the correct values
 for the contour plot. Also, there were sevral bugs in the plotting code as well, 
 but these have now been fixed in the 'plot_mixC.py' code.

 Second important feature is the intoduction is the addition of structs. These
 help eliminate boilerplate and introduce matrices that are dynamically allocated.
 This also means that N is no longer capped at '512' and can go to giher powers of
 2, although though the plots don't necessarily require it. Nevertheless, the code 
 does not cause any memory leaks anymore.

 Unfortunately, there is no good function for fftshift, which is why that is 
 handled in the python code for plotting. This at least asures that the plots are
 correct. Although, there is the function fft2() which uses functions from the 
 fftw3 library and works exaclty like the Python one. 

 There are still ways to improve the code. Namely all matrices can be vectorized 
 with the MatrixVector sturcts, which would improve the performance of the code.
 I will prolly make another version of this in a new file.
===============================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>
#include <fftw3.h>
#define N 512
#define pi 3.14159265358979323846

typedef struct {
    double **data;
    int rows;
    int cols;
} Matrix;

// Define a struct for a complex matrix
typedef struct {
    int rows;
    int cols;
    double complex **data;
} ComplexMatrix;

Matrix  p_diffuse(Matrix *matrix, bool in_x_direction);
Matrix  iterate(int no_of_times, Matrix *matrix);
void    printMatrix(Matrix mat);
Matrix  matZeros(int rows, int cols);
void    matFree(Matrix mat);
double  rand_range(double min, double max);
Matrix  matRand(int rows, int cols);
void    matrix_csv(char *filename, Matrix *matrix);
ComplexMatrix   matZerosComplex(int rows, int cols);
void            printComplexMatrix(ComplexMatrix mat);
void            matFreeComplex(ComplexMatrix mat);
ComplexMatrix convertToComplex(Matrix *mat);
ComplexMatrix fft2(Matrix *input);
ComplexMatrix fftshift(ComplexMatrix *mat);

int main() {
    struct timespec start, stop;
    //Declaraion of initial variables
    double a = 200.0;
    Matrix c1 = matZeros(N,N);
    int i, j, k;

    for (i = 0; i < c1.rows; i++) { //same code that was in lecture: initial sine
        for (j = 0; j < c1.cols; j++) {
            c1.data[i][j] = sin(j * pi * 2 / (double) N); //sine on i-axis
        }
    }

    //printMatrix(c1);

    // Get starting time
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    int numCycles;
    printf("Input the number of mixing cycles:\n");
    scanf("%d", &numCycles);

    c1 = iterate(numCycles, &c1);

    // Get end time
    clock_gettime(CLOCK_MONOTONIC, &stop);

    // Print matrix after 11 diffusion iterations
    printf("c1 = ");
    printMatrix(c1);

    // Calculate the elapsed time in seconds
    double time_taken = (stop.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (stop.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("\nDiffusion iterations took %.3e seconds.\n", time_taken);

    // Write contour plot data into a CSV file
    char str[14] = "colsRevamp";
    matrix_csv(str, &c1);

    /////////////////////////////////////////////////////////////////////

    // Initialize a complex matrix with zeros
    //ComplexMatrix mat_fft2 = matZerosComplex(N, N);

    ComplexMatrix mat_fft2 = fft2(&c1);

    printf("mat_fft2 = ");
    printComplexMatrix(mat_fft2);

    ComplexMatrix mat_fftshift = fftshift(&mat_fft2);
    
    /* Can't get fftshift to work :( */
    printf("mat_fftshift = ");
    printComplexMatrix(mat_fftshift);
    /*
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            c1.data[i][i] = sqrt(creal(mat_fftshift.data[i][j]*mat_fftshift.data[i][j]) + 
                                 cimag(mat_fftshift.data[i][j]*mat_fftshift.data[i][j]) );
        }
    }

    printMatrix(c1);

    // Write
    char str2[14] = "fourRevamp";
    matrix_csv(str2, &c1);
    */

    //Free allocated memory
    matFree(c1);
    matFreeComplex(mat_fft2);
    matFreeComplex(mat_fftshift);
    
    return 0;
}
//================================================================================
// Start of the diffusion function
//================================================================================
Matrix p_diffuse(Matrix *matrix, bool in_x_direction) {
    // Initialize new matrices with zeros
    Matrix new_matrix = matZeros(matrix->rows, matrix->cols);
    Matrix c_init     = matZeros(matrix->rows, matrix->cols);
    
    // Set the diffusion coefficient
    double a = 200.0;

    // Loop through each element in the matrix
    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < matrix->cols; j++) {
            double w1, w2;
            c_init.data[i][j] = sin(j * pi * 2 / (double) N); //sine on i-axis
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

    // Return the resulting matrix
    return new_matrix;
    
    // Note: The following lines will not be executed as they are after the return statement.
    // Still, I run them just in case as a sanity check
    matFree(new_matrix);
    matFree(c_init);
}

//================================================================================
// Start of the iterate function
//================================================================================
Matrix iterate(int no_of_times, Matrix *matrix) {
    int k, i, j;
    double a = 200;
    Matrix c = matZeros(matrix->rows, matrix->rows);

    for (k = 0; k < no_of_times; k++) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                if (k % 2 == 0) { //even => mixes one way, odd => another way
                    c.data[i][j] = matrix->data[i][(j + (int)(a * cos(i * pi * 2 / N))) & N - 1];
                } else {
                    c.data[i][j] = matrix->data[(i + (int)(a * cos(j * pi * 2 / N))) & N - 1][j];
                }
            }
        }
        c = p_diffuse(&c, k % 2 == 0);
        for (i = 0; i < N; i++) { //put elements of new matrix into the old one, continue changing matrix c.
            for (j = 0; j < N; j++) {
                matrix->data[i][j] = c.data[i][j];
            }
        }
    }
    for (i = 0; i < N; i++) { //put elements of new matrix into the old one, continue changing matrix c.
        for (j = 0; j < N; j++) {
            c.data[i][j] = matrix->data[i][j];
        }
    }
    return c;
    matFree(c);
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Function to print a matrix
void printMatrix(Matrix mat) {
    int max_print_size = 5;

    printf("[\n");

    if (mat.rows <= max_print_size && mat.cols <= max_print_size) {
        for (int i = 0; i < mat.rows; i++) {
            printf("  [");
            for (int j = 0; j < mat.cols; j++) {
                printf("%f", mat.data[i][j]);
                if (j < mat.cols - 1) printf(", ");
            }
            printf("]");
            if (i < mat.rows - 1) printf(",\n");
        }
    } else {
        for (int i = 0; i < max_print_size; i++) {
            printf("  [");
            for (int j = 0; j < max_print_size; j++) {
                printf("%f", mat.data[i][j]);
                if (j < max_print_size - 1) printf(", ");
            }
            printf(", ...");
            printf("]");
            if (i < max_print_size - 1) printf(",\n");
        }
        printf(",\n\t\t\t  ...\n");
        //printf("  ...");
        //printf("  ...");
    }

    printf("]\n");
}

Matrix matZeros(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;

    mat.data = (double **)malloc(rows * sizeof(double *));
    if (mat.data == NULL) {
        fprintf(stderr, "Memory allocation failed for matrix rows.\n");
        mat.rows = 0; // Reset rows in case of failure
        return mat;
    }

    for (int i = 0; i < rows; i++) {
        mat.data[i] = (double *)calloc(cols, sizeof(double));
        if (mat.data[i] == NULL) {
            fprintf(stderr, "Memory allocation failed for matrix columns.\n");

            // Clean up memory allocated so far
            for (int j = 0; j < i; j++) {
                free(mat.data[j]);
            }
            free(mat.data);

            mat.rows = 0; // Reset rows in case of failure
            return mat;
        }
    }

    return mat;
}

// Free the memory allocated for the matrix.
void matFree(Matrix mat) {
    for (int i = 0; i < mat.rows; i++) {
        free(mat.data[i]);    // Free memory for each row.
    }
    free(mat.data);           // Free memory for the array of pointers
}

// Generate a random floating-point number within the given range.
double rand_range(double min, double max) {
    return min + (max-min) * ((double)rand()/RAND_MAX);
}

Matrix matRand(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;

    // Allocate memory for an array of pointers to double arrays (rows).
    mat.data = (double **)malloc(rows * sizeof(double *));

    if (mat.data == NULL) {
        perror("Memory allocation failed");
        mat.rows = 0; // Reset rows in case of failure
        return mat;
    }

    // Allocate memory for each row of the matrix.
    for (int i = 0; i < rows; i++) {
        mat.data[i] = (double *)malloc(cols * sizeof(double));

        if (mat.data[i] == NULL) {
            perror("Memory allocation failed");

            // Clean up memory allocated so far
            for (int j = 0; j < i; j++) {
                free(mat.data[j]);
            }
            free(mat.data);

            mat.rows = 0; // Reset rows in case of failure
            return mat;
        }
    }

    // Seed the random number generator based on the current time.
    srand((unsigned int)time(NULL));

    // Fill the matrix with random numbers.
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            mat.data[i][j] = rand_range(0.0, 1.0);
        }
    }

    return mat;
}

void matrix_csv(char *filename, Matrix *matrix) {
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
            fprintf(fp, "%10.14e, ", matrix->data[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("\n%s file created.\n", filename);
}

///////////////////////////////////////////////////////////////////////////

// Function to initialize a complex matrix with zeros
ComplexMatrix matZerosComplex(int rows, int cols) {
    ComplexMatrix mat;
    mat.rows = rows;
    mat.cols = cols;

    mat.data = (double complex **)malloc(rows * sizeof(double complex *));
    if (mat.data == NULL) {
        fprintf(stderr, "Memory allocation failed for matrix rows.\n");
        mat.rows = 0; // Reset rows in case of failure
        return mat;
    }

    for (int i = 0; i < rows; i++) {
        mat.data[i] = (double complex *)calloc(cols, sizeof(double complex));
        if (mat.data[i] == NULL) {
            fprintf(stderr, "Memory allocation failed for matrix columns.\n");

            // Clean up memory allocated so far
            for (int j = 0; j < i; j++) {
                free(mat.data[j]);
            }
            free(mat.data);

            mat.rows = 0; // Reset rows in case of failure
            return mat;
        }
    }

    return mat;
}

// Function to print a complex matrix
void printComplexMatrix(ComplexMatrix mat) {
    int max_print_size = 5;

    printf("[\n");

    if (mat.rows <= max_print_size && mat.cols <= max_print_size) {
        for (int i = 0; i < mat.rows; i++) {
            printf("  [");
            for (int j = 0; j < mat.cols; j++) {
                printf("(%.3e + %.3ei)", creal(mat.data[i][j]), cimag(mat.data[i][j]));
                if (j < mat.cols - 1) printf(", ");
            }
            printf("]");
            if (i < mat.rows - 1) printf(",\n");
        }
    } else {
        for (int i = 0; i < max_print_size; i++) {
            printf("  [");
            for (int j = 0; j < max_print_size; j++) {
                printf("(%.3e + %.3ei)", creal(mat.data[i][j]), cimag(mat.data[i][j]));
                if (j < max_print_size - 1) printf(", ");
            }
            printf(", ...");
            printf("]");
            if (i < max_print_size - 1) printf(",\n");
        }
        printf(",\n\t\t\t  ...\n");
        //printf("  ...");
    }

    printf("]\n");
}

// Function to free the memory allocated for a complex matrix
void matFreeComplex(ComplexMatrix mat) {
    for (int i = 0; i < mat.rows; i++) {
        free(mat.data[i]);    // Free memory for each row.
    }
    free(mat.data);           // Free memory for the array of pointers
}

// Function to convert a matrix to a matrix of complex numbers
ComplexMatrix convertToComplex(Matrix *mat) {
    ComplexMatrix matComplex = matZerosComplex(mat->rows, mat->cols);

    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            matComplex.data[i][j] = mat->data[i][j];
        }
    }

    return matComplex;
}

// Function to perform 2D FFT on a Matrix and return a ComplexMatrix
ComplexMatrix fft2(Matrix *input) {
    ComplexMatrix result = matZerosComplex(input->rows, input->cols);

    // Allocate memory for FFTW input and output arrays
    double *in_data = (double *)malloc(input->rows * input->cols * sizeof(double));
    fftw_complex *out_data = (fftw_complex *)malloc(input->rows * (input->cols / 2 + 1) * sizeof(fftw_complex));

    // Flatten the input matrix into a 1D array
    for (int i = 0; i < input->rows; ++i) {
        for (int j = 0; j < input->cols; ++j) {
            in_data[i * input->cols + j] = input->data[i][j];
        }
    }

    // Create FFTW plan
    fftw_plan plan = fftw_plan_dft_r2c_2d(input->rows, input->cols, in_data, out_data, FFTW_ESTIMATE);

    // Execute the plan (perform 2D FFT)
    fftw_execute(plan);    

    // Copy the complex data from FFTW output to the result
    for (int i = 0; i < result.rows; ++i) {
        for (int j = 0; j < result.cols; ++j) {
            result.data[i][j] = out_data[i * result.cols + j];
        }
    }

    // Clean up
    fftw_destroy_plan(plan);
    free(in_data);
    free(out_data);

    return result;
}

ComplexMatrix fftshift(ComplexMatrix *mat) {
    
    int rows = mat->rows;
    int cols = mat->cols;
    complex double tmp;
    ComplexMatrix matFinal = matZerosComplex(rows, cols);

    /* Swap elements in the first and second halves of the rows */
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols / 2; j++) {
            tmp = mat->data[i][j];
            mat->data[i][j] = mat->data[i][j + cols / 2];
            mat->data[i][j + cols / 2] = tmp;
        }
    }

    /* Swap elements in the first and second halves of the columns */
    for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows / 2; i++) {
            tmp = mat->data[i][j];
            mat->data[i][j] = mat->data[i + rows / 2][j];
            mat->data[i + rows / 2][j] = tmp;
        }
    }

    for (int j = 0; j < rows; j++) {
        for (int i = 0; i < cols; i++) {
            matFinal.data[i][j] = mat->data[i][j];
        }
    }

    return matFinal;
    matFreeComplex(matFinal);
}
