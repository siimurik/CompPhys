/*
===============================================================================
 Complile and execute with:
    $ gcc mixing.c -o mix -lm -lfftw3
    $ ./mix
===============================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>
#include <fftw3.h>
#define N 512
#define M_PI 3.14159265358979323846

//void p_diffuse(double c_init[N][N], double matrix[N][N], bool in_x_direction);
//void iterate(int no_of_times, double c_old[N][N], double c_init[N][N]);
void p_diffuse(double matrix[N][N], bool in_x_direction);
void iterate(int no_of_times, double c_old[N][N]);
void matrix_csv(char *filename,double a[][N],int n,int m);
void write_array_to_csv(const char *filename, const complex double data[][N], int rows, int cols);
void fft2(complex double (*data)[N], int rows, int cols);
void fft1(complex double *data, int n);
void fftshift(complex double data[][N], int rows, int cols);

int main() {
    double a = 200;

    static double c1[N][N] = {0};
    static double c_final[N][N] = {0};
    static double c_init[N][N] = {0};
    int i, j;


    for (i = 0; i < N; i++) { //same code that was in lecture: initial sine
        for (j = 0; j < N; j++) {
            c1[i][j] = sin(j * M_PI * 2 / N); //sine on i-axis
            //c_init[i][j] = sin(j * M_PI * 2 / N); //same sine for the iterate() command, stays constant.
        }
    }

    //signum(int(N / 2), c1);
    //diffuse(amount, c1, True);

    //iterate(11, c1, c_init);
    iterate(11, c1);
    // print the result (only the first few elements)
    for (i = 0; i < 5; i++){
        for (j = 0; j < 5; j++){
            printf("%f ", c1[i][j]);
        }
        printf("\n");
    }
    char str[10] = "colors";
    matrix_csv(str,c1,N,N);

    // Convert the c1 array to a fftw_complex array
    // static fftw_complex *input;
    // init_input(input, c1);

    // Allocate memory for the input array
    fftw_complex *input = fftw_malloc(sizeof(fftw_complex) * N * N);
    if (input == NULL) {
        fprintf(stderr, "Error allocating memory for input array\n");
        return 1;
    }

    // Initialize the elements of the input array
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            input[i * N + j] = c1[i][j];
        }
    }

    // Allocate memory for the FFT of the input array
    //fftw_complex *fft = fftw_malloc(sizeof(fftw_complex) * N * N);
    //if (fft == NULL) {
    //fprintf(stderr, "Error allocating memory for FFT\n");
    //return 1;
    //}

    // Allocate memory for the FFT of the input array
    static complex double fft[N][N];

    // Create a plan for computing the FFT of the input array
    fftw_plan plan = fftw_plan_dft_2d(N, N, input, &fft[0][0], FFTW_FORWARD, FFTW_ESTIMATE);

    // Compute the FFT of the input array
    fftw_execute(plan);

    // Shift the FFT output so that the origin is at the center of the array
    fftshift(fft, N, N);

    // Print the first few elements of the FFT output
    for (int i = 0; i < 5; i++){
        for (int j = 0; j < 5; j++){
            printf("%f ", fft[i][j]);
        }
    printf("\n");
    }

    // Write the values in the data array to a CSV file
    write_array_to_csv("fourier.csv", fft, N, N);

    // Clean up
    fftw_destroy_plan(plan);
    fftw_free(input);
    //fftw_free(fft);

    return 0;
}
//---------------------------------------------------------------------------------
void p_diffuse(double matrix[N][N], bool in_x_direction) {
    double new_matrix[N][N] = {0};
    //double c_init[N][N] = {0};
    double c_init[N][N];
    int i, j;
    double a = 200;

    for (i = 0; i < N; i++) { //same code that was in lecture: initial sine
        for (j = 0; j < N; j++) {
            c_init[i][j] = sin(j * M_PI * 2 / N); //same sine for the iterate() command, stays constant.
        }
    }
    for (i = 0; i < N; i++) { 
        for (j = 0; j < N; j++) {
            if (in_x_direction) { //if transformation was made horizontally
                double w1 = a * cos(j * M_PI * 2 / N) - floor(a * cos(j * M_PI * 2 / N));
                double w2 = 1 - w1;
                new_matrix[i][j] = w1 * matrix[i][j - 1] + w2 * matrix[i][j] + c_init[i][j]; //here we add original sine again
            } else { //else => transformation was made vertically
                double w1 = a * cos(i * M_PI * 2 / N) - floor(a * cos(i * M_PI * 2 / N));
                double w2 = 1 - w1;
                new_matrix[i][j] = w1 * matrix[i - 1][j] + w2 * matrix[i][j] + c_init[i][j]; //here we add original sine again
            }
        }
    }
}
//---------------------------------------------------------------------------------
void iterate(int no_of_times, double c_old[N][N]) {
    double c[N][N] = {0};
    int k, i, j;
    double a = 200;
    //double c_init[N][N];

    //for (i = 0; i < N; i++) { //same code that was in lecture: initial sine
    //    for (j = 0; j < N; j++) {
    //        c_init[i][j] = sin(j * M_PI * 2 / N); //same sine for the iterate() command, stays constant.
    //    }
    //}

    for (k = 0; k < no_of_times; k++) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                if (k % 2 == 0) { //even => mixes one way, odd => another way
                    c[i][j] = c_old[i][(j + (int)(a * cos(i * M_PI * 2 / N))) & N - 1];
                } else {
                    c[i][j] = c_old[(i + (int)(a * cos(j * M_PI * 2 / N))) & N - 1][j];
                }
            }
        }
        //diffuse(amount, c, k % 2 == 0);
        p_diffuse(c, k % 2 == 0);
        //p_diffuse(c_init, c, k % 2 == 0);

        for (i = 0; i < N; i++) { //put elements of new matrix into the old one, continue changing matrix c.
            for (j = 0; j < N; j++) {
                c_old[i][j] = c[i][j];
            }
        }
    }
    for (i = 0; i < N; i++) { //put elements of new matrix into the old one, continue changing matrix c.
        for (j = 0; j < N; j++) {
            c[i][j] = c_old[i][j];
        }
    }
}
//---------------------------------------------------------------------------------
// Function that creates and stores the temperature matrix values
void matrix_csv(char *filename,double a[][N],int n,int m){
    printf("\n Creating %s.csv file for matrix.",filename);
    FILE *fp;
    int i,j;
    filename=strcat(filename,".csv");
    fp=fopen(filename,"w+");
    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++){
            fprintf(fp,"%20.14f, ",a[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    printf("\n %s file created.\n",filename);
}

/*
================================================================================
Documentation for fftshift():
 This function takes a two-dimensional array of double-precision floating-point 
 values, along with the number of rows and columns in the array, and shifts the 
 zero-frequency component to the center of the spectrum. It does this by swapping 
 the values in the first and second half of each row and column, which has the 
 effect of shifting the zero-frequency component to the center of the array.

 This function shifts the elements of a two-dimensional array by swapping 
 elements in the first and second halves of the rows and columns. This is 
 equivalent to the np.fft.fftshift() function in NumPy.

 To use this function, you would call it with the two-dimensional array data and
 the number of rows and columns in the array as arguments, e.g., 
 fftshift((complex double **) data, rows, cols). This will shift the elements of
 the array in-place (i.e., the original array will be overwritten with the 
 shifted array).
================================================================================
*/
void fftshift(complex double data[][N], int rows, int cols)
{
    complex double tmp;

    /* Swap elements in the first and second halves of the rows */
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols / 2; j++) {
            tmp = data[i][j];
            data[i][j] = data[i][j + cols / 2];
            data[i][j + cols / 2] = tmp;
        }
    }

    /* Swap elements in the first and second halves of the columns */
    for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows / 2; i++) {
            tmp = data[i][j];
            data[i][j] = data[i + rows / 2][j];
            data[i + rows / 2][j] = tmp;
        }
    }
}

/*
================================================================================
Documentation for fft2():
 This function takes a two-dimensional array of double-precision floating-point 
 values, along with the number of rows and columns in the array, and computes 
 the two-dimensional discrete Fourier transform (DFT) of the input data. It does
 this by using the FFTW library to perform the actual FFT computation, and then
 scaling and computing the magnitude of the complex output values to obtain the 
 final result.
================================================================================
*/
void fft2(complex double (*data)[N], int rows, int cols)
{
    complex double tmp[N];

    /* Allocate memory for temporary array */
    //tmp = (complex double **) malloc(rows * cols * sizeof(complex double));

    /* Compute 1D FFT for each row */
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            tmp[j] = data[i][j];
        }
        fft1(tmp, cols);
        for (int j = 0; j < cols; j++) {
            data[i][j] = tmp[j];
        }
    }

    /* Compute 1D FFT for each column */
    for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            tmp[i] = data[i][j];
        }
        fft1(tmp, rows);
        for (int i = 0; i < rows; i++) {
            data[i][j] = tmp[i];
        }
    }

    /* Free memory for temporary array */
    //free(tmp);
}

/*
This function computes the one-dimensional discrete Fourier transform 
(DFT) of the input array data using the Cooley-Tukey algorithm. It 
takes in the array data and its size n as arguments, and it returns 
the DFT of the array in-place (i.e., the original array is overwritten 
with the DFT of the array).

The Cooley-Tukey algorithm is a fast Fourier transform (FFT) algorithm 
that divides the input array into even and odd indices, computes the DFT 
of each subarray, and then combines the results to obtain the DFT of the 
entire array. This allows the computation of the DFT to be performed in 
O(n log n) time, which is much faster than the O(n^2) time complexity 
of the naive DFT algorithm.
*/
void fft1(complex double *data, int n)
{
    complex double *tmp;
    complex double omega, omega_k, omega_n;

    /* Allocate memory for temporary array */
    tmp = (complex double *) malloc(n * sizeof(complex double));

    /* Compute FFT using Cooley-Tukey algorithm */
    for (int k = 0; k < n; k++) {
        tmp[k] = 0;
        omega_n = cexp(-2 * M_PI * I / n);
        omega_k = pow(omega_n, k);
        for (int t = 0; t < n; t++) {
            omega = pow(omega_n, t * k);
            tmp[k] += data[t] * omega;
        }
    }

    /* Copy results from temporary array to input array */
    for (int i = 0; i < n; i++) {
        data[i] = tmp[i];
    }

    /* Free memory for temporary array */
    free(tmp);
}

// Function to write the values in a 2D array to a CSV file
void write_array_to_csv(const char *filename, const complex double data[][N], int rows, int cols) {
    printf("\n Creating %s.csv file for matrix.",filename);
    // Open a file for writing the CSV data
    FILE *fp = fopen(filename, "w");
    // Check if the file was opened successfully
    if (fp == NULL) {
        fprintf(stderr, "Error opening file\n");
        return;
    }

    // Write the values in the array to the file
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            fprintf(fp, "%f", data[i][j]);
            if (j < cols - 1) {
                fprintf(fp, ",");
            }
        }
        fprintf(fp, "\n");
    }
    
    // Close the file
    fclose(fp);
    printf("\n %s file created.\n",filename);
}