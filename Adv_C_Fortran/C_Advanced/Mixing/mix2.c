/* Compile and execute with
    $ gcc mix2.c -o mix2 -lm -lfftw3
    $ ./mix2
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

void p_diffuse(double matrix[N][N], bool in_x_direction);
void iterate(int no_of_times, double c_old[N][N]);
void matrix_csv(char *filename,double a[][N],int n,int m);
void fft2(complex double (*data)[N], int rows, int cols);
void fft1(complex double *data, int n);
void fftshift(complex double data[][N], int rows, int cols);

int main() {
    double a = 200;

    static double c1[N][N];
    static double c_final[N][N];
    static double c_init[N][N];
    static double cnew[N][N];
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            c1[i][j] = sin(j * M_PI * 2 / N);
        }
    }

    iterate(11, c1);

    for (int i = 0; i < 5; i++){
        for (int j = 0; j < 5; j++){
            printf("%f ", c1[i][j]);
        }
        printf("\n");
    }
    char str[10] = "colors";
    matrix_csv(str,c1,N,N);

    // Fourier transform
    complex double (*cdata)[N] = (complex double (*)[N]) c1;

    fft2(cdata, N, N);

    fftshift(cdata, N, N);

    printf("\nShifted FFT array:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%f + %fi \t", creal(cdata[i][j]), cimag(cdata[i][j]));
        }
        printf("\n");
    }
    printf("\n");

    // Loop through the rows and columns of the output of fftshift() and compute the absolute value of each complex number
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cnew[i][j] = log(cabs(cdata[i][j])*cabs(cdata[i][j]));
        }
    }

    char nstr[20] = "fourier";
    matrix_csv(nstr,cnew,N,N);

    return 0;
}

/*
 This function applies a diffusion transformation to a matrix. It takes a 
 two-dimensional array of double-precision floating-point values, a 
 boolean value indicating whether the transformation should be applied 
 in the x-direction or not, and an integer specifying the dimensions of 
 the array. The function then calculates the diffusion coefficients for 
 each element in the matrix, applies the transformation, and returns the 
 transformed matrix.
*/
void p_diffuse(double matrix[N][N], bool in_x_direction) {
//void p_diffuse(double c_init[N][N], double matrix[N][N], bool in_x_direction) {
    double new_matrix[N][N] = {0};
    //double c_init[N][N] = {0};
    double c_init[N][N] = {0};
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

void iterate(int no_of_times, double c_old[N][N]) {
    // Declare and initialize new matrix to hold updated values
    double c[N][N] = {0};
    int k, i, j;
    double a = 200;

    // Loop over the specified number of iterations
    for (k = 0; k < no_of_times; k++) {
        // Loop over all elements in the matrix
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                // Check if the current iteration is even or odd
                // and mix the matrix in the corresponding direction
                if (k % 2 == 0) { 
                    // Even iteration: mix horizontally
                    // Use bitwise AND operator to perform modulo N operation
                    // and avoid costly calls to the modulo operator
                    c[i][j] = c_old[i][(j + (int)(a * cos(i * M_PI * 2 / N))) & N - 1];
                } else {
                    // Odd iteration: mix vertically
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

