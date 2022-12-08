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
#define N 1024
#define M_PI 3.14159265358979323846

// Changed the function prototypes to match the definitions below.
void p_diffuse(double **matrix, int rows, int cols, bool in_x_direction);
void iterate(int no_of_times, double c_old[][N]);
void matrix_csv(char *filename,double a[][N],int n,int m);

int main() {
    double a = 200;

    static double c1[N][N]; // Changed the type to static.
    int i, j;

    for (i = 0; i < N; i++) { //same code that was in lecture: initial sine
        for (j = 0; j < N; j++) {
            c1[i][j] = sin(j * M_PI * 2 / N); //sine on i-axis
        }
    }

    iterate(11, c1); // Cast c1 to a double** before passing it as an argument.
    // print the result (only the first few elements)
    //for (i = 0; i < 5; i++){
    //    for (j = 0; j < 5; j++){
    //        printf("%f ", c1[i][j]);
    //    }
    //    printf("\n");
    //}
    char str[20] = "colors_test";
    matrix_csv(str, c1, N, N); // Cast c1 to a double** before passing it as an argument.

    return 0;
}


//---------------------------------------------------------------------------------
void p_diffuse(double **matrix, int rows, int cols, bool in_x_direction) {
    double c[rows][cols]; // Changed the type to double.

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double sum = 0; // Changed the type to double.
            double count = 0; // Changed the type to double.
            if (i > 0) {
                sum += matrix[i - 1][j];
                count++;
            }
            if (i < rows - 1) {
                sum += matrix[i + 1][j];
                count++;
            }
            if (in_x_direction) {
                if (j > 0) {
                    sum += matrix[i][j - 1];
                    count++;
                }
                if (j < cols - 1) {
                    sum += matrix[i][j + 1];
                    count++;
                }
            }
            c[i][j] = sum / count;
        }
    }
    memcpy(matrix, c, sizeof(double) * rows * cols); // Added a call to memcpy.
}

//---------------------------------------------------------------------------------
void iterate(int no_of_times, double c_old[][N]) {
    // Declare and initialize new matrix to hold updated values
    double c[N][N];
    double a = 200; // Moved the declaration and initialization outside the for loop.
    int k, i, j;

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
                    c[i][j] = c_old[(i + (int)(a * sin(j * M_PI * 2 / N))) & N - 1][j];
                }
            }
        }
        // Copy updated values back to original matrix
        memcpy(c_old, c, sizeof(double) * N * N);
    }
}


void matrix_csv(char *filename,double a[][N],int n,int m){
    printf("\n Creating %s.csv file for matrix.",filename);
    FILE *fp;
    int i,j;
    // Create a local variable to store the concatenated string.
    char filename_with_ext[100];
    sprintf(filename_with_ext, "%s.csv", filename);
    fp=fopen(filename_with_ext, "w+");
    // Cast the a argument to double**.
    for(i = 0; i < n; i++){
        for(j = 0; j < m; j++){
            // Simply print the value of the element as a double.
            fprintf(fp,"%f,", a[i][j]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    printf("\n Done.");
}

