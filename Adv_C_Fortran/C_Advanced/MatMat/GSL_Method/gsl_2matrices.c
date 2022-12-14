/*
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
*/
/*
    Original Code @ https://www.gnu.org/software/gsl/doc/html/blas.html
*/
/*
Compile with:
    $ gcc -I/usr/include/gsl -lgsl -lgslcblas -lm gsl_2matrices.c -o gsl_2matrices
    $ ./gsl_2matrices
*/
#include <stdio.h>
#include <gsl/gsl_blas.h>
#define rows 2  // NB! User must specify!
#define cols 3  // NB! User must specify!

int main(void){

    double matrix[rows][rows];
    int i, j;
    
    double a[] = {  0.11, 0.12, 0.13,
                    0.21, 0.22, 0.23 };

    double b[] = {  1011, 1012,
                    1021, 1022,
                    1031, 1032 };
    
    double c[] = {  0.00, 0.00,
                    0.00, 0.00 };

    // Initialize 'matrix' elements as zeros.
    // Necessary for easier formating & printing.
    for (i = 0; i < rows; i++){ // for ROWS
        for (j = 0; j < rows; j++){ // for COLUMNS
            matrix[i][j] = 0.0;
        }
    }
/*
    double a[] = {0.11, 0.12, 0.13,
                  0.21, 0.22, 0.23,
                  0.22, 0.23, 0.33};

    double b[] = {1011, 1012, 1013,
                  1021, 1022, 1023,
                  1031, 1032, 1033};

    double c[] = {0.00, 0.00, 0.00,
                  0.00, 0.00, 0.00,
                  0.00, 0.00, 0.00};
*/  
    // NB! User must specify ROWS and COLUMNS    
    gsl_matrix_view A = gsl_matrix_view_array(a, rows, cols);
    gsl_matrix_view B = gsl_matrix_view_array(b, cols, rows);
    gsl_matrix_view C = gsl_matrix_view_array(c, rows, rows);

    /* Compute C = A B 
    The 'd' in 'dgemm' means 'double precision', 
    The 'ge' in 'dgemm' means 'general matrix', 
    The 'mm' in 'dgemm' means 'matrix-matrix' operation. 
    */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                   1.0, &A.matrix, &B.matrix,
                   0.0, &C.matrix);

    // Stores values in "matrix format".
    // P.S. Before, elements were stored in a vector.
    // c[2*dim] -> matrix[dim][dim]
    for (i = 0; i < rows; i++){
        for(j = 0; j < rows; j++){
            matrix[i][j] = c[rows*i+j];
        }
    }
    printf("\n");

    // Prints elements of matrix
    printf("Elements of matrix multiplication:\n");
    for (i = 0; i < rows; i++){
        printf("[");
        for (j = 0; j < rows; j++){
            printf(" %g ", matrix[i][j]);
        }
        printf("]\n");
    }
    
    printf("\n");
    return 0;
}
