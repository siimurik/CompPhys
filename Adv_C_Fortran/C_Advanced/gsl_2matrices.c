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
/* Compile and execute with:
    $ gcc -Wall -I/home/siim/gsl/include -c gsl_2matrices.c
    $ gcc -L/home/siim/gsl/lib gsl_2matrices.o -lgsl -lgslcblas -lm
    $ export LD_LIBRARY_PATH="/home/siim/gsl/lib"
    $ ./a.out 
*/
#include <stdio.h>
#include <gsl/gsl_blas.h>
#define rows 2
#define cols 3
int
main (void)
{
  	double matrix[rows][rows];
	int i, j;

	double a[] = { 	0.11, 0.12, 0.13,
        	       	0.21, 0.22, 0.23 };

  	double b[] = { 	1011, 1012,
        	       	1021, 1022,
               		1031, 1032 };

  	double c[] = { 	0.00, 0.00,
        	        0.00, 0.00 };

  	gsl_matrix_view A = gsl_matrix_view_array(a, rows, cols);
  	gsl_matrix_view B = gsl_matrix_view_array(b, cols, rows);
  	gsl_matrix_view C = gsl_matrix_view_array(c, rows, rows);

 /* Compute C = A B */

  	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
        	        1.0, &A.matrix, &B.matrix,
                	0.0, &C.matrix);

  	//printf ("[ %g, %g\n", c[0], c[1]);
  	//printf ("  %g, %g ]\n", c[2], c[3]);
	
	for(i = 0; i < rows; i++){
		for(j = 0; j < rows; j++){
		matrix[i][j] = c[i*rows+j];
		}
	}

	printf("\nResult of matrix multiplication:\n");	
	for(i = 0; i < rows; i++){
		printf("[");
		for(j = 0; j < rows; j++){
			printf(" %g ", matrix[i][j]);
		}
		printf("]\n");
	}
	printf("\n");

  	return 0;
}
