/*
=================================================
 Compile and execute with:
    $ gcc -fopenmp mat_omp.c -o mat
    $ time OMP_NUM_THREADS=4 ./mat
=================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N 1000

void matrix_multiply(double **a, double **b, double **c) {
    #pragma omp parallel for collapse(2)
    for(int k = 0; k < N; k++){
        for(int i = 0; i < N; i++){
            c[i][k]=0;
            for(int j = 0; j < N; j++){
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

int main() {
    double **a, **b, **c;

    a = (double **)malloc(N * sizeof(double *));
    b = (double **)malloc(N * sizeof(double *));
    c = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        a[i] = (double *)malloc(N * sizeof(double));
        b[i] = (double *)malloc(N * sizeof(double));
        c[i] = (double *)malloc(N * sizeof(double));
    }

    // initialize matrices a and b with some values
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            //a[i][j] = i + j;
            //b[i][j] = i - j;
            a[i][j] = drand48();
            b[i][j] = drand48();
        }
    }

    matrix_multiply(a, b, c);

    // print the result matrix
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            printf("%lf ", c[i][j]);
        }
        printf("\n");
    }

    for (int i = 0; i < N; i++) {
        free(a[i]);
        free(b[i]);
        free(c[i]);
    }
    free(a);
    free(b);
    free(c);

    return 0;
}

/*
 Extra explanation for double allocation:
This double allocation is necessary because the matrix is 
represented as a two-dimensional array in C, but the 
memory is allocated in a one-dimensional manner.

The first allocation (malloc(N * Nof(int *))) is used 
to allocate memory for the rows of the matrix. This creates 
an array of N pointers to integers, each pointer 
representing a row of the matrix.

The second allocation (malloc(N * Nof(int))) is used 
to allocate memory for the elements of the matrix. This 
creates an array of N*N integers, representing the 
elements of the matrix.

The first allocation creates an array of pointers to 
the second allocation, where the element of the first 
allocation are pointers to the elements of the second 
allocation, that's why this sort of double allocation 
is necessary.

This way, the matrix can be accessed using the syntax 
a[i][j] where i represents the row and j represents 
the column, just like a two-dimensional array.
*/