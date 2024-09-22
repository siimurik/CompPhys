/*
If the files 'nrutil.c' and 'nrutil.h' are in the same folder
Compile and execute with
    $ gcc nr_ludecomp.c nrutil.c -o ludecomp
    $ ./ludecomp 

If not, then use the systemic linker in /usr/local/lib
    $ gcc nr_ludecomp.c -o ludecomp -lnrutil
    $ ./invmat
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "nrutil.h"
#define TINY 1.0e-20 // A small number.

void print_matrix(float **mat, int rows, int cols);
void print_vector(float *vec, int size);
void print_solution_vector(float *vec, int size);
void ludcmp(float **a, int n, int *indx, float *d);
void lubksb(float **a, int n, int *indx, float b[]);

int main() {
    int n = 3; // Example: solving a 3x3 system
    float **a, *b, d;
    int *indx;

    // Allocate memory for the matrix and vectors
    a = matrix(1, n, 1, n);    // Matrix a[1..n][1..n]
    b = vector(1, n);          // Vector b[1..n]
    indx = ivector(1, n);      // Permutation vector indx[1..n]

    // Example system of equations:
    // 3x + 2y - z = 1
    // 2x - 2y + 4z = -2
    // -x + 0.5y - z = 0

    /* NOTE: Memory is allocated at runtime, and the compiler 
    doesn't know the memory layout ahead of time. Hence, 
    it can't directly initialize the matrix with braces. */
    // Method 1: Initialize the elements manually
    // a[1][1] = 3;  a[1][2] = 2;  a[1][3] = -1;
    // a[2][1] = 2;  a[2][2] = -2; a[2][3] = 4;
    // a[3][1] = -1; a[3][2] = 0.5; a[3][3] = -1;

    // Method 2: Using Pointer Arithmetic
    static float values[] = {
         3.0,  2.0, -1.0,
         2.0, -2.0,  4.0,
        -1.0,  0.5, -1.0
    };

    // Copy the values into the matrix
    int k = 0;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            a[i][j] = values[k++];
        }
    }

    b[1] = 1;
    b[2] = -2;
    b[3] = 0;

    printf("Matrix a:\n");
    print_matrix(a, n, n);

    printf("Vector b:\n");
    print_vector(b, n);

    // Perform LU decomposition
    ludcmp(a, n, indx, &d);

    // Solve the system using back substitution
    lubksb(a, n, indx, b);
    // The answer x will be given back in b. The 
    // original matrix A will have been destroyed.

    // Output the solution
    print_solution_vector(b, n);

    // Free allocated memory
    free_matrix(a, 1, n, 1, n);
    free_vector(b, 1, n);
    free_ivector(indx, 1, n);

    return 0;
}

void print_matrix(float **mat, int rows, int cols) {
    int max_print_size = 6;

    printf("[");

    if (rows <= max_print_size && cols <= max_print_size) {
        for (int i = 1; i <= rows; i++) {  // Start from 1
            printf("[");
            for (int j = 1; j <= cols; j++) {  // Start from 1
                printf("%f", mat[i][j]);
                if (j < cols) printf(", ");
            }
            printf("]");
            if (i < rows) printf(",\n");
        }
    } else {
        for (int i = 1; i <= max_print_size; i++) {  // Start from 1
            printf("[");
            for (int j = 1; j <= max_print_size; j++) {  // Start from 1
                printf("%f", mat[i][j]);
                if (j < max_print_size) printf(", ");
            }
            printf(", ...");
            printf("]");
            if (i < max_print_size) printf(",\n");
        }
        printf(",\n  ...\n");
        printf("  ...\n");
        printf("  ...");
    }

    printf("]\n");
}

void print_vector(float *vec, int size) {
    int max_print_size = 6;

    printf("[ ");
    for (int i = 1; i <= (size <= max_print_size ? size : max_print_size); i++) {  // Start from 1 for NR-style vectors
        printf("%f", vec[i]);
        if (i < size && i < max_print_size) {
            printf(", ");
        }
    }
    if (size > max_print_size) {
        printf(", ...");  // Indicate that more elements are available but not printed
    }
    printf(" ]\n");
}

void print_solution_vector(float *vec, int size) {
    int max_print_size = 6;

    printf("Solution vector x:\n");
    for (int i = 1; i <= (size <= max_print_size ? size : max_print_size); i++) {  // Start from 1 for NR-style vectors
        printf("x[%d] = %f\n", i, vec[i]);
    }
    if (size > max_print_size) {
        printf("... (only the first %d elements printed)\n", max_print_size);  // Indicate partial print
    }
}


void ludcmp(float **a, int n, int *indx, float *d)
// Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
// permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
// indx[1..n] is an output vector that records the row permutation effected by the partial
// pivoting; d is output as ±1 depending on whether the number of row interchanges was even
// or odd, respectively. This routine is used in combination with lubksb to solve linear equations
// or invert a matrix.
{
    int i, imax, j, k;
    float big, dum, sum, temp;
    float *vv; // vv stores the implicit scaling of each row.

    vv = vector(1,n);
    *d = 1.0;               // No row interchanges yet.
    for (i=1; i<=n; i++) {  // Loop over rows to get the implicit scaling informa-
        big=0.0;            // tion.  
        for (j=1; j<=n; j++) {
            if ((temp=fabs(a[i][j])) > big) big=temp;
        }
        if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
        // No nonzero largest element.
        vv[i] = 1.0/big; // Save the scaling.
    }
    for (j=1; j<=n; j++) {      // This is the loop over columns of Crout’s method.
        for (i=1; i<j; i++) {   // This is equation (2.3.12) except for i = j.
            sum = a[i][j];
            for (k=1; k<i; k++) sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
        }
        big=0.0;                // Initialize for the search for largest pivot element.
        for (i=j;i<=n;i++) {    // This is i = j of equation (2.3.12) and i = j + 1 . . . N
            sum = a[i][j];      // of equation (2.3.13).
            for (k=1;k<j;k++) { 
                sum -= a[i][k]*a[k][j];
            }
            a[i][j]=sum;
            if ( (dum = vv[i]*fabs(sum)) >= big) {
                // Is the figure of merit for the pivot better than the best so far?
                big=dum;
                imax=i;
            }
        }
        if (j != imax) {            // Do we need to interchange rows?
            for (k=1;k<=n;k++) {    // Yes, do so...
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            *d = -(*d);     // ...and change the parity of d.
            vv[imax] = vv[j]; // Also interchange the scale factor.
        }
        indx[j] = imax;
        if (a[j][j] == 0.0) a[j][j]=TINY;
        // If the pivot element is zero the matrix is singular (at least to the precision of the
        // algorithm). For some applications on singular matrices, it is desirable to substitute
        // TINY for zero.
        if (j != n) {   // Now, finally, divide by the pivot element.
            dum=1.0/(a[j][j]);
            for (i=j+1; i<=n; i++) a[i][j] *= dum;
        }
    }                   // Go back for the next column in the reduction.
    free_vector(vv,1,n);
}

void lubksb(float **a, int n, int *indx, float b[])
// Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, not as the matrix
// A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
// as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
// B, and returns with the solution vector X. a, n, and indx are not modified by this routine
// and can be left in place for successive calls with different right-hand sides b. This routine takes
// into account the possibility that b will begin with many zero elements, so it is efficient for use
// in matrix inversion.
{
    int i, ii = 0, ip, j;
    float sum;
    for (i=1; i<=n; i++) {  // When ii is set to a positive value, it will become the
        ip    = indx[i];    // index of the first nonvanishing element of b. We now
        sum   = b[ip];      // do the forward substitution, equation (2.3.6). The
        b[ip] = b[i];       // only new wrinkle is to unscramble the permutation
        if (ii) {           // as we go.
            for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
        }
        else if (sum) ii=i; // A nonzero element was encountered, so from now on we
        b[i] = sum;         // will have to do the sums in the loop above.b[i]=sum;
    }
    for (i=n;i>=1;i--) {    // Now we do the backsubstitution, equation (2.3.7).
        sum = b[i];
        for (j=i+1; j<=n; j++) sum -= a[i][j]*b[j];
        b[i] = sum/a[i][i]; // Store a component of the solution vector X.
    } // All done!
}
