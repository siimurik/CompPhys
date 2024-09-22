#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
    float *v;
    v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
    if (!v) nrerror("allocation failure in vector()");
    return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
    int *v;
    v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
    if (!v) nrerror("allocation failure in ivector()");
    return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
    unsigned char *v;
    v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
    if (!v) nrerror("allocation failure in cvector()");
    return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
    unsigned long *v;
    v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
    if (!v) nrerror("allocation failure in lvector()");
    return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
    double *v;
    v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
    if (!v) nrerror("allocation failure in dvector()");
    return v-nl+NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    float **m;
    /* allocate pointers to rows */
    m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;
    /* allocate rows and set pointers to them */
    m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
    if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    double **m;
    /* allocate pointers to rows */
    m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;
    /* allocate rows and set pointers to them */
    m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
    if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    int **m;
    /* allocate pointers to rows */
    m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;
    /* allocate rows and set pointers to them */
    m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
    if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
    long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
    float **m;
    /* allocate array of pointers to rows */
    m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
    if (!m) nrerror("allocation failure in submatrix()");
    m += NR_END;
    m -= newrl;
    /* set pointers to rows */
    for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
    long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
    float **m;
    /* allocate pointers to rows */
    m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
    if (!m) nrerror("allocation failure in convert_matrix()");
    m += NR_END;
    m -= nrl;
    /* set pointers to rows */
    m[nrl]=a-ncl;
    for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
    float ***t;
    /* allocate pointers to pointers to rows */
    t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
    if (!t) nrerror("allocation failure 1 in f3tensor()");
    t += NR_END;
    t -= nrl;
    /* allocate pointers to rows and set pointers to them */
    t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
    if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;
    /* allocate rows and set pointers to them */
    t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
    if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;
    for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
    for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
    }
    /* return pointer to array of pointers to rows */
    return t;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
    free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
    free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
    free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
    free((FREE_ARG) (t[nrl]+ncl-NR_END));
    free((FREE_ARG) (t+nrl-NR_END));
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
    float TINY = 1.0e-20;

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

void inv(float **a, int n, float **y, int *indx, float d, float *col) 
// The matrix y will now contain the inverse of the 
// original matrix a, which will have been destroyed.
{
    ludcmp(a, n, indx, &d);         // Decompose the matrix just once.
    for (int j = 1; j <= n; j++) {  // Find inverse by columns.
        for (int i = 1; i <= n; i++) col[i] = 0.0;
        col[j] = 1.0;
        lubksb(a, n, indx, col);
        for (int i = 1; i <= n; i++) y[i][j] = col[i];
    }
}

float det(float **a, int n, int *indx) {
    float d = 1.0;

    ludcmp(a, n, indx, &d);  // This returns d as ±1.

    for (int j = 1; j <= n; j++) {  // Start from 1 instead of 0
        d *= a[j][j];
    }

    return d;
}