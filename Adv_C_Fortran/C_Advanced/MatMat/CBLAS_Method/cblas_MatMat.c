/*
================================================================================================
 Best possible version at the moment with multithreading enabled
    $ gcc -o mat cblas_MatMat.c -lopenblas -lpthread -lm
================================================================================================
 This code needs the BLAS library to be installed:
    $ sudo apt install libblas-dev
 On Mac:
    $ brew install openblas
================================================================================================
 Compile and execute
    $ gcc cblas_MatMat.c -o mat -lblas
    or
    $ gcc cblas_MatMat.c -o mat -lcblas -lblas
    $ ./mat
 Also it might be a good idea to use a different
 optimizer as well. Use it together with 'time' like so:
    $ gcc cblas_MatMat.c -O2 -lblas
    $ time ./a.out
================================================================================================
 Compiling on Mac:
 ---
 To copile on Mac, the main C file need two header files
 to be linked: "cblas.h" and "openblas_config.h".
 Copying them to "/usr/local/include"  with the command
    $ sudo cp /opt/homebrew/Cellar/openblas/0.3.21/include/cblas.h /usr/local/include
    $ sudo cp /opt/homebrew/Cellar/openblas/0.3.21/include/openblas_config.h /usr/local/include
 will allow you to to just use the command:
    $ gcc -o mat cblas_MatMat.c -lblas
 Otherwise you just have to use:
    $ gcc -o mat cblas_MatMat.c -I/opt/homebrew/Cellar/openblas/0.3.21/include -lblas
================================================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <time.h>

#define dim 5000
int main()
{
    // Define the dimensions of the matrices
    int m = dim, n = dim, k = dim;
    struct timespec start, stop;

    int num_threads = 8; // Specify the number of threads

    // Set the number of threads for OpenBLAS
    openblas_set_num_threads(num_threads);

    printf("Dimension of matrix matrix multiplicaton: %d x %d.\n", m,n);

    // Allocate memory for the matrices
    double *A = malloc(m * k * sizeof(double)); // 1000x1000 matrix
    double *B = malloc(k * n * sizeof(double)); // 1000x1000 matrix
    double *C = malloc(m * n * sizeof(double)); // 1000x1000 matrix
    // fill the matrices with some values
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < k; j++)
        {
            //A[i * 1000 + j] = (double) i * j;
            //B[i * 1000 + j] = (double) i + j;
	        A[i * m + j] = ((double)rand()/(double)(RAND_MAX));
	        B[i * n + j] = ((double)rand()/(double)(RAND_MAX));
	    }
    }

    // print the result (only the first few elements)
    printf("\nPrinting out the first 5x5 elements of matrix A:\n");
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            printf("%f ", A[i * m + j]);
        }
        printf("\n");
    }
    printf("\nPrinting out the first 5x5 elements of matrix B:\n");
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            printf("%f ", B[i * n + j]);
        }
        printf("\n");
    }

    // Get starting time
    clock_gettime(CLOCK_MONOTONIC, &start);

    // perform matrix multiplication using dgemm()
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, A, k, B, n, 0.0, C, n);

    // Get end time
    clock_gettime(CLOCK_MONOTONIC, &stop);
    

    // print the result (only the first few elements)
    printf("\nPrinting out the first 5x5 elements of the result matrix:\n");
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            printf("%f ", C[i * k + j]);
        }
        printf("\n");
    }
    
    // Calculate the elapsed time in seconds
    double time_taken = (stop.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (stop.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("\nMatrix multiplication took %lf seconds.\n", time_taken);

    // free the allocated memory
    free(A);
    free(B);
    free(C);

    return 0;
}

/*
In this example, dgemm() is used to multiply a m x k matrix A with a kxn matrix B 
to obtain a mxn matrix C. The function takes the following arguments:

    CblasRowMajor: This specifies the storage order of the matrices. In this case, 
    it indicates that the matrices are stored in row-major order (i.e., each row 
    is stored in contiguous memory locations).
    CblasNoTrans: This specifies whether to transpose the matrices. In this case, 
    it indicates that the matrices should not be transposed.
    'm': This is the number of rows in the result matrix C.
    'n': This is the number of columns in the result matrix C.
    'k': This is the number of columns in the matrix A (or the number of rows in the matrix B).
    '1': This is the scalar value by which to multiply the matrices A and B before adding them.
    '*A': This is a pointer to the first element of the matrix A.
    'k': This is the leading dimension (i.e., the number of columns) of the matrix A.
    '*B': This is a pointer to the first element of the matrix B.
    'n': This is the leading dimension (i.e., the number of columns) of the matrix B.
    '0': This is the scalar value to which to add the result of the matrix multiplication.
    '*C': This is a pointer to the first element of the matrix C.
    'n': This is the leading dimension (i.e., the number of columns) of the matrix C.

The matrices are allocated dynamically using malloc() instead of being declared as local 
variables on the stack. This allows the program to allocate as much memory as needed for 
the matrices, without being limited by the size of the stack. It also allows the program 
to free the memory when it is no longer needed, using free().

When using malloc(), it is important to ensure that the memory is allocated successfully 
before using it. If the allocation fails, malloc() will return NULL, and the program 
should check for this and handle the error appropriately (e.g., by printing an error 
message and exiting). It is also important to free the memory when it is no longer 
needed, to avoid memory leaks.
*/
