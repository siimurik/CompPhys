/*
 Compile and execute with:
    $ gcc mtrxmult.c -o mtrx -fopenmp -lblas -lm
 On Mac:
    $ clang -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp mtrxmult.c -o mat -lblas -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <cblas.h>

// Define the tile size for inner tiling
#define ROW_COL_PARALLEL_INNER_TILING_TILE_SIZE 64

// Define matrix dimesions
#define dim 5000

#define CLOCK_MONOTONIC 1

// Define a custom matrix structure
typedef struct {
    double **data;
    int rows;
    int cols;
} Matrix;

typedef struct {
    double *data;
    int size;
} Vector;

typedef struct {
    double *data;
    int rows;
    int cols;
} MatrixVector;


double      rand_range(double min, double max);
double    **random_matrix(int rows, int cols);
void        free_matrix(double **matrix, int rows);
void        print_matrix(double **matrix, int rows, int cols);
double    **create_zeros(int rows, int cols);

void    matFree(Matrix mat);
Matrix  matRand(int rows, int cols);
Matrix  matEmpty(int rows, int cols);
Matrix  matZeros(int rows, int cols);
void    printMatrix(Matrix mat);
Matrix  matMul(const Matrix *A, const Matrix *B);
Matrix  matCreate(int rows, int cols, const double *data);
void matmulImplRowColParallelInnerTiling(const Matrix *left,
                                         const Matrix *right,
                                         Matrix *result);

Vector  vecRand(int size);
void    printVector(Vector vec);
Vector  vecZeros(int size);
void    vecFree(Vector vec);
double  vecDot(const Vector *v1, const Vector *v2);
Vector  vecAdd(const Vector *v1, const Vector *v2);
Matrix  matMulFromVectors(const Vector *v1, const Vector *v2);
Vector  vecCreate(int size, const double *data);

MatrixVector matvecZeros(int rows, int cols);
MatrixVector matvecRand(int rows, int cols);
void         matvecFree(MatrixVector matvec);
MatrixVector matMulVector(const MatrixVector *A, const MatrixVector *B);
void         printMatrixVector(MatrixVector matvec);
MatrixVector matvecCreate(int rows, int cols, const double *data);
void         matvecMatmulParallel(const MatrixVector *left, 
                                  const MatrixVector *right, 
                                  MatrixVector *result);
void matmul(const MatrixVector *A, const MatrixVector *B, MatrixVector *C);

//=========================================================================
int main() {
    struct timespec start, stop;
    int nRows = dim; int nCols = dim;
    // Create a matrix
    /*
    Matrix matA = matRand(nRows, nCols); // Replace with your desired dimensions
    Matrix matB = matRand(nCols, nRows);
    Matrix matC = matZeros(nRows,nRows);
    printf("Printing random matrices using matRand()\n");
    printMatrix(matA);
    printMatrix(matB);
    */

    // Insert matrix values manually into a data vector
/*
    static double dataA[] = {
        0.154930,   0.590628,   0.312663,   0.363340,   0.848290,
        0.682167,   0.945112,   0.116302,   0.902817,   0.762891,
        0.270557,   0.548013,   0.847965,   0.776792,   0.093645,
        0.330188,   0.103298,   0.853670,   0.069543,   0.752094,
        0.440037,   0.519449,   0.358802,   0.617297,   0.088619
    }; 

    static double dataB[] = {
        2.8165e-01,   5.7706e-01,   3.4104e-01,   5.7045e-01,   6.4831e-01,
        7.4746e-01,   8.6455e-03,   2.2676e-01,   3.3502e-02,   8.8636e-01,
        6.2228e-01,   1.3133e-01,   8.9364e-01,   3.0900e-01,   5.1634e-01,
        7.0856e-01,   3.2182e-02,   7.4357e-01,   2.1385e-01,   1.3731e-01,
        1.5752e-02,   1.4248e-01,   5.1852e-01,   7.8274e-01,   3.1304e-01
    };

    Matrix A = matCreate(5, 5, dataA);
    Matrix B = matCreate(5, 5, dataB);

    printf("Manually inserted matrix values:\n");
    printMatrix(A);
    printMatrix(B);
*/
    /*
    // Create a vector
    Vector vec1 = vecRand(5);
    Vector vec2 = vecRand(4);
    //Vector vec3 = vecZeros(5, 4);
    printf("\nPrinting random vectors using vecRand()\n");
    printVector(vec1);
    printVector(vec2);
    */

/*
    const double vectorData[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    Vector vec = vecCreate(5, vectorData);
    printf("\nManually inserted vector values:\n");
    printVector(vec);
*/
    
    // Create a special matrix, that is 1D
    MatrixVector matvecA = matvecRand(nRows, nCols);
    MatrixVector matvecB = matvecRand(nCols, nRows);
    MatrixVector matvecC = matvecZeros(nRows, nCols);
    printf("\nPrinting random vectors using matvecRand()\n");
    printMatrixVector(matvecA);
    printMatrixVector(matvecB);
    
/*
    double matvecData[] = {
        0.11, 0.12, 0.13,
        0.21, 0.22, 0.23,
        0.31, 0.32, 0.33
    };

    MatrixVector matvec = matvecCreate(3, 3, matvecData);
    printf("\nManually inserted matrix values:\n");
    printMatrixVector(matvec);
*/
    // Populate the matrix with values
    /*
    for (int i = 0; i < mat.rows; i++) {
        for (int j = 0; j < mat.cols; j++) {
            mat.data[i][j] = i * mat.cols + j;
        }
    }
    */
    
    // Get starting time
    clock_gettime(CLOCK_MONOTONIC, &start);

    // Perform matrix multiplication
    //Matrix       matC         = matMul(&matA, &matB);
    
    // Call the matrix multiplication function
    //matmulImplRowColParallelInnerTiling(&matA, &matB, &matC);
    
    //Matrix       matResult    = matMulFromVectors(&vec1, &vec2);
    //matvecC = matMulVector(&matvecA, &matvecB);

    // Call the matrix multiplication function
    //matvecMatmulParallel(&matvecA, &matvecB, &matvecC);   // Uses the OpenMP parallel computing capabilities
    matmul(&matvecA, &matvecB, &matvecC);                   // Uses the optimised BLAS library function DGEMM()

    // Get end time
    clock_gettime(CLOCK_MONOTONIC, &stop);


    // Print the matrix
    printf("\nResult matrix:\n");
    //printMatrix(matC);
    //printMatrix(matResult);
    printMatrixVector(matvecC);

    // Calculate the elapsed time in seconds
    double time_taken = (stop.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (stop.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("\nMatrix multiplication took %.3e seconds.\n", time_taken);

    // Clean up memory
    /*
    matFree(matA);
    matFree(matB);
    matFree(matC);
    */
    /*
    vecFree(vec1);
    vecFree(vec2);
    matFree(matResult);
    */
    matvecFree(matvecA);
    matvecFree(matvecB);
    matvecFree(matvecC);
    
    return 0;
}
//=========================================================================

// Generate a random floating-point number within the given range.
double rand_range(double min, double max) {
    return min + (max-min) * ((double)rand()/RAND_MAX);
}

// Generate a random matrix with the specified number of rows and columns.
double **random_matrix(int rows, int cols) {
    // Allocate memory for an array of pointers to double arrays (rows).
    double **matrix = (double **)malloc(rows * sizeof(double *));

    if (matrix == NULL) {
        perror("Memory allocation failed");
        return NULL;
    }

    // Allocate memory for each row of the matrix.
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double *)malloc(cols * sizeof(double));

        if (matrix[i] == NULL) { // Check matrix[i], not matrix
            perror("Memory allocation failed");
            return NULL;
        }
    }
    
    // Seed the random number generator based on the current time.
    srand((unsigned int)time(NULL));

    // Fill the matrix with random numbers.
    for (int i = 0; i < rows; i++) { // Use i, not j
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = rand_range(0.0, 1.0);
        }
    }

    return matrix;
}

// Free the memory allocated for the matrix.
void free_matrix(double **matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);    // Free memory for each row.
    }
    free(matrix);           // Free memory for the array of pointers 
}

double** create_zeros(int rows, int cols) {
    double** matrix = (double**)malloc(rows * sizeof(double*));
    if (matrix == NULL) {
        fprintf(stderr, "Memory allocation failed for matrix rows.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*)calloc(cols, sizeof(double));
        if (matrix[i] == NULL) {
            fprintf(stderr, "Memory allocation failed for matrix columns.\n");
            exit(EXIT_FAILURE);
        }
    }

    return matrix;
}

// Function for printing a 2D matrix
void print_matrix(double **matrix, int rows, int cols) {
    int max_print_size = 6;

    printf("[\n");

    if (rows <= max_print_size && cols <= max_print_size) {
        for (int i = 0; i < rows; i++) {
            printf("  [");
            for (int j = 0; j < cols; j++) {
                printf("%f", matrix[i][j]);
                if (j < cols - 1) printf(", ");
            }
            printf("]");
            if (i < rows - 1) printf(",\n");
        }
    } else {
        for (int i = 0; i < max_print_size; i++) {
            printf("  [");
            for (int j = 0; j < max_print_size; j++) {
                printf("%f", matrix[i][j]);
                if (j < max_print_size - 1) printf(", ");
            }
            printf(", ...");
            printf("]");
            if (i < max_print_size - 1) printf(",\n");
        }
        printf(",\n  ...\n");
        printf("  ...\n");
        printf("  ...");
    }

    printf("\n]\n");
}

// Free the memory allocated for the matrix.
void matFree(Matrix mat) {
    for (int i = 0; i < mat.rows; i++) {
        free(mat.data[i]);    // Free memory for each row.
    }
    free(mat.data);           // Free memory for the array of pointers
}


// Generate a random matrix with the specified number of rows and columns.
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

// Function to create a new matrix
Matrix matEmpty(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++) {
        mat.data[i] = (double *)malloc(cols * sizeof(double));
    }
    return mat;
}

// Function to print a matrix
void printMatrix(Matrix mat) {
    int max_print_size = 6;

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
        printf(",\n  ...\n");
        printf("  ...\n");
        printf("  ...");
    }

    printf("\n]\n");
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


Matrix matMul(const Matrix *A, const Matrix *B) {
    if (A->cols != B->rows) {
        fprintf(stderr, "Matrix dimensions are not compatible for multiplication.\n");
        exit(EXIT_FAILURE);
    }

    Matrix result = matZeros(A->rows, B->cols);

    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->cols; j++) {
            double sum = 0.0;
            for (int k = 0; k < A->cols; k++) {
                sum += A->data[i][k] * B->data[k][j];
            }
            result.data[i][j] = sum;
        }
    }

    return result;
}

Matrix matCreate(int rows, int cols, const double *data) {
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
        mat.data[i] = (double *)malloc(cols * sizeof(double));
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

        for (int j = 0; j < cols; j++) {
            mat.data[i][j] = data[i * cols + j];
        }
    }

    return mat;
}

// Function for matrix multiplication with OpenMP parallelization and tiling
void matmulImplRowColParallelInnerTiling(const Matrix *left,
                                         const Matrix *right,
                                         Matrix *result) {
// Matrix multiplication with 5000x5000 matrices
// "> Matrix multiplication took 1.762e+02 seconds."
// Not bad...  
#pragma omp parallel for shared(result, left, right) default(none) \
    collapse(2) num_threads(8)
    for (int rowTile = 0; rowTile < left->rows; rowTile += 256) {
        for (int columnTile = 0; columnTile < right->cols; columnTile += 256) {
            for (int innerTile = 0; innerTile < left->cols; innerTile += ROW_COL_PARALLEL_INNER_TILING_TILE_SIZE) {
                for (int row = rowTile; row < rowTile + 256 && row < left->rows; row++) {
                    int innerTileEnd = innerTile + ROW_COL_PARALLEL_INNER_TILING_TILE_SIZE;
                    if (innerTileEnd > left->cols) {
                        innerTileEnd = left->cols;
                    }
                    for (int inner = innerTile; inner < innerTileEnd; inner++) {
                        for (int col = columnTile; col < columnTile + 256 && col < right->cols; col++) {
                            result->data[row][col] +=
                                left->data[row][inner] * right->data[inner][col];
                        }
                    }
                }
            }
        }
    }
}

//============================================================
Vector vecRand(int size) {
    Vector vec;
    vec.size = size;
    vec.data = (double *)malloc(size * sizeof(double));

    if (vec.data == NULL) {
        perror("Memory allocation failed");
        vec.size = 0; // Reset size in case of failure
        return vec;
    }

    for (int i = 0; i < size; i++) {
        vec.data[i] = rand_range(0.0, 1.0);
    }

    return vec;
}

void printVector(Vector vec) {
    printf("[");
    for (int i = 0; i < vec.size; i++) {
        printf("%f", vec.data[i]);
        if (i < vec.size - 1) printf(", ");
    }
    printf("]\n");
}

Vector vecZeros(int size) {
    Vector vec;
    vec.size = size;
    vec.data = (double *)calloc(size, sizeof(double));

    if (vec.data == NULL) {
        fprintf(stderr, "Memory allocation failed for vector.\n");
        vec.size = 0; // Reset size in case of failure
    }

    return vec;
}

void vecFree(Vector vec) {
    free(vec.data);
}

double vecDot(const Vector *v1, const Vector *v2) {
    if (v1->size != v2->size) {
        fprintf(stderr, "Vector dimensions are not compatible for dot product.\n");
        exit(EXIT_FAILURE);
    }

    double result = 0.0;
    for (int i = 0; i < v1->size; i++) {
        result += v1->data[i] * v2->data[i];
    }

    return result;
}

Vector vecAdd(const Vector *v1, const Vector *v2) {
    if (v1->size != v2->size) {
        fprintf(stderr, "Vector dimensions are not compatible for vector addition.\n");
        exit(EXIT_FAILURE);
    }

    Vector result = vecZeros(v1->size);
    for (int i = 0; i < v1->size; i++) {
        result.data[i] = v1->data[i] + v2->data[i];
    }

    return result;
}

// Function to perform matrix multiplication based on two vectors
Matrix matMulFromVectors(const Vector *v1, const Vector *v2) {
    Matrix result = matEmpty(v1->size, v2->size);

    for (int i = 0; i < v1->size; i++) {
        for (int j = 0; j < v2->size; j++) {
            result.data[i][j] = v1->data[i] * v2->data[j];
        }
    }

    return result;
}

Vector vecCreate(int size, const double *data) {
    Vector vec;
    vec.size = size;

    vec.data = (double *)malloc(size * sizeof(double));
    if (vec.data == NULL) {
        fprintf(stderr, "Memory allocation failed for vector.\n");
        vec.size = 0; // Reset size in case of failure
        return vec;
    }

    for (int i = 0; i < size; i++) {
        vec.data[i] = data[i];
    }

    return vec;
}

//====================================================================================
// Function to create a new MatrixVector and initialize with zeros
MatrixVector matvecZeros(int rows, int cols) {
    MatrixVector matvec;
    matvec.rows = rows;
    matvec.cols = cols;
    matvec.data = (double *)calloc(rows * cols, sizeof(double));

    if (matvec.data == NULL) {
        fprintf(stderr, "Memory allocation failed for MatrixVector.\n");
        matvec.rows = 0; // Reset rows in case of failure
    }

    return matvec;
}

// Function to generate a new MatrixVector with random values
MatrixVector matvecRand(int rows, int cols) {
    MatrixVector matvec = matvecZeros(rows, cols);

    for (int i = 0; i < rows * cols; i++) {
        matvec.data[i] = rand_range(0.0, 1.0);
    }

    return matvec;
}

// Function to free the memory allocated for MatrixVector
void matvecFree(MatrixVector matvec) {
    free(matvec.data);
}

// Function to perform matrix multiplication using MatrixVector
MatrixVector matMulVector(const MatrixVector *A, const MatrixVector *B) {
    if (A->cols != B->rows) {
        fprintf(stderr, "Matrix dimensions are not compatible for multiplication.\n");
        exit(EXIT_FAILURE);
    }

    MatrixVector result = matvecZeros(A->rows, B->cols);

    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->cols; j++) {
            double sum = 0.0;
            for (int k = 0; k < A->cols; k++) {
                sum += A->data[i * A->cols + k] * B->data[k * B->cols + j];
            }
            result.data[i * result.cols + j] = sum;
        }
    }

    return result;
}


// Function to print a MatrixVector
void printMatrixVector(MatrixVector matvec) {
    int max_print_size = 6;

    printf("[\n");

    if (matvec.rows <= max_print_size) {
        for (int i = 0; i < matvec.rows; i++) {
            printf("  [");
            for (int j = 0; j < matvec.cols; j++) {
                printf("%f", matvec.data[i * matvec.cols + j]);
                if (j < matvec.cols - 1) printf(", ");
            }
            printf("]");
            if (i < matvec.rows - 1) printf(",\n");
        }
    } else {
        for (int i = 0; i < max_print_size; i++) {
            printf("  [");
            for (int j = 0; j < max_print_size; j++) {
                printf("%f", matvec.data[i * matvec.cols + j]);
                if (j < matvec.cols - 1) printf(", ");
            }
            printf(" ...");
            printf("]");
            if (i < max_print_size - 1) printf(",\n");
        }
        printf(",\n  ...\n");
        printf("  ...\n");
        printf("  ...");
    }

    printf("\n]\n");
}

MatrixVector matvecCreate(int rows, int cols, const double *data) {
    MatrixVector matvec;
    matvec.rows = rows;
    matvec.cols = cols;

    matvec.data = (double *)malloc(rows * cols * sizeof(double));
    if (matvec.data == NULL) {
        fprintf(stderr, "Memory allocation failed for matrix vector.\n");
        matvec.rows = 0; // Reset rows in case of failure
        matvec.cols = 0; // Reset cols in case of failure
        return matvec;
    }

    for (int i = 0; i < rows * cols; i++) {
        matvec.data[i] = data[i];
    }

    return matvec;
}

// Matrix multiplication with 5000x5000 matrices
// "> Matrix multiplication took 1.732e+02 seconds."
// Not bad...  at all...
void matvecMatmulParallel(const MatrixVector *left, const MatrixVector *right, MatrixVector *result) {
    // Parallelize the matrix multiplication loop using OpenMP
    #pragma omp parallel for shared(result, left, right) default(none) \
        collapse(2) num_threads(8)
    // Divide the computation into tiles for row, column, and inner dimension
    for (int rowTile = 0; rowTile < left->rows; rowTile += 256) {
        for (int columnTile = 0; columnTile < right->cols; columnTile += 256) {
            for (int innerTile = 0; innerTile < left->cols; innerTile += ROW_COL_PARALLEL_INNER_TILING_TILE_SIZE) {
                // Perform matrix multiplication within the specified tiles
                for (int row = rowTile; row < rowTile + 256 && row < left->rows; row++) {
                    // Determine the end of the inner tile
                    int innerTileEnd = innerTile + ROW_COL_PARALLEL_INNER_TILING_TILE_SIZE;
                    if (innerTileEnd > left->cols) {
                        innerTileEnd = left->cols;
                    }
                    // Perform the inner loop of matrix multiplication
                    for (int inner = innerTile; inner < innerTileEnd; inner++) {
                        for (int col = columnTile; col < columnTile + 256 && col < right->cols; col++) {
                            // Accumulate the result by multiplying elements and summing up
                            result->data[row * result->cols + col] +=
                                left->data[row * left->cols + inner] * right->data[inner * right->cols + col];
                        }
                    }
                }
            }
        }
    }
}

void matmul(const MatrixVector *A, const MatrixVector *B, MatrixVector *C) {
    // Check if matrix dimensions are compatible for multiplication
    if (A->cols != B->rows || A->rows != C->rows || B->cols != C->cols) {
        fprintf(stderr, "ERROR: Incorrect matrix sizes for multiplication!\n");
        return;
    }

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                A->rows, B->cols, A->cols,
                1.0, A->data, A->cols, B->data, B->cols,
                0.0, C->data, C->cols);
}
