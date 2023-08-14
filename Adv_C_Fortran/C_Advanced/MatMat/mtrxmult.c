#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

Vector  vecRand(int size);
void    printVector(Vector vec);
Vector  vecZeros(int size);
void    vecFree(Vector vec);
double  vecDot(const Vector *v1, const Vector *v2);
Vector  vecAdd(const Vector *v1, const Vector *v2);
Matrix  matMulFromVectors(const Vector *v1, const Vector *v2);

MatrixVector matVecZeros(int rows, int cols);
MatrixVector matVecRand(int rows, int cols);
void         matVecFree(MatrixVector matVec);
MatrixVector matMulVector(const MatrixVector *A, const MatrixVector *B);
void         printMatrixVector(MatrixVector matVec);
// 
// For deallocating use vecFree()

int main() {
    // Create a matrix
    Matrix matA = matRand(4, 5); // Replace with your desired dimensions
    Matrix matB = matRand(5, 4);
    Matrix matC = matZeros(4,4);

    // Create a vector
    Vector vec1 = vecRand(5);
    Vector vec2 = vecRand(4);
    //Vector vec3 = vecZeros(5, 4);

    // Create a special matrix, that is 1D
    MatrixVector matVecA = matVecRand(3, 3);
    MatrixVector matVecB = matVecRand(3, 3);
    printMatrixVector(matVecA);
    printMatrixVector(matVecB);

    // Populate the matrix with values
    /*
    for (int i = 0; i < mat.rows; i++) {
        for (int j = 0; j < mat.cols; j++) {
            mat.data[i][j] = i * mat.cols + j;
        }
    }
    */
    
    // Perform matrix multiplication
    matC                      = matMul(&matA, &matB);
    Matrix       matResult    = matMulFromVectors(&vec1, &vec2);
    MatrixVector matVecResult = matMulVector(&matVecA, &matVecB);

    // Print the matrix
    //printMatrix(matC);
    //printMatrix(matResult);
    printMatrixVector(matVecResult);


    // Clean up memory
    matFree(matA);
    matFree(matB);
    matFree(matC);

    vecFree(vec1);
    vecFree(vec2);
    matFree(matResult);

    matVecFree(matVecA);
    matVecFree(matVecB);
    matVecFree(matVecResult);

    return 0;
}


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
//====================================================================================
// Function to create a new MatrixVector and initialize with zeros
MatrixVector matVecZeros(int rows, int cols) {
    MatrixVector matVec;
    matVec.rows = rows;
    matVec.cols = cols;
    matVec.data = (double *)calloc(rows * cols, sizeof(double));

    if (matVec.data == NULL) {
        fprintf(stderr, "Memory allocation failed for MatrixVector.\n");
        matVec.rows = 0; // Reset rows in case of failure
    }

    return matVec;
}

// Function to generate a new MatrixVector with random values
MatrixVector matVecRand(int rows, int cols) {
    MatrixVector matVec = matVecZeros(rows, cols);

    for (int i = 0; i < rows * cols; i++) {
        matVec.data[i] = rand_range(0.0, 1.0);
    }

    return matVec;
}

// Function to free the memory allocated for MatrixVector
void matVecFree(MatrixVector matVec) {
    free(matVec.data);
}

// Function to perform matrix multiplication using MatrixVector
MatrixVector matMulVector(const MatrixVector *A, const MatrixVector *B) {
    if (A->cols != B->rows) {
        fprintf(stderr, "Matrix dimensions are not compatible for multiplication.\n");
        exit(EXIT_FAILURE);
    }

    MatrixVector result = matVecZeros(A->rows, B->cols);

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
void printMatrixVector(MatrixVector matVec) {
    int max_print_size = 6;

    printf("[\n");

    if (matVec.rows <= max_print_size) {
        for (int i = 0; i < matVec.rows; i++) {
            printf("  [");
            for (int j = 0; j < matVec.cols; j++) {
                printf("%f", matVec.data[i * matVec.cols + j]);
                if (j < matVec.cols - 1) printf(", ");
            }
            printf("]");
            if (i < matVec.rows - 1) printf(",\n");
        }
    } else {
        for (int i = 0; i < max_print_size; i++) {
            printf("  [");
            for (int j = 0; j < matVec.cols; j++) {
                printf("%f", matVec.data[i * matVec.cols + j]);
                if (j < matVec.cols - 1) printf(", ");
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
