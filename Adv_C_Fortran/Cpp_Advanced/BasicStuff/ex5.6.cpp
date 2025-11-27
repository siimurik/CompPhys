#include <iostream>
#include <iomanip>
#include <cassert>

/*
Ex 5.6. Overload the function Multiply written in the previous exercise so that it may
be used to multiply:
1. a vector and a matrix of given sizes;
2. a matrix and a vector of given sizes;
3. a scalar and a matrix of a given size; and
4. a matrix of a given size and a scalar.
*/

// Function overloading is pretty based tbh
void Multiply(double** A, double** B, double** C, int rowsA, int colsA, int colsB);
void Multiply(double* vec, double** mat, double* result, int vecSize, int matCols);
void Multiply(double** mat, double* vec, double* result, int matRows, int vecSize);
void Multiply(double scalar, double** mat, double** result, int rows, int cols);
void Multiply(double** mat, double scalar, double** result, int rows, int cols);

void printMat(double** M, int n,   const std::string& name);
void printVec(double* v, int size, const std::string& name);
void initHilbertMat(double** H, int n);
void freeMatrix(double** M, int n);
void freeVector(double* v);
double** allocateMatrix(int n);
double*  allocateVector(int size);

int main()
{
    const int n = 3;
    double scalar = 2.0;

    // Allocate matrices A, B, C
    double** A = allocateMatrix(n);
    double** B = allocateMatrix(n);
    double** C = allocateMatrix(n);

    // Allocate vectors
    double* u = allocateVector(n);
    double* v = allocateVector(n);

    // Initialize vectors
    for (int i = 0; i < n; i++){
        u[i] = (i+1) * 2.0;
    }

    // Initialize matrices A and B
    A[0][0] = 1.0; A[0][1] = 2.0; A[0][2] = 3.0;
    A[1][0] = 4.0; A[1][1] = 5.0; A[1][2] = 6.0;
    A[2][0] = 7.0; A[2][1] = 8.0; A[2][2] = 9.0;

    initHilbertMat(B, n);

    // Print original matrices
    printMat(A, n, "A");
    printMat(B, n, "B");
    
    // Multiply matrices A and B, store result in C
    Multiply(A, B, C, n, n, n);

    // Print result matrix C
    printMat(C, n, "C = A * B");

    // 1. Test vector × matrix multiplication
    printVec(u, n, "u");
    Multiply(u, B, v, n, n); // v = u × B
    printVec(v, n, "v = u * B");

    // 2. Test matrix × vector multiplication
    Multiply(A, u, v, n, n); // v = A × u
    printVec(v, n, "v = A * u");

    // 3. Test scalar × matrix multiplication
    Multiply(scalar, A, C, n, n); // C = scalar × A
    printMat(C, n, "C = scalar * A");

    // 4. Test matrix × scalar multiplication
    Multiply(B, scalar, C, n, n); // C = B × scalar
    printMat(C, n, "C = B * scalar");

    // Free memory
    freeMatrix(A, n);
    freeMatrix(B, n);
    freeMatrix(C, n);
    freeVector(u);
    freeVector(v);

    return 0;
}

// Function to properly allocate matrix memory
double** allocateMatrix(int n)
{
    double** M = new double*[n];
    for (int i = 0; i < n; i++) {
        M[i] = new double[n];  
    }
    return M;
}

// Function to allocate vector memory
double* allocateVector(int size)
{
    return new double[size];
}

// Matrix matrix multiplication C = A × B
void Multiply(double** A, double** B, double** C, int rowsA, int colsA, int colsB)
{
    // Assert that matrices have compatible dimensions for multiplication
    // For A × B to be valid: columns of A must equal rows of B
    // Since we're using square matrices in this example, we check colsA == rowsA
    // In general, we would want: assert(colsA == rowsB);
    assert(colsA == rowsA && "Matrix dimensions must be compatible for multiplication");
    // Perform matrix multiplication: C[i][j] = Σ A[i][k] * B[k][j]
    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < colsB; j++) {
            for (int k = 0; k < colsA; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// 1. Vector × Matrix multiplication: result = vector × matrix
void Multiply(double* vec, double** mat, double* result, int vecSize, int matCols)
{
    assert(vecSize > 0 && matCols > 0 && "Invalid dimensions for vector-matrix multiplication");
    for (int j = 0; j < matCols; j++) {
        result[j] = 0.0;
        for (int i = 0; i < vecSize; i++) {
            result[j] += vec[i] * mat[i][j];
        }
    }
}

// 2. Matrix × Vector multiplication: result = matrix × vector
void Multiply(double** mat, double* vec, double* result, int matRows, int vecSize)
{
    assert(matRows > 0 && vecSize > 0 && "Invalid dimensions for matrix-vector multiplication");
    for (int i = 0; i < matRows; i++) {
        result[i] = 0.0;
        for (int j = 0; j < vecSize; j++) {
            result[i] += mat[i][j] * vec[j];
        }
    }
}

// 3. Scalar × Matrix multiplication: result = scalar × matrix
void Multiply(double scalar, double** mat, double** result, int rows, int cols)
{
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[i][j] = scalar * mat[i][j];
        }
    }
}

// 4. Matrix × Scalar multiplication: result = matrix × scalar
void Multiply(double** mat, double scalar, double** result, int rows, int cols)
{
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[i][j] = mat[i][j] * scalar;
        }
    }
}

// Function to initialize a Hilbert matrix
void initHilbertMat(double** H, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            H[i][j] = 1.0 / (i + j + 1);
        }
    }
}

// Function to print a matrix in a pretty and formatted to 2 decimal places way
void printMat(double** M, int n, const std::string& name)
{
    std::cout << name << ":\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << std::fixed << std::setprecision(2) << std::setw(6) << M[i][j] << " ";
        }
        std::cout << '\n';
    }
    std::cout << '\n';
}

// Function to print a vector in a pretty and formatted to 2 decimal places way
void printVec(double* v, int size, const std::string& name)
{
    std::cout << name << ":\n";
    for (int i = 0; i < size; i++)
    {
        std::cout << std::fixed << std::setprecision(2) << std::setw(6) << v[i] << " ";
    }
    std::cout << "\n\n";
}

// Free dynamically allocated memory for matrices
void freeMatrix(double** M, int n)
{
    for (int i = 0; i < n; i++)
    {
        delete[] M[i];
    }
    delete[] M;
}

// Free dynamically allocated memory for vectors
void freeVector(double* v)
{
    delete[] v;
}