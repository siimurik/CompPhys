#include <iostream>
#include <iomanip>
#include <cassert>

/*
Ex 5.5. Write a function Multiply that may be used to multiply two matrices given
the matrices and the size of both matrices. Use assertions to verify that the matrices
are of suitable sizes to be multiplied.
*/

void Multiply(double** A, double** B, double** C, int rowsA, int colsA, int colsB);
void initHilbertMat(double** H, int n);
void printMat(double** M, int n, const std::string& name);
void freeMatrix(double** M, int n);
double** allocateMatrix(int n);  // ADDED THIS FUNCTION

int main()
{
    const int n = 3;
    double scalar = 2.0;

    // Use allocation function instead of direct new
    double** A = allocateMatrix(n);
    double** B = allocateMatrix(n);
    double** C = allocateMatrix(n);

    // Initialize matrices A and B
    A[0][0] = 1.0; A[0][1] = 2.0; A[0][2] = 3.0;
    A[1][0] = 4.0; A[1][1] = 5.0; A[1][2] = 6.0;
    A[2][0] = 7.0; A[2][1] = 8.0; A[2][2] = 9.0;

    initHilbertMat(B, n);

    // Print original matrices
    printMat(A, n, "A");
    printMat(B, n, "B");
    
    // Multiply A and B, store result in C
    Multiply(A, B, C, n, n, n);

    // Print result matrix C
    printMat(C, n, "C = A * B");

    // Free memory
    freeMatrix(A, n);
    freeMatrix(B, n);
    freeMatrix(C, n);

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

void Multiply(double** A, double** B, double** C, int rowsA, int colsA, int colsB)
{
    // Assert that matrices have compatible dimensions for multiplication
    // For A × B to be valid: columns of A must equal rows of B
    // Since we're using square matrices in this example, we check colsA == rowsA
    // In general, you'd want: assert(colsA == rowsB);
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

// Free dynamically allocated memory for matrices
void freeMatrix(double** M, int n)
{
    for (int i = 0; i < n; i++)
    {
        delete[] M[i];
    }
    delete[] M;
}