/*
5.9 Write a module for solving the 3 × 3 linear system Au = b where A is nonsingular.
*/

#include <iostream>
#include <iomanip>
#include <cmath>

double** allocateMat(int size);
double*  allocateVec(int size);
double** readMat(double* data, int n);
void SolveLinearSys(double** A, double* x, double* b, int n);
void printVec(double* v, int size, const std::string& name);
void printMat(double** M, int n, const std::string& name);
void freeMat(double** mat, int size);
void freeVec(double* v);


int main(){

    const int n = 4;
    double** A = allocateMat(n);
    double*  x = allocateVec(n);
    double*  b = allocateVec(n);

    // Initialize matrix A (symmetric, positive-definite)
    double* a_data = new double[n*n]{
        4.0,  2.0,  1.0,  1.0,
        2.0,  5.0,  1.0,  2.0,
        1.0,  1.0,  6.0,  2.0,
        1.0,  2.0,  2.0,  7.0
    };
    
    A = readMat(a_data, n);
    printMat(A, n, "Matrix A");

    // Initialize vector b
    b[0] = 8.0;
    b[1] = 10.0;
    b[2] = 10.0;
    b[3] = 12.0;

    printVec(b, n, "Vector b");

    // Solve the linear system Ax = b
    SolveLinearSys(A, x, b, n);
    printVec(x, n, "Solution x");

    // Free memory
    freeMat(A, n);
    freeVec(x);
    freeVec(b);
    delete[] a_data;

    return 0;
}

double** readMat(double* data, int n){
    double** mat = new double*[n];
    for (int i = 0; i < n; i++){
        mat[i] = new double[n];
        for (int j = 0; j < n; j++){
            mat[i][j] = data[i*n + j];
        }
    }
    return mat;
}

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

double** allocateMat(int size){
    double** mat = new double*[size]; // Explanation: a vector of pointers
    for (int i = 0; i < size; i++){
        mat[i] = new double[size]; // Explanation: each pointer points to a vector of doubles
    }
    return mat;
}

void freeMat(double** mat, int size){
    for (int i = 0; i < size; i++){
        delete[] mat[i]; // Explanation: free each row
    }
    delete[] mat; // Explanation: free the array of pointers
}

// Function to allocate vector memory
double* allocateVec(int size)
{
    return new double[size];
}

// Free dynamically allocated memory for vectors
void freeVec(double* v)
{
    delete[] v;
}

// Print vector
void printVec(double* v, int size, const std::string& name)
{
    std::cout << name << ":\n";
    for (int i = 0; i < size; i++)
    {
        std::cout << std::fixed << std::setprecision(2) << std::setw(6) << v[i] << " ";
    }
    std::cout << "\n\n";
}

// Cholesky decomposition: A = L * L^T
// Returns true if decomposition successful, false if matrix is not positive definite
bool CholeskyDecomposition(double** A, double** L, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            L[i][j] = 0.0;
        }
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            
            if (j == i) {
                // Diagonal elements
                for (int k = 0; k < j; k++) {
                    sum += L[j][k] * L[j][k];
                }
                double val = A[j][j] - sum;
                if (val <= 0.0) {
                    std::cerr << "Matrix is not positive definite!" << "\n";
                    return false;
                }
                L[j][j] = sqrt(val);
            } else {
                // Off-diagonal elements
                for (int k = 0; k < j; k++) {
                    sum += L[i][k] * L[j][k];
                }
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
    }
    return true;
}

// Forward substitution: Solve L*y = b for y
void ForwardSubstitution(double** L, double* y, double* b, int n) {
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / L[i][i];
    }
}

// Backward substitution: Solve L^T*x = y for x
void BackwardSubstitution(double** L, double* x, double* y, int n) {
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += L[j][i] * x[j];
        }
        x[i] = (y[i] - sum) / L[i][i];
    }
}

// Solve linear systems of the form Ax=b, where A is a symmetric positive-definite matrix
void SolveLinearSys(double** A, double* x, double* b, int n) {
    // Allocate memory for L matrix
    double** L = allocateMat(n);
    
    // Allocate memory for intermediate vector y
    double* y = allocateVec(n);
    
    // Perform Cholesky decomposition: A = L * L^T
    if (!CholeskyDecomposition(A, L, n)) {
        // Cleanup and return if decomposition fails
        freeMat(L, n);
        freeVec(y);
        return;
    }
    
    // Solve L*y = b using forward substitution
    ForwardSubstitution(L, y, b, n);
    
    // Solve L^T*x = y using backward substitution
    BackwardSubstitution(L, x, y, n);
    
    // Cleanup
    freeMat(L, n);
    freeVec(y);
}