/*
5.10 Write a module for solving the n × n linear system Au = b using Gaussian
 elimination with pivoting, where A is nonsingular. See Sect.A.2.1.3 for details of
 this algorithm.
*/

#include <iostream>
#include <iomanip>
#include <cmath>

double** allocateMat(int size);
double*  allocateVec(int size);
double** readMat(double* data, int n);
double*  readVec(double* data, int n);
void printMat(double** M, int n, const std::string& name);
void printVec(double* v, int size, const std::string& name);
void GaussElimPivot(double** A, double* x, double* b, int n);
void freeMat(double** mat, int size);
void freeVec(double* v);

int main(){

    const int n = 3;
    double** A = allocateMat(n);
    double*  x = allocateVec(n);
    double*  b = allocateVec(n);

    double* A_data = new double[n*n]{
        3.0,  9.0,  6.0,
        2.0,  1.0, -1.0,
        1.0,  1.0,  1.0
    };
    A = readMat(A_data, n);
    printMat(A, n, "Matrix A");

    double* b_data = new double[n]{
        3.0,
        2.0,
        2.0
    };
    b = readVec(b_data, n);
    printVec(b, n, "Vector b");

    // Solve the linear system Ax = b
    GaussElimPivot(A, x, b, n);
    printVec(x, n, "Solution x");

    // Free memory
    freeMat(A, n);
    freeVec(x);
    freeVec(b);
    delete[] A_data;
    delete[] b_data;

    return 0;
}

double** allocateMat(int size){
    double** mat = new double*[size]; // Explanation: a vector of pointers
    for (int i = 0; i < size; i++){
        mat[i] = new double[size]; // Explanation: each pointer points to a vector of doubles
    }
    return mat;
}

// Function to allocate vector memory
double* allocateVec(int size)
{
    return new double[size];
}

void freeMat(double** mat, int size){
    for (int i = 0; i < size; i++){
        delete[] mat[i]; // Explanation: free each row
    }
    delete[] mat; // Explanation: free the array of pointers
}

void freeVec(double* v){
    delete[] v;
}

double** readMat(double* data, int n){
    double** mat = allocateMat(n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            mat[i][j] = data[i*n + j];
        }
    }
    return mat;
}

double* readVec(double* data, int n){
    double* vec = allocateVec(n);
    for (int i = 0; i < n; i++){
        vec[i] = data[i];
    }
    return vec;
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

void printVec(double* v, int size, const std::string& name)
{
    std::cout << name << ":\n";
    for (int i = 0; i < size; i++)
    {
        std::cout << std::fixed << std::setprecision(2) << std::setw(6) << v[i] << " ";
    }
    std::cout << "\n\n";
}

void GaussElimPivot(double** A, double* x, double* b, int n){
    const double EPSILON = 1e-10;
    
    // Forward elimination with partial pivoting
    for (int k = 0; k < n - 1; k++) {
        // Find pivot: row with largest absolute value in column k
        int pivotRow = k;
        double maxVal = fabs(A[k][k]);
        
        for (int i = k + 1; i < n; i++) {
            if (fabs(A[i][k]) > maxVal) {
                maxVal = fabs(A[i][k]);
                pivotRow = i;
            }
        }
        
        // Check for singular matrix
        if (maxVal < EPSILON) {
            std::cerr << "Error: Matrix is singular or nearly singular!" << std::endl;
            return;
        }
        
        // Swap rows k and pivotRow if necessary
        if (pivotRow != k) {
            std::swap(A[k], A[pivotRow]);
            std::swap(b[k], b[pivotRow]);
        }
        
        // Eliminate column entries below pivot
        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            
            // Update row i
            for (int j = k; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }
    
    // Check if last diagonal element is non-zero
    if (fabs(A[n-1][n-1]) < EPSILON) {
        std::cerr << "Error: Matrix is singular or nearly singular!" << std::endl;
        return;
    }
    
    // Back substitution
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }
}