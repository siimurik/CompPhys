/*
5.8 The determinant of a square matrix may be defined recursively: see Sect.A.1.3.
 Write a recursive function that may be used to calculate the determinant of a square
 matrix of a given size. Check the accuracy of your code by comparison with the
 known formulae for square matrices of size 2 and 3:
    det([[A00, A01], [A10, A11]]) = A00*A11 - A01*A10
 and
    det([[A00, A01, A02],
         [A10, A11, A12],
         [A20, A21, A22]]) = A00*(A11*A22 - A12*A21)
                           - A01*(A10*A22 - A12*A20)
                           + A02*(A10*A21 - A11*A20)
*/

#include <iostream>
#include <iomanip>


double   recDet(double** mat, int n);
double** allocateMat(int size);
void     freeMat(double** mat, int size);
void     printMat(double** mat, int size, const std::string& name);
double** readMat(double* data, int n);


int main(){

    // Test for 2x2 matrix 
    const int n = 2;
    double** A = allocateMat(n);
    A[0][0] = 2.0; A[0][1] = -5.0;
    A[1][0] = 7.0; A[1][1] =  3.0;

    std::cout << "Determinant of 2x2 matrix A:\n";
    printMat(A, n, "A");

    double det_A = recDet(A, n);
    std::cout << "det(A) = " << det_A << "\n";

    // Test for 3x3 matrix 
    const int m = 3;
    double** B = allocateMat(m);
    double*  b_data = new double[m*m]{
        6.0,  1.0, 1.0,
        4.0, -2.0, 5.0,
        2.0,  8.0, 7.0
    };
    
    std::cout << "\nDeterminant of 3x3 matrix B:\n";
    B = readMat(b_data, m);
    printMat(B, m, "B");
    double det_B = recDet(B, m);
    std::cout << "det(B) = " << det_B << "\n";

    // Test for 6x6 matrix
    const int p = 6;
    double** C = allocateMat(p);
    double*  c_data = new double[p*p]{
        85.0,  24.0, -14.0,  8.0, -5.0,  2.0,
        24.0,  61.0,   8.0, -4.0,  3.0, -1.0,
        -14.0,  8.0,  39.0, -2.0,  1.0,  0.0,
        8.0,   -4.0,  -2.0, 16.0, -1.0,  0.0,
        -5.0,   3.0,   1.0, -1.0,  9.0,  0.0,
        2.0,   -1.0,   0.0,  0.0,  0.0,  4.0
    };
    std::cout << "\nDeterminant of 6x6 matrix C:\n";
    C = readMat(c_data, p);
    printMat(C, p, "C");
    double det_C = recDet(C, p);
    std::cout << "det(C) = " << det_C << "\n";

    // Free memory
    freeMat(A, n);
    freeMat(B, m);
    delete[] b_data;
    freeMat(C, p);
    delete[] c_data;

    return 0;
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

double recDet(double** mat, int n){
    if (n == 1){
        return mat[0][0];
    } else if (n == 2){
        return mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
    } else {
        double det = 0.0;
        for (int p = 0; p < n; p++){
            // Create submatrix
            double** subMat = allocateMat(n - 1);
            for (int i = 1; i < n; i++){
                int subCol = 0;
                for (int j = 0; j < n; j++){
                    if (j == p) continue;
                    subMat[i - 1][subCol] = mat[i][j];
                    subCol++;
                }
            }
            // Recursive call
            det += mat[0][p] * recDet(subMat, n - 1) * 
                ((p % 2 == 0) ? 1 : -1); // Explanation: if p is even, sign is +1, else -1
            freeMat(subMat, n - 1);
        }
        return det;
    }
}

// Function to read in an array of values to create a n x n matrix
double** readMat(double* data, int n){
    double** mat = allocateMat(n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            mat[i][j] = data[i*n + j];
        }
    }
    return mat;
}