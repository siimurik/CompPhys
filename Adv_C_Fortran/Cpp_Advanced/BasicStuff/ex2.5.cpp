#include <algorithm>
#include <iostream>
#include <cstring>  // for memcpy
#include <iomanip>
#include <cmath>
#include <omp.h>    // for OpenMP

/*g++ ex2.5.cpp -o ex25 -fopenmp*/

// Fast matrix inverse using Gauss-Jordan elimination with partial pivoting
template<size_t N>
bool inverse_matrix(const double (&A)[N][N], double (&A_inv)[N][N]) {
    // Create augmented matrix [A | I]
    double augmented[N][2*N] = {0.0};
    
    // Initialize augmented matrix: left half is A, right half is identity
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            augmented[i][j] = A[i][j];
        }
        augmented[i][i + N] = 1.0;  // Identity matrix
    }
    
    // Perform Gauss-Jordan elimination
    for (size_t pivot = 0; pivot < N; pivot++) {
        // Partial pivoting: find row with largest absolute value in current column
        size_t max_row = pivot;
        double max_val = std::abs(augmented[pivot][pivot]);
        
        for (size_t row = pivot + 1; row < N; row++) {
            double abs_val = std::abs(augmented[row][pivot]);
            if (abs_val > max_val) {
                max_val = abs_val;
                max_row = row;
            }
        }
        
        // Check for singular matrix
        if (max_val < 1e-12) {
            return false;  // Matrix is singular or nearly singular
        }
        
        // Swap rows if necessary
        if (max_row != pivot) {
            for (size_t col = 0; col < 2*N; col++) {
                std::swap(augmented[pivot][col], augmented[max_row][col]);
            }
        }
        
        // Normalize the pivot row
        double pivot_val = augmented[pivot][pivot];
        for (size_t col = pivot; col < 2*N; col++) {
            augmented[pivot][col] /= pivot_val;
        }
        
        // Eliminate other rows
        for (size_t row = 0; row < N; row++) {
            if (row != pivot) {
                double factor = augmented[row][pivot];
                for (size_t col = pivot; col < 2*N; col++) {
                    augmented[row][col] -= factor * augmented[pivot][col];
                }
            }
        }
    }
    
    // Extract the inverse from the right half of augmented matrix
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            A_inv[i][j] = augmented[i][j + N];
        }
    }
    return true;  // Success
}

template<size_t Rows, size_t Cols>
void print_matrix(const double (&mat)[Rows][Cols], const std::string& name = "mat") {
    std::cout << name << ":\n";
    for (int i = 0; i < Rows; i++) {
        std::cout << "  ";
        for (int j = 0; j < Cols; j++) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(4) << mat[i][j];
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

// Ultra-fast matrix multiplication with cache blocking and OpenMP
template<size_t M, size_t N, size_t P>
void fast_matmul_blocked(const double (&A)[M][N], const double (&B)[N][P], double (&C)[M][P], 
                         size_t block_size = 64) {
    // Initialize result matrix to zero
    #pragma omp parallel for collapse(2)
    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < P; j++) {
            C[i][j] = 0.0;
        }
    }
    
    // Cache-blocked matrix multiplication
    #pragma omp parallel for collapse(2)
    for (size_t i0 = 0; i0 < M; i0 += block_size) {
        for (size_t k0 = 0; k0 < N; k0 += block_size) {
            for (size_t j0 = 0; j0 < P; j0 += block_size) {
                // Process block
                size_t i_end = std::min(i0 + block_size, M);
                size_t k_end = std::min(k0 + block_size, N);
                size_t j_end = std::min(j0 + block_size, P);
                
                for (size_t i = i0; i < i_end; i++) {
                    for (size_t k = k0; k < k_end; k++) {
                        double a_ik = A[i][k];
                        for (size_t j = j0; j < j_end; j++) {
                            C[i][j] += a_ik * B[k][j];
                        }
                    }
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////

int main(){
    double A[2][2] = {
        {4.0, 10.0},
        {1.0, 1.0}
    };

    double A_inv[2][2] = {0.0};
    inverse_matrix(A, A_inv); 
    print_matrix(A_inv);
    
    double AA[2][2] = {0.0};
    fast_matmul_blocked(A, A_inv, AA);
    print_matrix(AA);

    return 0;
}