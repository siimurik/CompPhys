#include <iostream>

int main() {
    const long long ITERATIONS = 1000000000;
    
    for (long long iter = 0; iter < ITERATIONS; iter++) {
        // Allocate as 1D arrays of 4 elements (2x2 = 4)
        double* A = new double[4];
        double* B = new double[4];
        double* C = new double[4];
        
        // Assign values using row-major order
        // A = [1 2]   B = [5 6]
        //     [3 4]       [7 8]
        A[0] = 1.0; A[1] = 2.0; A[2] = 3.0; A[3] = 4.0;
        B[0] = 5.0; B[1] = 6.0; B[2] = 7.0; B[3] = 8.0;
        
        // Calculate C = A + B
        for (int i = 0; i < 4; i++) {
            C[i] = A[i] + B[i];
        }
        
        // Print C (uncomment to see results)
        /*
        std::cout << "Matrix C = A + B:" << "\n";
        std::cout << C[0] << " " << C[1] << "\n";
        std::cout << C[2] << " " << C[3] << "\n";
        std::cout << "\n";
        */
        
        // De-allocate memory
        delete[] A;
        delete[] B;
        delete[] C;
    }
    
    std::cout << "Completed " << ITERATIONS 
    << " iterations without memory leaks!" 
    << '\n';
    return 0;
}