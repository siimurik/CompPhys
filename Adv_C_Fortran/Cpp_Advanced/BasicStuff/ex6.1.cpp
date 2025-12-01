// ============================================================================
// Helper functions for matrix operations
// ============================================================================
// Compile with: g++ ex6.1.cpp ComplexNumber.cpp -o ex61.exe -Wall
// ============================================================================
#include "ComplexNumber.hpp"
#include <iostream>
#include <cmath>
#include <windows.h>

// Function to allocate a 3x3 matrix of complex numbers
ComplexNumber** AllocateComplexMatrix(int n) {
    ComplexNumber** matrix = new ComplexNumber*[n];
    for (int i = 0; i < n; i++) {
        matrix[i] = new ComplexNumber[n];
    }
    return matrix;
}

// Function to deallocate a matrix
void DeallocateComplexMatrix(ComplexNumber** matrix, int n) {
    for (int i = 0; i < n; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

// Function to multiply two complex matrices
void MatrixMultiply(ComplexNumber** A, ComplexNumber** B, ComplexNumber** result, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = ComplexNumber(0.0, 0.0);
            for (int k = 0; k < n; k++) {
                result[i][j] = result[i][j] + A[i][k] * B[k][j];
            }
        }
    }
}

// Function to copy a matrix
void CopyMatrix(ComplexNumber** source, ComplexNumber** dest, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dest[i][j] = source[i][j];
        }
    }
}

// Function to create an identity matrix
void CreateIdentityMatrix(ComplexNumber** I, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                I[i][j] = ComplexNumber(1.0, 0.0);
            } else {
                I[i][j] = ComplexNumber(0.0, 0.0);
            }
        }
    }
}

// Function to compute factorial
double Factorial(int n) {
    if (n <= 1) return 1.0;
    double result = 1.0;
    for (int i = 2; i <= n; i++) {
        result *= i;
    }
    return result;
}

// Function to scale a matrix by a scalar
void ScaleMatrix(ComplexNumber** matrix, double scalar, ComplexNumber** result, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = matrix[i][j] * ComplexNumber(scalar, 0.0);
        }
    }
}

// Function to add two matrices
void AddMatrices(ComplexNumber** A, ComplexNumber** B, ComplexNumber** result, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
}

// Function to compute matrix exponential exp(A) = sum(A^n / n!)
void MatrixExponential(ComplexNumber** A, ComplexNumber** expA, int n, int maxTerms = 20) {
    // Initialize expA as identity matrix (term for n=0)
    CreateIdentityMatrix(expA, n);
    
    // Allocate temporary matrices
    ComplexNumber** powerA = AllocateComplexMatrix(n);
    ComplexNumber** temp = AllocateComplexMatrix(n);
    ComplexNumber** scaledTerm = AllocateComplexMatrix(n);
    
    // Initialize powerA as identity (A^0 = I)
    CreateIdentityMatrix(powerA, n);
    
    // Sum the series: exp(A) = I + A + A^2/2! + A^3/3! + ...
    for (int k = 1; k <= maxTerms; k++) {
        // Compute A^k = A^(k-1) * A
        MatrixMultiply(powerA, A, temp, n);
        CopyMatrix(temp, powerA, n);
        
        // Scale by 1/k!
        ScaleMatrix(powerA, 1.0 / Factorial(k), scaledTerm, n);
        
        // Add to result
        AddMatrices(expA, scaledTerm, temp, n);
        CopyMatrix(temp, expA, n);
    }
    
    // Clean up temporary matrices
    DeallocateComplexMatrix(powerA, n);
    DeallocateComplexMatrix(temp, n);
    DeallocateComplexMatrix(scaledTerm, n);
}

// Function to print a complex matrix
void PrintComplexMatrix(ComplexNumber** matrix, int n, const std::string& name) {
    std::cout << name << ":\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << matrix[i][j] << "  ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

int main(){

    // Set console to UTF-8 if the OS is Windows or Linux
    #ifdef _WIN32
        SetConsoleOutputCP(CP_UTF8);
    #elif __linux__
        std::cout << "\033%G"; // Switch to UTF-8 in Linux terminal
    #endif

    std::cout << "===========================================\n";
    std::cout << "      Testing ComplexNumber class          \n";
    std::cout << "===========================================\n";
    
    // Test 1: Constructor with real and imaginary parts + Getter methods
    std::cout << "\nTest 1: Constructor and Getter Methods\n";
    std::cout << "---------------------------------------\n";
    ComplexNumber z1(3.0, 4.0);
    std::cout << "z1 = " << z1 << "\n";
    std::cout << "Real part (using GetRealPart): " << z1.GetRealPart() << "\n";
    std::cout << "Imaginary part (using GetImaginaryPart): " << z1.GetImaginaryPart() << "\n";
    
    // Test 2: Friend functions RealPart and ImaginaryPart
    std::cout << "\nTest 2: Friend Functions\n";
    std::cout << "------------------------\n";
    ComplexNumber z2(5.0, -2.0);
    std::cout << "z2 = " << z2 << "\n";
    std::cout << "Real part (using RealPart function): " << RealPart(z2) << "\n";
    std::cout << "Imaginary part (using ImaginaryPart function): " << ImaginaryPart(z2) << "\n";
    std::cout << "Verify: z2.GetRealPart() == RealPart(z2)? " 
              << (z2.GetRealPart() == RealPart(z2) ? "YES" : "NO") << "\n";
    
    // Test 3: Copy constructor
    std::cout << "\nTest 3: Copy Constructor\n";
    std::cout << "------------------------\n";
    ComplexNumber z3(2.0, 3.0);
    ComplexNumber z4(z3);  // Using copy constructor
    std::cout << "Original z3 = " << z3 << "\n";
    std::cout << "Copy z4 (created with copy constructor) = " << z4 << "\n";
    std::cout << "Are they equal? " 
              << (z3.GetRealPart() == z4.GetRealPart() && 
                  z3.GetImaginaryPart() == z4.GetImaginaryPart() ? "YES" : "NO") << "\n";
    
    // Test 4: Single argument constructor (real number)
    std::cout << "\nTest 4: Single Argument Constructor\n";
    std::cout << "------------------------------------\n";
    ComplexNumber z5(7.5);  // Should create 7.5 + 0i
    std::cout << "z5(7.5) = " << z5 << "\n";
    std::cout << "Real part: " << z5.GetRealPart() << "\n";
    std::cout << "Imaginary part: " << z5.GetImaginaryPart() << "\n";
    std::cout << "Is imaginary part zero? " 
              << (z5.GetImaginaryPart() == 0.0 ? "YES" : "NO") << "\n";
    
    // Test 5: Conjugate method
    std::cout << "\nTest 5: Conjugate Method\n";
    std::cout << "------------------------\n";
    ComplexNumber z6(3.0, 4.0);
    ComplexNumber z6_conj = z6.CalculateConjugate();
    std::cout << "z6 = " << z6 << "\n";
    std::cout << "Conjugate of z6 = " << z6_conj << "\n";
    std::cout << "Real parts equal? " 
              << (z6.GetRealPart() == z6_conj.GetRealPart() ? "YES" : "NO") << "\n";
    std::cout << "Imaginary parts opposite? " 
              << (z6.GetImaginaryPart() == -z6_conj.GetImaginaryPart() ? "YES" : "NO") << "\n";
    
    // Test conjugate of conjugate equals original
    ComplexNumber z6_conj_conj = z6_conj.CalculateConjugate();
    std::cout << "Conjugate of conjugate = " << z6_conj_conj << "\n";
    std::cout << "Equals original? " 
              << (z6.GetRealPart() == z6_conj_conj.GetRealPart() && 
                  z6.GetImaginaryPart() == z6_conj_conj.GetImaginaryPart() ? "YES" : "NO") << "\n";
    
    // Test 6: Mathematical property - z + conjugate(z) = 2*Re(z)
    std::cout << "\nTest 6: Mathematical Property\n";
    std::cout << "------------------------------\n";
    ComplexNumber z7(6.0, 8.0);
    ComplexNumber z7_conj = z7.CalculateConjugate();
    ComplexNumber sum = z7 + z7_conj;
    std::cout << "z7 = " << z7 << "\n";
    std::cout << "conjugate(z7) = " << z7_conj << "\n";
    std::cout << "z7 + conjugate(z7) = " << sum << "\n";
    std::cout << "Imaginary part should be 0: " 
              << (std::abs(sum.GetImaginaryPart()) < 1e-10 ? "YES" : "NO") << "\n";
    std::cout << "Real part should be " << 2*z7.GetRealPart() << ": " 
              << (std::abs(sum.GetRealPart() - 2*z7.GetRealPart()) < 1e-10 ? "YES" : "NO") << "\n";
    
    // Test 7: Integration test with all features
    std::cout << "\nTest 7: Integration Test\n";
    std::cout << "------------------------\n";
    ComplexNumber a(1.0, 2.0);
    ComplexNumber b(3.0);  // Using single-argument constructor
    ComplexNumber c = a + b;  // Addition
    std::cout << "a = " << a << "\n";
    std::cout << "b = " << b << " (real number)\n";
    std::cout << "c = a + b = " << c << "\n";
    std::cout << "RealPart(c) = " << RealPart(c) << "\n";
    std::cout << "ImaginaryPart(c) = " << ImaginaryPart(c) << "\n";
    std::cout << "Conjugate(c) = " << c.CalculateConjugate() << "\n";
    
    // Test 8: Copy constructor with operations
    std::cout << "\nTest 8: Copy Constructor in Operations\n";
    std::cout << "---------------------------------------\n";
    ComplexNumber z8(4.0, 5.0);
    ComplexNumber z9 = z8;  // Copy using assignment (calls copy constructor)
    ComplexNumber z10(z8);  // Direct copy construction
    std::cout << "Original z8 = " << z8 << "\n";
    std::cout << "Copy z9 (via assignment) = " << z9 << "\n";
    std::cout << "Copy z10 (via constructor) = " << z10 << "\n";
    
    // Test 9: Feature 6 - SetToConjugate method
    std::cout << "\nTest 9: SetToConjugate Method\n";
    std::cout << "------------------------------\n";
    ComplexNumber z11(7.0, 3.0);
    std::cout << "Before: z11 = " << z11 << "\n";
    z11.SetToConjugate();
    std::cout << "After SetToConjugate(): z11 = " << z11 << "\n";
    std::cout << "Imaginary part changed sign? " 
              << (z11.GetImaginaryPart() == -3.0 ? "YES" : "NO") << "\n";
    
    // Test 10: Feature 7 - Matrix exponential
    std::cout << "\nTest 10: Matrix Exponential (Feature 7)\n";
    std::cout << "----------------------------------------\n";
    const int n = 3;
    
    // Allocate 3x3 matrix
    ComplexNumber** A = AllocateComplexMatrix(n);
    ComplexNumber** expA = AllocateComplexMatrix(n);
    
    // Initialize a simple test matrix
    A[0][0] = ComplexNumber(0.0, 0.0);
    A[0][1] = ComplexNumber(1.0, 0.0);
    A[0][2] = ComplexNumber(0.0, 0.0);
    
    A[1][0] = ComplexNumber(0.0, 0.0);
    A[1][1] = ComplexNumber(0.0, 0.0);
    A[1][2] = ComplexNumber(1.0, 0.0);
    
    A[2][0] = ComplexNumber(0.0, 0.0);
    A[2][1] = ComplexNumber(0.0, 0.0);
    A[2][2] = ComplexNumber(0.0, 0.0);
    
    PrintComplexMatrix(A, n, "Matrix A");
    
    // Compute exp(A)
    MatrixExponential(A, expA, n, 20);
    PrintComplexMatrix(expA, n, "exp(A)");
    
    // IMPORTANT: Deallocate memory (as required in the exercise)
    std::cout << "Deallocating matrix memory...\n";
    DeallocateComplexMatrix(A, n);
    DeallocateComplexMatrix(expA, n);
    std::cout << "Memory successfully deallocated!\n\n";
    
    // Test 11: Feature 8 - Special cases
    std::cout << "\nTest 11: Special Cases (Feature 8)\n";
    std::cout << "-----------------------------------\n";
    
    // Test (0+0i)^n should equal 0 for n > 0
    ComplexNumber zero(0.0, 0.0);
    std::cout << "Testing (0+0i)^n:\n";
    for (int power = 1; power <= 3; power++) {
        ComplexNumber result = zero.CalculatePower(power);
        std::cout << "  (0+0i)^" << power << " = " << result;
        double modulus = result.CalculateModulus();
        std::cout << " (modulus = " << modulus << ", should be ~0)";
        std::cout << (modulus < 1e-10 ? " ✓" : " ✗") << "\n";
    }
    
    // Test any number^0 should equal 1
    std::cout << "\nTesting z^0 should equal 1:\n";
    ComplexNumber testCases[] = {
        ComplexNumber(0.0, 0.0),
        ComplexNumber(5.0, 3.0),
        ComplexNumber(-2.0, 4.0),
        ComplexNumber(1.0, 0.0)
    };
    
    for (int i = 0; i < 4; i++) {
        ComplexNumber result = testCases[i].CalculatePower(0);
        std::cout << "  " << testCases[i] << "^0 = " << result;
        bool isOne = (std::abs(result.GetRealPart() - 1.0) < 1e-10 && 
                      std::abs(result.GetImaginaryPart()) < 1e-10);
        std::cout << (isOne ? " ✓" : " ✗") << "\n";
    }
    
    // Test (1+0i)^n should equal 1
    std::cout << "\nTesting (1+0i)^n should equal 1:\n";
    ComplexNumber one(1.0, 0.0);
    for (int power = 1; power <= 5; power++) {
        ComplexNumber result = one.CalculatePower(power);
        bool isOne = (std::abs(result.GetRealPart() - 1.0) < 1e-10 && 
                      std::abs(result.GetImaginaryPart()) < 1e-10);
        std::cout << "  (1+0i)^" << power << " = " << result;
        std::cout << (isOne ? " ✓" : " ✗") << "\n";
    }
    
    // Test i^4 should equal 1
    std::cout << "\nTesting powers of i:\n";
    ComplexNumber i(0.0, 1.0);
    std::cout << "  i^1 = " << i.CalculatePower(1) << "\n";
    std::cout << "  i^2 = " << i.CalculatePower(2) << " (should be -1+0i)\n";
    std::cout << "  i^3 = " << i.CalculatePower(3) << " (should be 0-1i)\n";
    std::cout << "  i^4 = " << i.CalculatePower(4) << " (should be 1+0i)\n";
    
    // Test conjugate special cases
    std::cout << "\nTesting conjugate special cases:\n";
    ComplexNumber realOnly(5.0, 0.0);
    ComplexNumber imagOnly(0.0, 3.0);
    std::cout << "  Conjugate of real number " << realOnly << " = " 
              << realOnly.CalculateConjugate() << "\n";
    std::cout << "  Conjugate of imaginary number " << imagOnly << " = " 
              << imagOnly.CalculateConjugate() << "\n";
    
    // Test SetToConjugate twice returns original
    std::cout << "\nTesting SetToConjugate twice:\n";
    ComplexNumber z12(4.0, -2.0);
    ComplexNumber original = z12;
    std::cout << "  Original: " << original << "\n";
    z12.SetToConjugate();
    std::cout << "  After 1st conjugate: " << z12 << "\n";
    z12.SetToConjugate();
    std::cout << "  After 2nd conjugate: " << z12 << "\n";
    bool backToOriginal = (std::abs(original.GetRealPart() - z12.GetRealPart()) < 1e-10 && 
                           std::abs(original.GetImaginaryPart() - z12.GetImaginaryPart()) < 1e-10);
    std::cout << "  Back to original? " << (backToOriginal ? "YES ✓" : "NO ✗") << "\n";
    
    std::cout << "\n===========================================\n";
    std::cout << "      All tests completed successfully!     \n";
    std::cout << "===========================================\n";
    
    return 0;
}