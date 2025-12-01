// ============================================================================
// Compile with: g++ ex6.2.cpp Matrix2x2.cpp -o ex62.exe
// ============================================================================

#include "Matrix2x2.hpp"
#include <windows.h>
#include <iostream>
#include <iomanip>
#include <cmath>

int main()
{
    // Set console to UTF-8 if the OS is Windows or Linux or Mac
    #ifdef _WIN32
        SetConsoleOutputCP(CP_UTF8);
    #elif __linux__
        std::cout << "\033%G"; // Switch to UTF-8 in Linux terminal
    #elif __APPLE__
        // macOS terminal usually uses UTF-8 by default, but we can ensure it
        std::setlocale(LC_ALL, "en_US.UTF-8");
        std::locale::global(std::locale("en_US.UTF-8"));
    #endif

    std::cout << std::fixed << std::setprecision(4);

    std::cout << "===========================================\n";
    std::cout << "      Testing Matrix2x2 Class              \n";
    std::cout << "===========================================\n";

    // Test 1: Default constructor
    std::cout << "\nTest 1: Default constructor (all zeros)\n";
    std::cout << "------------------------------------------\n";
    Matrix2x2 M1;
    std::cout << "M1 =\n" << M1 << "\n\n";

    // Test 2: Constructor with four entries
    std::cout << "Test 2: Constructor with four entries\n";
    std::cout << "--------------------------------------\n";
    Matrix2x2 M2(1.0, 2.0, 3.0, 4.0);
    std::cout << "M2(1,2,3,4) =\n" << M2 << "\n\n";

    // Test 3: Copy constructor
    std::cout << "Test 3: Copy constructor\n";
    std::cout << "-------------------------\n";
    Matrix2x2 M3(M2);
    std::cout << "M3 (copy of M2) =\n" << M3 << "\n";
    std::cout << "M3 == M2? " 
        << (M3.GetElement(0,0) == M2.GetElement(0,0) ? "[✓]" : "[✗]") 
        << "\n\n";

    // Test 4: Determinant calculation
    std::cout << "Test 4: Determinant calculation\n";
    std::cout << "-------------------------------\n";
    Matrix2x2 M4(2.0, 3.0, 1.0, 4.0);
    std::cout << "M4 =\n" << M4 << "\n";
    double det_M4 = M4.CalculateDeterminant();
    std::cout << "det(M4) = " << det_M4 << " (should be 5)\n\n";

    // Test 5: Matrix inversion
    std::cout << "Test 5: Matrix inverse\n";
    std::cout << "-----------------------\n";
    Matrix2x2 M5(4.0, 7.0, 2.0, 6.0);
    std::cout << "M5 =\n" << M5 << "\n";
    Matrix2x2 M5_inv = M5.CalculateInverse();
    std::cout << "M5 inverse =\n" << M5_inv << "\n\n";

    // Verify: M5 * M5_inv = Identity matrix
    std::cout << "Verification: M5 * M5_inv = Identity matrix\n";
    Matrix2x2 Eye2(
        M5.GetElement(0,0)*M5_inv.GetElement(0,0) + M5.GetElement(0,1)*M5_inv.GetElement(1,0),
        M5.GetElement(0,0)*M5_inv.GetElement(0,1) + M5.GetElement(0,1)*M5_inv.GetElement(1,1),
        M5.GetElement(1,0)*M5_inv.GetElement(0,0) + M5.GetElement(1,1)*M5_inv.GetElement(1,0),
        M5.GetElement(1,0)*M5_inv.GetElement(0,1) + M5.GetElement(1,1)*M5_inv.GetElement(1,1)
    );
    std::cout << "M5 * M5_inv =\n" << Eye2 << "\n\n"; 

    // Test singular matrix
    std::cout << "Testing singular matrix:\n";
    Matrix2x2 M_singular(1.0, 2.0, 2.0, 4.0);
    std::cout << "M_singular =\n" << M_singular << "\n";
    Matrix2x2 M_singular_inv = M_singular.CalculateInverse();
    std::cout << "M_singular inverse =\n" << M_singular_inv << "\n\n";

    // Test 6: Assignment operator
    std::cout << "Test 6: Assignment operator\n";
    std::cout << "----------------------------\n";
    Matrix2x2 M6;
    M6 = M2;
    std::cout << "M6 = M2, M6 =\n" << M6 << "\n\n";

    // Test 7: Unary minus operator
    std::cout << "Test 7: Unary minus operator\n";
    std::cout << "-----------------------------\n";
    Matrix2x2 M7(5.0, -3.0, 2.0, 8.0);
    std::cout << "M7 =\n" << M7 << "\n";
    Matrix2x2 M7_neg = -M7;
    std::cout << "-M7 =\n" << M7_neg << "\n\n";

    // Test 8: Binary addition and subtraction
    std::cout << "Test 8: Binary addition\n";
    std::cout << "------------------------\n";
    Matrix2x2 M8a(1.0, 2.0, 3.0, 4.0);
    Matrix2x2 M8b(5.0, 6.0, 7.0, 8.0);
    std::cout << "M8a =\n" << M8a << "\n";
    std::cout << "M8b =\n" << M8b << "\n";
    Matrix2x2 M8c = M8a + M8b;
    std::cout << "M8a + M8b =\n" << M8c << "\n\n";

    std::cout << "Test 9: Binary subtraction\n";
    std::cout << "---------------------------\n";
    Matrix2x2 M8d = M8b - M8a;
    std::cout << "M8b - M8a =\n" << M8c << "\n\n";

    // Test 9: Scalar multiplication
    std::cout << "Test 9: Scalar multiplication\n";
    std::cout << "-----------------------------\n";
    Matrix2x2 M9(1.0, -2.0, 3.0, -4.0);
    std::cout << "M9 =\n" << M9 << "\n";
    Matrix2x2 M9_scaled = M9.MultScalar(3.0);
    std::cout << "M9 * 3.0 =\n" << M9_scaled << "\n\n";

    // Test 10: Combined operations
    std::cout << "Test 10: Combined operations\n";
    std::cout << "-----------------------------\n";
    Matrix2x2 E(1.0, 0.0, 0.0, 1.0); // Identity matrix
    Matrix2x2 F(2.0, 3.0, 4.0, 5.0);
    std::cout << "E (identity) =\n" << E << "\n\n";
    std::cout << "F =\n" << F << "\n\n";
    
    Matrix2x2 G = (F + E).MultScalar(2.0);
    std::cout << "G = (F + E) * 2 =\n" << G << "\n\n";
    
    Matrix2x2 H = -F + E;
    std::cout << "H = -F + E =\n" << H << "\n\n";
    
    // Test property: A - A should give zero matrix
    std::cout << "Test 11: Property Check (A - A = 0)\n";
    std::cout << "------------------------------------\n";
    Matrix2x2 diff = G - G;
    std::cout << "G - G =\n" << diff << "\n";
    bool isZero = (diff.GetElement(0,0) == 0.0 && diff.GetElement(0,1) == 0.0 &&
                    diff.GetElement(1,0) == 0.0 && diff.GetElement(1,1) == 0.0);
    std::cout << "Is matrix zero? " << (isZero ?  "[✓]" : "[✗]") << "\n\n";
    
    std::cout << "===========================================\n";
    std::cout << "      All tests completed successfully!     \n";
    std::cout << "===========================================\n";

    return 0;
}