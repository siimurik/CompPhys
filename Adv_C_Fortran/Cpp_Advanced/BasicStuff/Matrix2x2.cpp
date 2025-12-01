// ============================================================================
// Matrix2x2.cpp
// ===========================================================================
#include "Matrix2x2.hpp"

// Feature 1: Default constructor - initializes all entries to zero
Matrix2x2::Matrix2x2()
{
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            mData[i][j] = 0.0;
        }
    }
}

// Feature 2: Copy constructor
Matrix2x2::Matrix2x2(const Matrix2x2& other)
{
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            mData[i][j] = other.mData[i][j];
        }
    }
}

// Feature 3: Constructor that specifies all four entries
// Matrix layout: [a00  a01]
//                [a10  a11]
Matrix2x2::Matrix2x2(double a00, double a01, double a10, double a11)
{
    mData[0][0] = a00;
    mData[0][1] = a01;
    mData[1][0] = a10;
    mData[1][1] = a11;
}

// Feature 4: Calculate determinant
// For 2x2 matrix: det = a00*a11 - a01*a10
double Matrix2x2::CalculateDeterminant() const
{
    return mData[0][0]*mData[1][1] - mData[0][1]*mData[1][0];
}

// Feature 5: Calculate inverse
// Inverse of [a b] = 1/(ad-bc) * [ d -b]
//            [c d]               [-c  a]
Matrix2x2 Matrix2x2::CalculateInverse() const
{
    double det = CalculateDeterminant();
    if (det == 0.0) {
        std::cerr << "Warning: Matrix is singular (det = 0), cannot invert.\n";
        // Return zero matrix if singular
        return Matrix2x2();
    }
    double invDet = 1.0 / det;
    return Matrix2x2(
         invDet * mData[1][1],
        -invDet * mData[0][1],
        -invDet * mData[1][0],
         invDet * mData[0][0]
    );
}

// Feature 6: Assignment operator
Matrix2x2& Matrix2x2::operator=(const Matrix2x2& other)
{
    // Check for self-assignment
    if (this != &other) {
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                mData[i][j] = other.mData[i][j];
            }
        }
    }
    return *this;
}

// Feature 7: Unary minus operator
Matrix2x2 Matrix2x2::operator-() const
{
    Matrix2x2 result;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            result.mData[i][j] = -mData[i][j];
        }
    }
    return result;
}

// Feature 8: Binary addition operator
Matrix2x2 Matrix2x2::operator+(const Matrix2x2& other) const
{
    Matrix2x2 result;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            result.mData[i][j] = mData[i][j] + other.mData[i][j];
        }
    }
    return result;
}

// Feature 8: Binary subtraction operator
Matrix2x2 Matrix2x2::operator-(const Matrix2x2& other) const
{
    Matrix2x2 result;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            result.mData[i][j] = mData[i][j] - other.mData[i][j];
        }
    }
    return result;
}

// Feature 9: Multiply matrix by a scalar
Matrix2x2 Matrix2x2::MultScalar(double scalar) const
{
    Matrix2x2 result;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            result.mData[i][j] = mData[i][j] * scalar;
        }
    }
    return result;
}

// Helper: Get element at position (i, j)
double Matrix2x2::GetElement(int i, int j) const
{
    if (i < 0 || i > 1 || j < 0 || j > 1) {
        std::cerr << "Error: Index out of bounds.\n";
        return 0.0;
    }
    return mData[i][j];
}

// Helper: Set element at position (i, j)
void Matrix2x2::SetElement(int i, int j, double value)
{
    if (i < 0 || i > 1 || j < 0 || j > 1) {
        std::cerr << "Error: Index out of bounds.\n";
        return;
    }
    mData[i][j] = value;
}

// Output operator overload
std::ostream& operator<<(std::ostream& output, const Matrix2x2& m)
{
    output << "[ " << m.mData[0][0] << "  " << m.mData[0][1] << " ]\n";
    output << "[ " << m.mData[1][0] << "  " << m.mData[1][1] << " ]\n";
    return output;
}

