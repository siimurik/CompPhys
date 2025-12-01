// ============================================================================
// Matrix2x2.hpp
// ============================================================================
#ifndef MATRIX2X2HEADERDEF
#define MATRIX2X2HEADERDEF

#include <iostream>

class Matrix2x2
{
private:
    double mData[2][2]; // 2x2 array to store matrix entries

public:
    // Feature 1: Default constructor -  initialize all entries to zero
    Matrix2x2();

    // Feature 2: Copy constructor
    Matrix2x2(const Matrix2x2& other);

    // Feature 3: Constructor with four entries
    Matrix2x2(double a00, double a01, double a10, double a11);

    // Feature 4: Method to calculate determinant
    double CalculateDeterminant() const;

    // Feature 5: Method to calculate inverse (returns false if singular)
    Matrix2x2 CalculateInverse() const;

    // Feature 6: Assignment operator overload
    Matrix2x2& operator=(const Matrix2x2& other);

    // Feature 7: Unary minus operator overload
    Matrix2x2 operator-() const;

    // Feature 8: Binary addition and subtraction operator overloads
    Matrix2x2 operator+(const Matrix2x2& other) const;
    Matrix2x2 operator-(const Matrix2x2& other) const;

    // Feature 9: Method to multiply matrix by a scalar
    Matrix2x2 MultScalar(double scalar) const;

    // Helper methods for accessing elements
    double GetElement(int i, int j) const;
    void SetElement(int i, int j, double value);

    // Friend function for output
    friend std::ostream& operator<<(std::ostream& output, const Matrix2x2& m);

};

#endif