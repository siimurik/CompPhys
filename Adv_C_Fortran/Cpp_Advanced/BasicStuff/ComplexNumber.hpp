
// ============================================================================
// ComplexNumber.hpp - UPDATED VERSION
// Save this as ComplexNumber.hpp (replaces your existing file)
// ============================================================================
#ifndef COMPLEXNUMBERHEADERDEF
#define COMPLEXNUMBERHEADERDEF

#include <iostream>

class ComplexNumber
{
private:
   double mRealPart;
   double mImaginaryPart;
public:
   ComplexNumber();
   ComplexNumber(double x, double y);
   ComplexNumber(double x);  // Feature 4: Constructor for real numbers
   ComplexNumber(const ComplexNumber& z);  // Feature 3: Copy constructor
   
   double CalculateModulus() const;
   double CalculateArgument() const;
   ComplexNumber CalculatePower(double n) const;
   ComplexNumber CalculateConjugate() const;  // Feature 5: Conjugate method
   
   // Feature 1: Getter methods
   double GetRealPart() const;
   double GetImaginaryPart() const;
   
   ComplexNumber& operator=(const ComplexNumber& z);
   ComplexNumber operator-() const;
   ComplexNumber operator+(const ComplexNumber& z) const;
   ComplexNumber operator-(const ComplexNumber& z) const;
   
   friend std::ostream& operator<<(std::ostream& output, 
                                   const ComplexNumber& z);
   
   // Feature 2: Friend functions for accessing parts
   friend double RealPart(const ComplexNumber& z);
   friend double ImaginaryPart(const ComplexNumber& z);

   // Add this with the other operator overloads
   ComplexNumber operator*(const ComplexNumber& z) const;

   // Add this with the other methods
   void SetToConjugate();
};

#endif