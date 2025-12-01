// ============================================================================
// AbstractOdeSolver.hpp
// ============================================================================
#ifndef ABSTRACTODESOLVERHEADERDEF
#define ABSTRACTODESOLVERHEADERDEF

#include <string>

class AbstractOdeSolver
{
private:
   // Don't allow default instantiation
   // (copy assignment operator and copy constructor)
   AbstractOdeSolver(const AbstractOdeSolver& otherOdeSolver) {}
   void operator=(const AbstractOdeSolver& otherOdeSolver) {}
   
protected:
   double mStepSize;                    // Step size h
   double mInitialTime;                 // Initial time T0
   double mFinalTime;                   // Final time T1
   double mInitialValue;                // Initial value Y0
   
   // Pure virtual function for right-hand side f(t,y)
   virtual double RightHandSide(double t, double y) = 0;
   
   // Pure virtual function for the specific numerical method
   virtual double ComputeNextY(double t, double y) = 0;
   
public:
   // Constructor
   AbstractOdeSolver(double stepSize, double initialTime, 
                     double finalTime, double initialValue);
   
   // Method to solve the equation and write to file
   void SolveEquation(const std::string& fileName);
   
   // Destructor (virtual to allow proper cleanup in derived classes)
   virtual ~AbstractOdeSolver() {}
};

#endif