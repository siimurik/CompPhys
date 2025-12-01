// ============================================================================
// ForwardEulerSolver.hpp
// ============================================================================

#ifndef FORWARDEULERSOLVERHDEF
#define FORWARDEULERSOLVERHDEF

#include "AbstractOdeSolver.hpp"
#include <string>

class ForwardEulerSolver : public AbstractOdeSolver
{
private:
   // Function pointer for user-specified right-hand side
   double (*mpRhsFunc)(double, double);
   
protected:
   // Implementation of RightHandSide using function pointer
   double RightHandSide(double t, double y);
   
   // Implementation of Forward Euler method
   double ComputeNextY(double t, double y);
   
public:
   // Constructor
   ForwardEulerSolver(double stepSize, double initialTime,
                      double finalTime, double initialValue,
                      double (*rhsFunc)(double, double));
};

#endif