// ============================================================================
// RungeKuttaSolver.hpp
// ============================================================================

#ifndef RUNGEKUTTASOLVERHDEF
#define RUNGEKUTTASOLVERHDEF

#include "AbstractOdeSolver.hpp"
#include <string>

class RungeKuttaSolver : public AbstractOdeSolver
{
private:
   // Function pointer for user-specified right-hand side
   double (*mpRhsFunc)(double, double);
   
protected:
   // Implementation of RightHandSide using function pointer
   double RightHandSide(double t, double y);
   
   // Implementation of 4th order Runge-Kutta method
   double ComputeNextY(double t, double y);
   
public:
   // Constructor
   RungeKuttaSolver(double stepSize, double initialTime,
                    double finalTime, double initialValue,
                    double (*rhsFunc)(double, double));
};

#endif