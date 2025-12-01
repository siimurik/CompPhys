// ============================================================================
// ForwardEulerSolver.cpp
// ============================================================================

#include "ForwardEulerSolver.hpp"

// Constructor
ForwardEulerSolver::ForwardEulerSolver(double stepSize, double initialTime,
                                       double finalTime, double initialValue,
                                       double (*rhsFunc)(double, double))
   : AbstractOdeSolver(stepSize, initialTime, finalTime, initialValue)
{
   mpRhsFunc = rhsFunc;
}

// Implementation of RightHandSide
double ForwardEulerSolver::RightHandSide(double t, double y)
{
   return (*mpRhsFunc)(t, y);
}

// Implementation of Forward Euler method: y_i = y_{i-1} + h*f(t_{i-1}, y_{i-1})
double ForwardEulerSolver::ComputeNextY(double t, double y)
{
   return y + mStepSize * RightHandSide(t, y);
}