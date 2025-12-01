// ============================================================================
// RungeKuttaSolver.cpp
// ============================================================================

#include "RungeKuttaSolver.hpp"

// Constructor
RungeKuttaSolver::RungeKuttaSolver(double stepSize, double initialTime,
                                   double finalTime, double initialValue,
                                   double (*rhsFunc)(double, double))
   : AbstractOdeSolver(stepSize, initialTime, finalTime, initialValue)
{
   mpRhsFunc = rhsFunc;
}

// Implementation of RightHandSide
double RungeKuttaSolver::RightHandSide(double t, double y)
{
   return (*mpRhsFunc)(t, y);
}

// Implementation of 4th order Runge-Kutta method
double RungeKuttaSolver::ComputeNextY(double t, double y)
{
   double k1 = mStepSize * RightHandSide(t, y);
   double k2 = mStepSize * RightHandSide(t + 0.5*mStepSize, y + 0.5*k1);
   double k3 = mStepSize * RightHandSide(t + 0.5*mStepSize, y + 0.5*k2);
   double k4 = mStepSize * RightHandSide(t + mStepSize, y + k3);
   
   return y + (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
}