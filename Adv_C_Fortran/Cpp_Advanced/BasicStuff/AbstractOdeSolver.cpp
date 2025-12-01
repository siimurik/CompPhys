// ============================================================================
// AbstractOdeSolver.cpp
// ============================================================================
#include "AbstractOdeSolver.hpp"
#include <fstream>
#include <cassert>
#include <iostream>

// Constructor implementation
AbstractOdeSolver::AbstractOdeSolver(double stepSize, double initialTime,
                                     double finalTime, double initialValue)
{
   mStepSize     = stepSize;
   mInitialTime  = initialTime;
   mFinalTime    = finalTime;
   mInitialValue = initialValue;
   
   // Verify that inputs make sense
   assert(mStepSize > 0.0);
   assert(mFinalTime > mInitialTime);
}

// Method to solve the ODE and write results to file
void AbstractOdeSolver::SolveEquation(const std::string& fileName)
{
   // Open output file
   std::ofstream output(fileName.c_str());
   if (!output.is_open())
   {
      std::cerr << "Error: Could not open file " << fileName << std::endl;
      return;
   }
   
   // Calculate number of steps
   int numSteps = static_cast<int>((mFinalTime - mInitialTime) / mStepSize);
   
   // Initialize
   double t = mInitialTime;
   double y = mInitialValue;
   
   // Write initial condition
   output << t << " " << y << "\n";
   
   // Solve using the specific method implemented in derived class
   for (int i = 0; i < numSteps; i++)
   {
      y = ComputeNextY(t, y);
      t += mStepSize;
      output << t << " " << y << "\n";
   }
   
   output.close();
   std::cout << "Solution written to " << fileName << std::endl;
}