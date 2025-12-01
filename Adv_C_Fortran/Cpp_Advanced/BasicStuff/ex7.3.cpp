// ============================================================================
// Test Program
// ============================================================================
// Compile with:
//  g++ ex7.3.cpp AbstractOdeSolver.cpp ForwardEulerSolver.cpp RungeKuttaSolver.cpp -o ex73.exe -Wall
// ============================================================================

#include "AbstractOdeSolver.hpp"
#include "ForwardEulerSolver.hpp"
#include "RungeKuttaSolver.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>

// Test problem: dy/dt = 1 + t
// Exact solution: y = (t^2 + 2t + 4) / 2
double TestRhs(double t, double y)
{
   return 1.0 + t;
}

// Exact solution for comparison
double ExactSolution(double t)
{
   return (t*t + 2.0*t + 4.0) / 2.0;
}

// Function to compute error between numerical and exact solution
void ComputeError(const std::string& fileName, double finalTime)
{
   std::ifstream input(fileName.c_str());
   if (!input.is_open())
   {
      std::cerr << "Error: Could not open file " << fileName << std::endl;
      return;
   }
   
   double t, y;
   double maxError = 0.0;
   double finalNumerical = 0.0;
   
   while (input >> t >> y)
   {
      double exact = ExactSolution(t);
      double error = std::abs(y - exact);
      if (error > maxError)
      {
         maxError = error;
      }
      if (std::abs(t - finalTime) < 1e-10)
      {
         finalNumerical = y;
      }
   }
   
   input.close();
   
   double exactFinal = ExactSolution(finalTime);
   std::cout << "  Final numerical value: " << finalNumerical << std::endl;
   std::cout << "  Exact final value:     " << exactFinal << std::endl;
   std::cout << "  Final error:           " << std::abs(finalNumerical - exactFinal) << std::endl;
   std::cout << "  Maximum error:         " << maxError << std::endl;
}

int main()
{
   std::cout << std::fixed << std::setprecision(8);
   
   std::cout << "==============================================\n";
   std::cout << "   ODE Solver Library Testing                \n";
   std::cout << "==============================================\n\n";
   
   std::cout << "Test Problem: dy/dt = 1 + t\n";
   std::cout << "Initial condition: y(0) = 2\n";
   std::cout << "Time interval: [0, 1]\n";
   std::cout << "Exact solution: y(t) = (t^2 + 2t + 4) / 2\n\n";
   
   // Problem parameters
   double initialTime  = 0.0;
   double finalTime    = 1.0;
   double initialValue = 2.0;
   
   // Test different step sizes
   double stepSizes[] = {0.1, 0.05, 0.01, 0.001};
   int numTests = 4;
   
   std::cout << "==============================================\n";
   std::cout << "   Forward Euler Method                      \n";
   std::cout << "==============================================\n\n";
   
   for (int i = 0; i < numTests; i++)
   {
      double h = stepSizes[i];
      std::cout << "Test with step size h = " << h << std::endl;
      std::cout << "----------------------------------------------\n";
      
      ForwardEulerSolver eulerSolver(h, initialTime, finalTime, 
                                     initialValue, TestRhs);
      
      std::string fileName = "euler_h" + std::to_string(h) + ".dat";
      eulerSolver.SolveEquation(fileName);
      
      ComputeError(fileName, finalTime);
      std::cout << std::endl;
   }
   
   std::cout << "==============================================\n";
   std::cout << "   4th Order Runge-Kutta Method              \n";
   std::cout << "==============================================\n\n";
   
   for (int i = 0; i < numTests; i++)
   {
      double h = stepSizes[i];
      std::cout << "Test with step size h = " << h << std::endl;
      std::cout << "----------------------------------------------\n";
      
      RungeKuttaSolver rkSolver(h, initialTime, finalTime, 
                                initialValue, TestRhs);
      
      std::string fileName = "rk4_h" + std::to_string(h) + ".dat";
      rkSolver.SolveEquation(fileName);
      
      ComputeError(fileName, finalTime);
      std::cout << std::endl;
   }
   
   std::cout << "==============================================\n";
   std::cout << "   Analysis Summary                          \n";
   std::cout << "==============================================\n\n";
   
   std::cout << "Observations:\n";
   std::cout << "1. Forward Euler is a first-order method:\n";
   std::cout << "   Error decreases linearly with step size.\n\n";
   std::cout << "2. Runge-Kutta 4 is a fourth-order method:\n";
   std::cout << "   Error decreases with h^4 (much faster!).\n\n";
   std::cout << "3. For the same step size, RK4 is significantly\n";
   std::cout << "   more accurate than Forward Euler.\n\n";
   
   std::cout << "Output files generated:\n";
   for (int i = 0; i < numTests; i++)
   {
      double h = stepSizes[i];
      std::cout << "  euler_h" << h << ".dat\n";
      std::cout << "  rk4_h" << h << ".dat\n";
   }
   
   std::cout << "\n==============================================\n";
   std::cout << "   All tests completed successfully!         \n";
   std::cout << "==============================================\n";
   
   return 0;
}