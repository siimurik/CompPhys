/*
 Compile and execute with:
    > gcc euler_ODE.c -o euler
    > ./euler
-----------------------------------------------------------
 To solve an initial value problem for a single first-order
 ordinary differential equation (ODE) in C programming, you 
 can use the Euler method. This method involves approximating 
 the solution to the ODE at a particular point in time by 
 using the slope of the function at that point and the 
 difference in time between the point and the initial point.

 In this example, the function f() represents the ODE that 
 we are trying to solve. The initial time, initial value, 
 step size, and number of steps are all defined at the 
 beginning of the main() function. The Euler method is then 
 used to iteratively approximate the solution to the ODE 
 at each time step. The current time and approximated 
 solution are printed at each step.
-----------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>

// Function representing the ODE
double f(double t, double y) {
    // Return the value of the ODE at the given t and y
    //return y + t;
    return 0.1*y - 2.0E-4*y*y;
}

int main(int argc, char** argv) {
    // Initialize variables
    double t0 = 0.0;  // Initial time
    double y0 = 1.0;  // Initial value
    double h = 0.1;   // Step size
    int n = 5;       // Number of steps

    // Loop through each step
    for (int i = 0; i < n; i++) {
        // Use Euler's method to approximate the solution at the next time step
        double y1 = y0 + h * f(t0, y0);

        // Print the current time and approximated solution
        printf("t = %f, y = %f\n", t0, y0);

        // Update the time and solution for the next iteration
        t0 += h;
        y0 = y1;
  }

    return 0;
}
