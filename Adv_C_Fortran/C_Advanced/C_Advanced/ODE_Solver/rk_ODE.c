/*
 Compile and execute with:
    > gcc rk_ODE.c -o rk
    > ./rk
--------------------------------------------------------------
 To solve an initial value problem for a single first-order 
 ordinary differential equation (ODE) using the Runge-Kutta
 method in C programming, you can use the fourth-order 
 Runge-Kutta (RK4) method. This method involves using the 
 slope of the function at four different points to approximate 
 the solution to the ODE at the next point in time.

 In this example, the function f() represents the ODE that we 
 are trying to solve. The initial time, initial value, step size, 
 and number of steps are all defined at the beginning of the 
 main() function. The RK4 method is then used to iteratively 
 approximate the solution to the ODE at each time step. The 
 current time and approximated solution are printed at each step.
--------------------------------------------------------------
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
    double t0 = 0.0;    // Initial time
    double y0 = 1.0;    // Initial value
    double h = 0.01;    // Step size
    int n = 2001;       // Number of steps

    // Loop through each step
    for (int i = 0; i < n; i++) {
        // Use RK4 to approximate the solution at the next time step
        double k1 = h * f(t0, y0);
        double k2 = h * f(t0 + 0.5 * h, y0 + 0.5 * k1);
        double k3 = h * f(t0 + 0.5 * h, y0 + 0.5 * k2);
        double k4 = h * f(t0 + h, y0 + k3);
        double y1 = y0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;

        // Print the current time and approximated solution
        printf("t = %f, y = %f\n", t0, y0);

        // Update the time and solution for the next iteration
        t0 += h;
        y0 = y1;
    }

    return 0;
}
