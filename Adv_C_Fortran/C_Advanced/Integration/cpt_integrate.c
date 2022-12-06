/*
 Compile and execute with:
    > gcc -lgsl -lgslcblas -lm cpt_integrate.c -o int
    > ./int
------------------------------------------------------------------------
 To solve an integral numerically in C programming, you can use the GNU 
 Scientific Library (GSL) library. This library provides a collection of
  numerical routines for scientific computing, including routines for 
  numerical integration.

 To use the GSL library, you will first need to install it. Here is an 
 example of how to do this on a Linux system:
    > $ sudo apt-get install libgsl-dev

 This code will solve the integral of x^2 from 0 to 1 numerically, 
 and print the result to the console. The gsl_integration_qag 
 function from the GSL library is used to perform the integration. 
 This function uses a Gaussian quadrature algorithm to compute 
 the integral with a specified error tolerance. 
------------------------------------------------------------------------
*/
#include <stdio.h>
#include <gsl/gsl_integration.h>

// Function to be integrated
double f(double x, void *params)
{
    return x*x;
}

int main()
{
    // Define the limits of the integration
    double a = 0.0;
    double b = 1.0;

    // Define the integration workspace
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    // Define the function to be integrated
    gsl_function F;
    F.function = &f;
    F.params = NULL;

    // Define the error tolerance
    double epsabs = 1e-8;
    double epsrel = 1e-8;

    // Variable to store the result of the integration
    double result;
    double error;

    // Perform the integration
    int status = gsl_integration_qag(&F, a, b, epsabs, epsrel, 1000,
                                     GSL_INTEG_GAUSS61, w, &result, &error);
    if (status != EXIT_SUCCESS) {
        fprintf(stderr, "GSL integration failed with error %d\n", status);
        return 1;
    }

    // Print the result of the integration
    printf("Integral = %.18f\n", result);

    // Free the integration workspace
    gsl_integration_workspace_free(w);

    return 0;
}
