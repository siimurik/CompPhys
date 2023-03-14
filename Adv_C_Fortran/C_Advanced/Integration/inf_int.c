/*======================================
 Compile and execute with
    > gcc int.c -o int -lgsl -lm
    > ./int
======================================*/
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

double f(double x, void *params)
{
    return exp(-x*x);   // int_0^+oo exp(-x^2) dx = sqrt(pi)/2 = 0.886227
    //return 3.0*x/((x*x-2.0)*log(x));
    //return 1/((x+1)*(3*x-2));
}

int main()
{
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);

    double result, error;
    double a = 0.0; // lower limit
    double epsabs = 0.0; // absolute error tolerance
    double epsrel = 1e-10; // relative error tolerance

    gsl_function F;
    F.function = &f;
    F.params = NULL;

    //gsl_integration_qagi(&F, epsabs, epsrel, 1000, workspace, &result, &error);   // for -oo to +oo
    gsl_integration_qagiu(&F, a, epsabs, epsrel, 1000, workspace, &result, &error); // for  a  to +oo

    printf("result = %g\n", result);
    printf("error = %g\n", error);

    gsl_integration_workspace_free(workspace);

    return 0;
}
