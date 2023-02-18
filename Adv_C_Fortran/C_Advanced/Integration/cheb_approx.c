/*
=======================================================
 Compile and execute with:
    $ gcc cheb_approx.c -o cheb -lgsl -lm
 or on Mac:
    $ gcc cheb_approx.c -o cheb -I/opt/homebrew/Cellar/gsl/2.7.1/include -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas
    $ ./cheb
=======================================================
*/
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#define M_PI 3.14159265358979323846
double func (double x, void *params) {
    // Define the function to be integrated here
    return 1/sqrt(2*M_PI)*exp(-x*x/2);  //CDF of Normal Dist.
}

int main() {
    gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc (1000);
    gsl_function F;
    F.function = &func;
    F.params = 0;

    double result, error;
    double a = 0, b = 2; // integration limits
    gsl_integration_cquad (&F, a, b, 0, 1e-8, w, &result, &error, NULL);

    printf ("result = % .18f\n", result);
    printf ("estimated error = % .18f\n", error);
    printf ("intervals = %ld\n", w->size);

    gsl_integration_cquad_workspace_free (w);
    return 0;
}
/*
========================================================
The gsl_integration_cquad function is a function for 
numerical integration that is part of the GNU Scientific 
Library (GSL). It uses a fixed-sample Chebyshev 
integration algorithm, which is a type of adaptive 
integration algorithm.

One of the key features of the gsl_integration_cquad 
function is its ability to automatically adapt the 
number of function evaluations used to achieve a given 
level of accuracy. This is done by recursively 
subdividing the integration interval into smaller 
regions, where the function is evaluated at a set of 
Chebyshev nodes. These nodes are chosen so that they 
are densely packed near the endpoints of the interval, 
where the function is typically more rapidly varying. 
This allows the algorithm to achieve high accuracy with 
fewer function evaluations than a naive fixed-sample 
integration algorithm.

Another advantage of this function is that it is able 
to return an estimate of the error in the integration 
result, which can be used to determine if the result 
is accurate enough for a given application. This is 
important in cases where the result of the integration 
is used as input to another computation, and it is 
crucial to know the error of this input.

Additionally, the gsl_integration_cquad is able to 
handle functions which are not smooth, or which have 
singularities, by using a technique called "integration 
by parts" which reduces the problem to a simpler one.

In summary, gsl_integration_cquad is an efficient and 
robust algorithm for numerical integration that can 
handle a wide range of functions and automatically 
adapts to achieve a desired level of accuracy with 
minimal computational effort.
========================================================
When approximating an integral using Chebyshev nodes,
one commonly used method is to interpolate the 
function being integrated using a polynomial that 
passes through the function's values at the Chebyshev 
nodes. The integral can then be approximated by 
evaluating the polynomial's definite integral over 
the integration interval.

One way to express this approximation in mathematical 
notation is to use the Clenshaw-Curtis quadrature 
formula:

∫ f(x) dx ≈ (b-a)/2 ∑ w_i f(x_i)

Where f(x) is the function being integrated, 
a and b are the integration limits, x_i are 
the Chebyshev nodes, and $w_i$ are the 
corresponding weights.

In this formula, the Chebyshev nodes are 
chosen as the roots of the Chebyshev 
polynomials. The weights are chosen to 
minimize the error in the approximation.

Alternatively, we can use the Gauss-Chebyshev 
quadrature formula, which is similar to the 
Clenshaw-Curtis formula, but the weights are 
chosen differently:

∫ f(x) dx ≈ (b-a)/2 ∑ w_i f( (b+a)/2 + (b-a)/2 cos(iπ/n))
========================================================
*/
