/*===================================================
 Compile and execute with:
    $ gcc integral_test.c -o integral_test -lgsl -lm
    $ ./integral_test
===================================================*/
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <math.h>

double f(double x, void* params) {
    //return 1.0/log(x);
    return 1.0/tgamma(x+1); // tgamma(x+1) = x!
    // return gsl_pow_int(2, x) / tgamma(x+1.0); // should actually converge
}

int main() {
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    double a = 4.0;
    //double b = INFINITY;
    double epsabs = 0.0, epsrel = 1e-8;
    double result, error;
    gsl_function F;
    F.function = &f;
    F.params = 0;

    // Adding my own custom error announcement, if the series is 
    // divergent. Otherwise it just says that the gsl solver ran
    // into an error like:
    // "gsl: qags.c:543: ERROR: number of iterations was insufficient"
    gsl_set_error_handler_off();

    // Integral solver with infinite upper bound
    int status = gsl_integration_qagiu(&F, a, epsabs, epsrel, 1000, w, &result, &error); // for a to +oo

    // If successful, show result
    if (status == EXIT_SUCCESS) {
        printf("∫ f(x) dx = % .18f\n", result);
        printf("Estimated error = %g\n", error);
        /*
        double sum = a, term = 1.0, tol = epsrel;
        int n = 0;
        while (fabs(term) > tol) {
            sum += term;
            n++;
            term = 1.0 / tgamma(n+1); // compute the next term
            //term = pow(2,n)/tgamma(n+1);
        }
        printf("Σ f(x_n) = %.16f\n", sum);
        printf("n = %i\n", n);
        */
    }
    else { // If divergent, let the user know.
        printf("GSL integration failed. The series is divergent.\n");
        printf("*In Uncle Rodger voice*\n You fucked uppp...\n");
    }
    gsl_integration_workspace_free(w);
    return 0;
}
