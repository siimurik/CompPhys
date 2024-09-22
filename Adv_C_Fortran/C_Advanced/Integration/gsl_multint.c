#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

// Typedefs for function pointers and parameters
typedef double (*integrand_func)(double, void *);

// Global error tolerance variables
const double ABS_ERROR = 1e-14;
const double REL_ERROR = 1e-14;  // Set relative tolerance here

// Struct to hold parameters for z-integration
typedef struct {
    double x;
    double y;
} xy_params;

// Function to integrate over z
double func_z(double z, void *params) {
    xy_params *p = (xy_params *)params;
    double x = p->x;
    double y = p->y;
    return x * x + y * y - z * z;
}

// Wrapper function to perform integration over z
double integrate_z(double x, double y, integrand_func fz) {
    xy_params p = {x, y};
    double z_min = -sqrt(1 - x * x - y * y);
    double z_max = sqrt(1 - x * x - y * y);

    // Perform integration over z using GSL
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);  // Increased workspace
    gsl_function F;
    F.function = fz;
    F.params = &p;

    double result, error;
    int status = gsl_integration_qags(&F, z_min, z_max, ABS_ERROR, REL_ERROR, 10000, w, &result, &error);
    gsl_integration_workspace_free(w);

    if (status != GSL_SUCCESS) {
        fprintf(stderr, "Warning: Integration over z failed at x=%.5f, y=%.5f (status: %d)\n", x, y, status);
        return 0;  // Return zero for failed integration
    }

    return result;
}

// Function to integrate over y
double integrate_y(double y, void *params) {
    double x = *(double *)params;
    return integrate_z(x, y, func_z);
}

// Wrapper function to perform integration over y
double integrate_y_wrapper(double x, integrand_func fy) {
    double y_min = -sqrt(1 - x * x);
    double y_max = sqrt(1 - x * x);

    // Perform integration over y using GSL
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);  // Increased workspace
    gsl_function F;
    F.function = fy;
    F.params = &x;

    double result, error;
    int status = gsl_integration_qags(&F, y_min, y_max, ABS_ERROR, REL_ERROR, 10000, w, &result, &error);
    gsl_integration_workspace_free(w);

    if (status != GSL_SUCCESS) {
        fprintf(stderr, "Warning: Integration over y failed at x=%.5f (status: %d)\n", x, status);
        return 0;  // Return zero for failed integration
    }

    return result;
}

// Function to integrate over x
double integrate_x(integrand_func fx, double x_min, double x_max) {
    // Perform integration over x using GSL
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);  // Increased workspace
    gsl_function F;
    F.function = fx;

    double result, error;
    int status = gsl_integration_qags(&F, x_min, x_max, ABS_ERROR, REL_ERROR, 10000, w, &result, &error);
    gsl_integration_workspace_free(w);

    if (status != GSL_SUCCESS) {
        fprintf(stderr, "Warning: Integration over x failed (status: %d)\n", status);
        return 0;  // Return zero for failed integration
    }

    return result;
}

// Function that wraps the x-integration, called by GSL
double integrate_x_wrapper(double x, void *params) {
    return integrate_y_wrapper(x, integrate_y);
}

int main() {
    // Turn off the GSL default error handler
    gsl_set_error_handler_off();

    // Define the limits of integration for x
    double x_min = -1.0;
    double x_max = 1.0;

    // Perform the triple integration over x, y, and z
    double result = integrate_x(integrate_x_wrapper, x_min, x_max);

    // Print the result with high precision
    printf("Triple Integral Result: %.15f\n", result);

    return 0;
}
