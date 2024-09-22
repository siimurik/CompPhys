#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_integration.h>
#include <string.h>
#include <math.h>

// Typedef for function pointer
typedef double (*func_ptr)(double, void*);

// Function declarations for user-defined functions
double func_1(double x, void *params);
double func_2(double x, void *params);

// Function to choose the correct function based on user input
func_ptr choose_function(const char *func_name);

int main(void) {
    // Variable declarations
    char func_name[50];
    double lower_bound, upper_bound;
    func_ptr selected_func;
    
    // Integration workspace and results
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function F;
    
    // Get user input
    printf("Enter function (func_1 or func_2): ");
    scanf("%49s", func_name);
    printf("Lower bound: ");
    scanf("%lf", &lower_bound);
    printf("Upper bound: ");
    scanf("%lf", &upper_bound);

    // Choose the function
    selected_func = choose_function(func_name);
    if (selected_func == NULL) {
        fprintf(stderr, "Error: Function not recognized.\n");
        gsl_integration_workspace_free(workspace);
        return EXIT_FAILURE;
    }

    // Set up the GSL function
    F.function = selected_func;
    F.params = NULL;

    // Perform the integration
    gsl_integration_qags(&F, lower_bound, upper_bound, 0, 1e-10, 1000, workspace, &result, &error);

    // Output the result
    printf("Result: %.10f\n", result);
    printf("Estimated error: %.6e\n", error);

    // Clean up
    gsl_integration_workspace_free(workspace);

    return EXIT_SUCCESS;
}

// Sample functions
double func_1(double x, void *params) {
    return sin(x);
}

double func_2(double x, void *params) {
    return exp(-x * x);
}

// Function chooser
func_ptr choose_function(const char *func_name) {
    if (strcmp(func_name, "func_1") == 0) {
        return &func_1;
    } else if (strcmp(func_name, "func_2") == 0) {
        return &func_2;
    } else {
        return NULL;
    }
}
