#include <stdio.h>
#include <stdlib.h>

// Gaussian elimination algorithm
// Given a system of linear equations in the form Ax = b,
// where A is an n x n matrix, b is a vector of length n,
// and x is the vector of unknowns, this function will
// solve for x and return the solution in x.
void gaussian_elimination(double A[][3], double *b, double *x, int n)
{
    int i, j, k;

    // Eliminate variables to obtain a reduced row-echelon form
    for (k = 0; k < n; k++)
    {
        // Find the pivot row (i.e. the row with the largest coefficient
        // for the current variable)
        int i_max = k;
        for (i = k + 1; i < n; i++)
        {
            if (fabs(A[i][k]) > fabs(A[i_max][k]))
            {
                i_max = i;
            }
        }

        // Swap the pivot row with the current row
        if (i_max != k)
        {
            for (j = k; j < n; j++)
            {
                double temp = A[k][j];
                A[k][j] = A[i_max][j];
                A[i_max][j] = temp;
            }
            double temp = b[k];
            b[k] = b[i_max];
            b[i_max] = temp;
        }

        // Perform row operations to eliminate the current variable
        for (i = k + 1; i < n; i++)
        {
            double factor = A[i][k] / A[k][k];
            for (j = k; j < n; j++)
            {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    // Back-substitute to solve for the unknowns
    for (i = n - 1; i >= 0; i--)
    {
        double sum = 0.0;
        for (j = i + 1; j < n; j++)
        {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }
}

int main(void)
{
    // Define a system of linear equations
    // 3x + 2y - z = 2
    // 2x - 2y + 4z = -3
    // -x + y/2 - z = z/2
    double A[3][3] = {{3, 2, -1}, {2, -2, 4}, {-1, 0.5, -1}};
    double b[3] = {2, -3, 0.5};
    double x[3];

    // Solve the system of linear equations
    gaussian_elimination(A, b, x, 3);

    // Print the solution
    printf("x = %lf\n", x[0]);
    printf("y = %lf\n", x[1]);
    printf("z = %lf\n", x[2]);

    return 0;
}
