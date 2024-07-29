/*
Compile and execute with:
    $ gcc SOR_param.c -o SOR_param -lm
    $ ./SOR_param
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define SIZE 61

void SOR_steady(double T[SIZE][SIZE], int nx, int ny, double w, double* err) {
    double Told[SIZE][SIZE];
    memcpy(Told, T, sizeof(Told)); // Copy T to Told
    *err = 0.0;

    // Space loop
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            T[i][j] = (1 - w) * Told[i][j] + w * 0.25 * (T[i - 1][j] + Told[i + 1][j] + T[i][j - 1] + Told[i][j + 1]);
        }
    }

    // Compute the error
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double diff = fabs(T[i][j] - Told[i][j]);
            if (diff > *err) {
                *err = diff;
            }
        }
    }
}

int calculate_steps(double start, double end, int precision) {
    return (int)((end - start) * pow(10, precision)) + 1;
}

int main() {
    double L = 1.0;
    int nx = SIZE, ny = SIZE;
    double x[SIZE], y[SIZE], h, tol = 1e-4, err = 1.0;
    int count = 0;
    double T[SIZE][SIZE] = {0};
    double T1[SIZE][SIZE];

    // Define omega range and precision
    double w_start = 1.81;
    double w_end = 1.82;
    int precision = 6;
    int w_steps = calculate_steps(w_start, w_end, precision);
    double w[w_steps];
    int niter[w_steps];

    // Define x and y vectors
    for (int i = 0; i < SIZE; i++) {
        x[i] = i * (L / (nx - 1));
        y[i] = i * (L / (ny - 1));
    }
    h = x[1] - x[0];

    // Initializing the temperature field with boundary conditions
    double T_L = 300.0; // Left boundary temperature
    double T_T = 500.0; // Top boundary temperature
    double T_R = 300.0; // Right boundary temperature
    double T_B = 500.0; // Bottom boundary temperature

    for (int i = 1; i < nx - 1; i++) {
        T[i][0]      = T_L;    // Left boundary
        T[i][nx - 1] = T_R;    // Right boundary
        T[0][i]      = T_T;    // Top boundary
        T[nx - 1][i] = T_B;    // Bottom boundary
    }

    // Calculate the average temperature at the corners
    T[0][0] = (T_T + T_L) / 2.0;          // Top-left corner
    T[0][nx - 1] = (T_T + T_R) / 2.0;     // Top-right corner
    T[ny - 1][0] = (T_L + T_B) / 2.0;     // Bottom-left corner
    T[ny - 1][nx - 1] = (T_R + T_B) / 2.0; // Bottom-right corner

    memcpy(T1, T, sizeof(T1)); // Copy T to T1

    for (int i = 0; i < w_steps; i++) {
        w[i] = w_start + i * pow(10, -precision);
        niter[i] = 0;
    }

    // For steady state
    for (int k = 0; k < w_steps; k++) {
        memcpy(T, T1, sizeof(T)); // Reset T to T1
        err = 1.0;
        count = 0;

        while (err >= tol) {
            SOR_steady(T, nx, ny, w[k], &err);
            count++;
        }

        niter[k] = count;
    }

    // Optimal parameter
    double min_niter = niter[0];
    double w_opt = w[0];
    for (int i = 1; i < w_steps; i++) {
        if (niter[i] < min_niter) {
            min_niter = niter[i];
            w_opt = w[i];
        }
    }

    // Output the results
    printf("Optimal omega (SOR parameter): %f\n", w_opt);
    printf("Number of iterations for optimal omega: %d\n", (int)min_niter);

    // Optionally: write results to a file for plotting
    FILE *fp = fopen("results.csv", "w");
    if (fp) {
        fprintf(fp, "omega, iterations\n");
        for (int i = 0; i < w_steps; i++) {
            fprintf(fp, "%f, %d\n", w[i], niter[i]);
        }
        fclose(fp);
    } else {
        perror("Error opening results file");
    }

    return 0;
}
