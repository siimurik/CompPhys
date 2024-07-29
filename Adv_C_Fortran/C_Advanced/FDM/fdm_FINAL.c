/*
Compile and execute with:
    $ gcc fdm_FINAL.c -o fdm
    $ ./fdm
*/
// Original article: 
// https://skill-lync.com/student-projects/week-5-mid-term-project-solving-the-steady-and-unsteady-2d-heat-conduction-problem-35
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define SIZE 61

// Function prototypes
void matrix_csv(const char *filename, double a[][SIZE], int n, int m);
void vector_csv(const char *filename, double v[SIZE], double u[SIZE], int n);
void initialize_domain(double *Lx, double *Ly, int *nx, int *ny, double *hx, double *hy);
void initialize_vectors(double x[SIZE], double y[SIZE], double hx, double hy, int nx);
void initialize_temperature(double Te[SIZE][SIZE], int nx, int ny);
double jacobi_method(double Te[SIZE][SIZE], double Told[SIZE][SIZE], int nx, int ny, int *iterations);
double gauss_seidel_method(double Te[SIZE][SIZE], double Told[SIZE][SIZE], int nx, int ny, int *iterations);
double sor_method(double Te[SIZE][SIZE], double Told[SIZE][SIZE], int nx, int ny, double w, int *iterations);

int main(void) {
    char filename[20];
    int nx, ny, num, iterations;
    double Lx, Ly, hx, hy, tol, err;
    double x[SIZE], y[SIZE];
    double Te[SIZE][SIZE], Told[SIZE][SIZE];
    clock_t start, end;
    double cpu_time_used;

    // Initialize domain and grid parameters
    initialize_domain(&Lx, &Ly, &nx, &ny, &hx, &hy);

    // Initialize vectors x and y
    initialize_vectors(x, y, hx, hy, nx);

    // Initialize temperature field
    initialize_temperature(Te, nx, ny);

    tol = 1e-12;
    err = 1.0;
    iterations = 0;

    printf("Input the method of approximation:");
    printf("\n[1] Jacobi method");
    printf("\n[2] Gauss-Seidel method");
    printf("\n[3] SOR method\n");
    scanf("%d", &num);

    // Start the clock
    start = clock();

    // Choose the method based on user input
    while (err >= tol) {
        switch (num) {
            case 1:
                err = jacobi_method(Te, Told, nx, ny, &iterations);
                break;
            case 2:
                err = gauss_seidel_method(Te, Told, nx, ny, &iterations);
                break;
            case 3:
                err = sor_method(Te, Told, nx, ny, 1.813958, &iterations); // "w" found in a separate program where 0 < w < 2
                break;
            default:
                printf("Invalid method selected.\n");
                return 1;
        }
    }

    // End the clock
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

    printf("err = %E\n", err);
    printf("Number of iterations = %d\n", iterations);

    // Output results to CSV files
    matrix_csv("temps", Te, SIZE, SIZE);
    vector_csv("xy", x, y, SIZE);

    printf("\n\nChosen method took %4e seconds.\n", cpu_time_used);
    return 0;
}

void matrix_csv(const char *filename, double a[][SIZE], int n, int m) {
    printf("\nCreating %s.csv file for matrix.", filename);
    
    FILE *fp;
    int i, j;
    char full_filename[256];
    
    // Create full filename with ".csv" extension
    snprintf(full_filename, sizeof(full_filename), "%s.csv", filename);
    
    // Open file for writing
    fp = fopen(full_filename, "w");
    if (fp == NULL) {
        perror("Error opening file");
        return;
    }
    
    // Write matrix to file
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            // Print element and check if it's the last in the row
            if (j == n - 1) {
                fprintf(fp, "%10.14f\n", a[i][j]);
            } else {
                fprintf(fp, "%10.14f, ", a[i][j]);
            }
        }
    }
    
    // Close file
    fclose(fp);
    printf("\n%s file created.\n", full_filename);
}

void vector_csv(const char *filename, double v[SIZE], double u[SIZE], int n) {
    printf("\nCreating vector %s.csv file for x & y.", filename);
    FILE *fp;
    char full_filename[100];
    snprintf(full_filename, sizeof(full_filename), "%s.csv", filename);
    fp = fopen(full_filename, "w+");
    if (fp == NULL) {
        printf("Error opening file %s\n", full_filename);
        return;
    }

    for (int i = 0; i < n; i++) {
        fprintf(fp, "%f, %f\n", v[i], u[i]);
    }
    fclose(fp);
    printf("\nVector %s file created.", full_filename);
}

void initialize_domain(double *Lx, double *Ly, int *nx, int *ny, double *hx, double *hy) {
    *Lx = 6.0;
    *Ly = 4.0;
    *nx = SIZE;
    *ny = SIZE;
    *hx = *Lx / (*nx - 1);
    *hy = *Ly / (*ny - 1);
}

void initialize_vectors(double x[SIZE], double y[SIZE], double hx, double hy, int nx) {
    for (int i = 0; i < nx; i++) {
        x[i] = i * hx;
        y[i] = i * hy;
    }
}

void initialize_temperature(double Te[SIZE][SIZE], int nx, int ny) {
    // Define boundary conditions
    double T_L = 300.0; // Left boundary temperature
    double T_T = 500.0; // Top boundary temperature
    double T_R = 300.0; // Right boundary temperature
    double T_B = 500.0; // Bottom boundary temperature

    // Initialize the temperature field to 0.0
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            Te[i][j] = 0.0;
        }
    }

    // Set the boundary conditions
    for (int i = 1; i < nx - 1; i++) {
        Te[i][0]      = T_L;    // Left boundary
        Te[i][nx - 1] = T_R;    // Right boundary
        Te[0][i]      = T_T;    // Top boundary
        Te[nx - 1][i] = T_B;    // Bottom boundary
    }

    // Calculate the average temperature at the corners
    Te[0][0] = (T_T + T_L) / 2.0;           // Top-left corner
    Te[0][nx - 1] = (T_T + T_R) / 2.0;      // Top-right corner
    Te[ny - 1][0] = (T_L + T_B) / 2.0;      // Bottom-left corner
    Te[ny - 1][nx - 1] = (T_R + T_B) / 2.0; // Bottom-right corner
}

double jacobi_method(double Te[SIZE][SIZE], double Told[SIZE][SIZE], int nx, int ny, int *iterations) {
    double err = 0.0;
    double matrix[SIZE][SIZE];

    do {
        err = 0.0;
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                Told[i][j] = Te[i][j];
            }
        }

        for (int i = 1; i < nx - 1; i++) {
            for (int j = 1; j < ny - 1; j++) {
                Te[i][j] = 0.25 * (Told[i - 1][j] + Told[i + 1][j] + Told[i][j - 1] + Told[i][j + 1]);
            }
        }

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                matrix[i][j] = Te[i][j] - Told[i][j];
                if (matrix[i][j] > err) {
                    err = matrix[i][j];
                }
            }
        }
        (*iterations)++;
    } while (err >= 1e-12);

    return err;
}

double gauss_seidel_method(double Te[SIZE][SIZE], double Told[SIZE][SIZE], int nx, int ny, int *iterations) {
    double err = 0.0;
    double matrix[SIZE][SIZE];

    do {
        err = 0.0;
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                Told[i][j] = Te[i][j];
            }
        }

        for (int i = 1; i < nx - 1; i++) {
            for (int j = 1; j < ny - 1; j++) {
                Te[i][j] = 0.25 * (Te[i - 1][j] + Told[i + 1][j] + Te[i][j - 1] + Told[i][j + 1]);
            }
        }

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                matrix[i][j] = Te[i][j] - Told[i][j];
                if (matrix[i][j] > err) {
                    err = matrix[i][j];
                }
            }
        }
        (*iterations)++;
    } while (err >= 1e-12);

    return err;
}

double sor_method(double Te[SIZE][SIZE], double Told[SIZE][SIZE], int nx, int ny, double w, int *iterations) {
    double err = 0.0;
    double matrix[SIZE][SIZE];

    do {
        err = 0.0;
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                Told[i][j] = Te[i][j];
            }
        }

        for (int i = 1; i < nx - 1; i++) {
            for (int j = 1; j < ny - 1; j++) {
                Te[i][j] = (1 - w) * Told[i][j] + w * 0.25 * (Te[i - 1][j] + Told[i + 1][j] + Te[i][j - 1] + Told[i][j + 1]);
            }
        }

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                matrix[i][j] = Te[i][j] - Told[i][j];
                if (matrix[i][j] > err) {
                    err = matrix[i][j];
                }
            }
        }
        (*iterations)++;
    } while (err >= 1e-12);

    return err;
}
