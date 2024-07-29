# Heat Conduction Problem Solver

This program solves the steady 2D heat conduction problem using three different numerical methods: Jacobi, Gauss-Seidel, and Successive Over-Relaxation (SOR). The results are outputted to CSV files for further analysis.

## Files
- `fdm_FINAL.c`: The main C source code file for the heat conduction solver.

## Compilation and Execution

To compile and run the program, use the following commands in the terminal:

```
gcc fdm_FINAL.c -o fdm
./fdm
```

## Program Workflow

1. **Initialization:**
   - The domain and grid parameters are initialized.
   - Vectors `x` and `y` representing the grid points are initialized.
   - The temperature field is initialized with boundary conditions.

2. **Method Selection:**
   - The user is prompted to select a numerical method for solving the heat conduction problem:
     1. Jacobi method
     2. Gauss-Seidel method
     3. Successive Over-Relaxation (SOR) method

3. **Iteration:**
   - The chosen method iteratively updates the temperature field until the error is less than the specified tolerance.
   - The number of iterations and the final error are recorded.

4. **Output:**
   - The final temperature distribution is saved to `temps.csv`.
   - The `x` and `y` vectors are saved to `xy.csv`.

## Solution Methods and Mathematics

### Heat Conduction Equation

The steady-state heat conduction equation in two dimensions is given by the Laplace equation:

$$ \frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0 $$

This partial differential equation is discretized using finite difference methods.

### Discretization

For a grid with spacing $ h $ in both $ x $ and $ y $ directions, the discretized form of the Laplace equation at a grid point $(i, j)$ is:

$$ T_{i,j} = \frac{1}{4} (T_{i-1,j} + T_{i+1,j} + T_{i,j-1} + T_{i,j+1}) $$

This equation forms the basis for the iterative solution methods.

### Jacobi Method

The Jacobi method updates the temperature at each grid point using the average of the temperatures from the previous iteration:

$$ T_{i,j}^{(k+1)} = \frac{1}{4} (T_{i-1,j}^{(k)} + T_{i+1,j}^{(k)} + T_{i,j-1}^{(k)} + T_{i,j+1}^{(k)}) $$

where $ k $ denotes the iteration number. The process is repeated until convergence.

### Gauss-Seidel Method

The Gauss-Seidel method improves upon the Jacobi method by using the latest available temperature values during the update:

$$ T_{i,j}^{(k+1)} = \frac{1}{4} (T_{i-1,j}^{(k+1)} + T_{i+1,j}^{(k)} + T_{i,j-1}^{(k+1)} + T_{i,j+1}^{(k)}) $$

This typically results in faster convergence compared to the Jacobi method.

### Successive Over-Relaxation (SOR) Method

The SOR method introduces a relaxation factor $ \omega $ to accelerate convergence:

$$ T_{i,j}^{(k+1)} = (1 - \omega) T_{i,j}^{(k)} + \frac{\omega}{4} (T_{i-1,j}^{(k+1)} + T_{i+1,j}^{(k)} + T_{i,j-1}^{(k+1)} + T_{i,j+1}^{(k)}) $$

The optimal choice of $ \omega $ depends on the specific problem and can significantly reduce the number of iterations required for convergence.

## Function Documentation

### `void matrix_csv(const char *filename, double a[][SIZE], int n, int m)`

Creates a CSV file with the specified filename and writes the matrix `a` to it.

### `void vector_csv(const char *filename, double v[SIZE], double u[SIZE], int n)`

Creates a CSV file with the specified filename and writes the vectors `v` and `u` to it.

### `void initialize_domain(double *Lx, double *Ly, int *nx, int *ny, double *hx, double *hy)`

Initializes the domain dimensions, number of grid points, and grid spacing.

### `void initialize_vectors(double x[SIZE], double y[SIZE], double hx, double hy, int nx)`

Initializes the `x` and `y` vectors representing the grid points.

### `void initialize_temperature(double Te[SIZE][SIZE], int nx, int ny)`

Initializes the temperature field `Te` with boundary conditions.

### `double jacobi_method(double Te[SIZE][SIZE], double Told[SIZE][SIZE], int nx, int ny, int *iterations)`

Solves the heat conduction problem using the Jacobi method. Returns the final error and updates the number of iterations.

### `double gauss_seidel_method(double Te[SIZE][SIZE], double Told[SIZE][SIZE], int nx, int ny, int *iterations)`

Solves the heat conduction problem using the Gauss-Seidel method. Returns the final error and updates the number of iterations.

### `double sor_method(double Te[SIZE][SIZE], double Told[SIZE][SIZE], int nx, int ny, double w, int *iterations)`

Solves the heat conduction problem using the SOR method with relaxation factor `w`. Returns the final error and updates the number of iterations.

## Usage Example

1. Compile and run the program.
2. Choose a method (1, 2, or 3) when prompted.
3. The program will iterate until convergence is achieved.
4. The final error, number of iterations, and CPU time will be displayed.
5. Results will be saved to `temps.csv` and `xy.csv`.

## Notes

- The grid size is defined by the `SIZE` macro, currently set to 61.
- The relaxation factor for the SOR method is pre-determined and set to 1.813958.
- The tolerance for convergence is set to `1e-12`.

## License

This code is provided under the MIT License. See the `LICENSE` file for details.

## References

- Original article: [Skill-Lync Project](https://skill-lync.com/student-projects/week-5-mid-term-project-solving-the-steady-and-unsteady-2d-heat-conduction-problem-35)
