#include <iostream>
#include <iomanip>

// Function to print any 1D array with two decimal places
void print_vector(const double vec[], int size, const std::string& name = "vec") {
    for (int i = 0; i < size; i++) {
        std::cout << name << "[" << i << "] = " 
                  << std::fixed << std::setprecision(2) << vec[i] << "\n";
    }
    std::cout << "\n";  // Add spacing between outputs
}

// Template version for any 3x3 matrix
template<size_t Rows, size_t Cols>
void print_matrix(const double (&mat)[Rows][Cols], const std::string& name = "mat") {
    std::cout << name << ":\n";
    for (int i = 0; i < Rows; i++) {
        std::cout << "  ";
        for (int j = 0; j < Cols; j++) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(2) << mat[i][j];
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

/////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    double u[3] = {1.0, 2.0, 3.0};
    double v[3] = {6.0, 5.0, 4.0};
    double A[3][3] = {
        {1.0, 5.0, 0.0},
        {7.0, 1.0, 2.0},
        {0.0, 0.0, 1.0}
    };
    double B[3][3] = {
        {-2.0, 0.0, 1.0},
        { 1.0, 0.0, 0.0},
        { 4.0, 1.0, 0.0}
    };

    double w[3] = {0.0};
    int i;
    for (i=0; i<3; i++) {
        w[i] = u[i] - 3.0*v[i];
    }

    // Use the function instead of repeated loops
    print_vector(w, 3, "w");

    //-------------------------------------------
    // Exercises:
    // x = u - v
    double x[3] = {0.0};
    for (i = 0; i < 3; i++) {
        x[i] = u[i] - v[i];
    }
    print_vector(x, 3, "x");

    // y = Au
    double y[3] = {0.0}; 
    int j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            y[i] += A[i][j] * u[j];  // Example: matrix-vector multiplication
        }
    }
    print_vector(y, 3, "y");

    //z = Au - v,
    //y = Au
    double z[3] = {0.0};
    for (i = 0; i < 3; i++){
        z[i] = y[i] - v[i];
    }
    print_vector(z, 3, "z");

    // C = 4A - 3B
    double C[3][3] = {0.0};
    for (i=0; i<3; i++){
        for(j=0; j<3; j++){
            C[i][j] = 4.0*A[i][j] - 3.0*B[i][j];
        }
    }
    // Write a function for printing a matrix 
    print_matrix(C, "C = 4A - 3B");

    // D = AB
    double D[3][3] = {0.0}; 
    int k;
    for (i = 0; i < 3; i++){
        for (j = 0; j < 3; j++){
            for (k = 0; k < 3; k++){
                D[i][j] += A[i][k] * B[k][j]; 
            }
        }
    }
    print_matrix(D, "D = A*B");
    
    return 0;
}