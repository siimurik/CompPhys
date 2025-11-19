#include <iostream>
#include <cassert>
#include <fstream>

/*
    $ g++ ex3.3.cpp -o ex33
    $ ./ex33 101
*/

double dy(double y)
{
    return -y;      // dy/dx = -y
}

int main(int argc, char* argv[]){

    // Check if command line arguments exist
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] 
        << "\t<num_of_grid_points> (e.g. 11, 101, ...)" << "\n";
        return 1;
    }

    // Convert command line argument to integer
    int N = std::atoi(argv[1]);

    // Assert that number of grid points is greater than 1
    assert(N > 1);

    // Calculate step size
    double h = 1.0 / (N - 1);

    // Initial conditions
    double x = 0.0;
    double y = 1.0; // y(0) = 1.0

    // Open output file
    std::ofstream output_file("xy.dat");
    assert(output_file.is_open());
    output_file.setf(std::ios::scientific);
    output_file.precision(6);

    // Write initial point to file
    output_file << x << "  " << y << "\n";

    // Implicit Euler method 
    for (int n = 1; n < N; n++)
    {
        x = n * h;
        // Implicit Euler: y_n = y_{n-1} / (1 + h)
        // From: (y_n - y_{n-1})/h = -y_n
        // Rearranged: y_n = y_{n-1} / (1 + h)
        y = y / (1.0 + h);

        output_file << x << "  " << y << "\n";
    }

    output_file.close();

    std::cout << "Data written to xy.dat\n";
    std::cout << "Step size h = " << h << "\n";
    std::cout << "Number of grid points N = " << N << "\n";

    return 0;
}