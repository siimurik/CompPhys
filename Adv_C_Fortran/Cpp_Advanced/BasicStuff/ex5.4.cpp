#include <iostream>
#include <cmath>

// 5.4 Write a function that can be used to calculate the mean and standard deviation of
// an array of double precision floating point numbers. Note that the standard deviation
// σ of a collection of numbers x j, j = 1, 2,..., N is given by
// σ = sqrt( (1/(N-1) * Σ (x j - x¯)² )
// where x¯ is the mean of the numbers.

void calcMeanAndStdDev(double* arr, int size, double* mean, double* stddev)
{
    if (arr == nullptr || size <= 0 || mean == nullptr || stddev == nullptr)
    {
        std::cout << "Invalid input.\n";
        return;
    }

    // Calculate mean
    double sum = 0.0;
    for (int i = 0; i < size; i++)
    {
        sum += arr[i];
    }
    *mean = sum / size;

    // Calculate standard deviation
    double sumSqDiff = 0.0;
    for (int j = 0; j < size; j++)
    {
        sumSqDiff += (arr[j] - *mean) * (arr[j] - *mean);
    }
    *stddev = std::sqrt(sumSqDiff / (size - 1));
}

int main()
{   
    double* x;
    x = new double[5];
    x[0] = 10.0;
    x[1] = 12.0;
    x[2] = 23.0;
    x[3] = 23.0;
    x[4] = 16.0;
    double* y = new double[5]{87.0, 100.0, 94.0, 91.0, 87.0};
    
    double mean = 0.0;
    double stddev = 0.0;

    calcMeanAndStdDev(x, 5, &mean, &stddev);
    std::cout << "Mean: " << mean << "\n";
    std::cout << "Standard Deviation: " << stddev << "\n";

    calcMeanAndStdDev(y, 5, &mean, &stddev);
    std::cout << "Mean: " << mean << "\n";
    std::cout << "Standard Deviation: " << stddev << "\n";

    delete[] x;
    delete[] y;

    return 0;
}