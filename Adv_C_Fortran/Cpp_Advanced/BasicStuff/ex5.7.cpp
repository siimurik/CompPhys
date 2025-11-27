#include <cmath>
#include <iostream>
#include <iomanip>

/*
5.7 The p-norm of a vector v of length n is given by
    ||v||_p = (Σ_i^n |v_i|^p)^(1/p) for p >= 1
where p is a positive integer. Extend the code in Sect.5.10 to calculate the p-norm
of a given vector, where p takes the default value 2.
*/

void    freeVector(double* v);
double* allocateVector(int size);
double  CalculateNorm(double* v, int n, int p);
void    printVec(double* v, int size, const std::string& name);

int main()
{

    const int n = 3;
    double* v = allocateVector(n);
    v[0] = 3.0;
    v[1] = -4.0;
    v[2] = 12.0;

    // Calculate and print 2-norm (default)
    double norm2 = CalculateNorm(v, n, 2);
    printVec(v, n, "v");
    std::cout << "2-norm of v: " << norm2 << "\n";

    // Free memory
    freeVector(v);

    return 0;
}

double CalculateNorm(double* v, int n, int p) 
{
    double a = 0.0;
    for (int i = 0; i < n; i++) 
    {
        double temp = fabs(v[i]);
        a += pow(temp, p);
    }
    return pow(a, 1.0/p);
}

// Print vector
void printVec(double* v, int size, const std::string& name)
{
    std::cout << name << ":\n";
    for (int i = 0; i < size; i++)
    {
        std::cout << std::fixed << std::setprecision(2) << std::setw(6) << v[i] << " ";
    }
    std::cout << "\n\n";
}

// Function to allocate vector memory
double* allocateVector(int size)
{
    return new double[size];
}

// Deallocate vector memory
void freeVector(double* v)
{
    delete[] v;
}

