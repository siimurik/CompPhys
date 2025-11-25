#include <iostream>

// 5.3 Write a function that swaps the values of two double precision floating point
// numbers, so that these changes are visible in the code that has called this function.
// 1. Write this function using pointers.
// 2. Write this function using references.

void swapUsingPointers(double *a, double *b)
{
    if (a != nullptr && b != nullptr)
    {
        double temp = *a;
        *a = *b;
        *b = temp;
    }
    else
    {
        std::cout << "Null pointer received.\n";
    }
}

void swapUsingReferences(double &a, double &b)
{
    double temp = a;
    a = b;
    b = temp;
}

int main()
{
    double x = 5.0;
    double y = 10.0;

    std::cout << "Before swapUsingPointers: x = " << x << ", y = " << y << "\n";
    swapUsingPointers(&x, &y);
    std::cout << "After swapUsingPointers: x = " << x << ", y = " << y << "\n";

    std::cout << "Before swapUsingReferences: x = " << x << ", y = " << y << "\n";
    swapUsingReferences(x, y);
    std::cout << "After swapUsingReferences: x = " << x << ", y = " << y << "\n";

    return 0;
}