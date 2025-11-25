#include <iostream>

// Write code that sends the address of an integer to a 
// function that prints out the value of the integer.

void printIntValue(int* ptr)
{
    if (ptr != nullptr)
    {
        std::cout << "The value of the integer is: " << *ptr << "\n";
    }
    else
    {
        std::cout << "Null pointer received.\n";
    }
}

int main()
{
    // Send the address of the integer to the function
    int myInt = 42;
    printIntValue(&myInt);

    // Example with dynamic allocation
    int* dynInt = new int(100);
    printIntValue(dynInt);

    delete dynInt; // Clean up dynamically allocated memory

    // Example with null pointer
    int* nullPtr = nullptr;
    printIntValue(nullPtr);

    return 0;
}