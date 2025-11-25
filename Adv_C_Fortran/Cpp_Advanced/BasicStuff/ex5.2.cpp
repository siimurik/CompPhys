#include <iostream>

// Write code that sends the address of an integer to a 
// function that changes the value of the integer.

void changeIntValue(int* ptr, int newValue)
{
    if (ptr != nullptr)
    {
        *ptr = newValue;
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
    std::cout << "Original value of the integer: " << myInt << "\n";
    changeIntValue(&myInt, 100);
    std::cout << "Changed value of the integer: " << myInt << "\n";

    // Example with dynamic allocation
    int* dynInt = new int(200);
    std::cout << "Original value of the dynamically allocated integer: " << *dynInt << "\n";
    changeIntValue(dynInt, 300);
    std::cout << "Changed value of the dynamically allocated integer: " << *dynInt << "\n";

    delete dynInt; // Clean up dynamically allocated memory

    return 0;
}