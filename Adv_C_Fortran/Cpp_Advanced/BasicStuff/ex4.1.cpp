#include <iostream>

int main() 
{
    // Declare an integer i and set it to 5
    int i = 5;

    // Declare a pointer to an integer p_j and 
    // store the address of i this pointer
    int* p_j = &i;
    
    // Multiply the value of the variable i by 
    // 5 using only the pointer variable
    *p_j = *p_j * 5;
    
    // Declare another pointer to another integer
    // p_k and use new to allocate memory
    int* p_k = new int;

    // Store the contents of the variable i in 
    // this location
    *p_k = i;

    // Change the vale pointed to vy p_j to 0
    *p_j = 0;

    // Output the values to check correctness
    std::cout << "Value of i: " << i << "\n";
    std::cout << "Value pointed to by p_j: " << *p_j << "\n";
    std::cout << "Value pointed to by p_k: " << *p_k << "\n";

    // Clean up dynamically allocated memory
    delete p_k;

    return 0;
}