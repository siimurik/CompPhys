#include <iostream>

int main()
{
    int a = 10;
    int b = 20;

    std::cout << "Before swap: \n";
    std::cout << "a = " << a << ", b = " << b << "\n";

    int* p_a = &a;
    int* p_b = &b;

    // XOR swap using only pointers
    *p_a = *p_a ^ *p_b;
    *p_b = *p_a ^ *p_b;
    *p_a = *p_a ^ *p_b;

    std::cout << "After swap: \n";
    std::cout << "a = " << a << ", b = " << b << "\n";

    return 0;
}