
#include <iostream>

/*  g++ -O3 ex4.3.cpp -o ex43  */

int main() 
{
    const long long ITERATIONS = 1000000000;

    for (long long i = 0; i < ITERATIONS; i++)
    {
        // Dynamically allocate memory for two vectors 
        // of length 3
        double* vec1 = new double[3];
        double* vec2 = new double[3];

        // Assign values to each entry
        vec1[0] = 1.0; vec1[1] = 2.0; vec1[2] = 3.0;
        vec2[0] = 4.0; vec2[1] = 5.0; vec2[2] = 6.0;

        // Calculate the scalar (dot) product
        double dot_prod = 0.0;
        for (int j = 0; j < 3; j++)
        {
            dot_prod += vec1[j] * vec2[j];
        }

        // Print the dot prod (careful... (ó﹏ò｡) )
        //std::cout << "Dot product: " << dot_prod << "\n";

        // De-allocate the memory properly
        delete[] vec1;
        delete[] vec2;
    }

    std::cout << "Completed " << ITERATIONS <<
    " iterations without memory leaks.\n";

    return 0;
}