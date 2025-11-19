#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{
    std::ifstream read_file("x_and_y.dat");
    if (!read_file.is_open())
    {
        return 1;
    }
    
    int number_of_rows = 0;
    double dummy1, dummy2;
    
    while(read_file >> dummy1 >> dummy2)
    {
        number_of_rows++;
    }

    std::cout << "Number of rows = "
              << number_of_rows << "\n";
    read_file.close();
    
    return 0;
}

/*
Run the code above. This code does not give the correct answer. Why is this?
Does the code give the correct answer if the final newline character is removed from
the file x_and_y.dat? Modify the code so that it gives the correct answer.
[Hint: You might investigate the use of read_file.fail() which may be used to
probe whether the last read on the file stream was unsuccessful.]
*/