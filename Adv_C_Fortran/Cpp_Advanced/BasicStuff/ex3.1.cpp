#include <iostream>
#include <fstream>
#include <cassert>

int main (int argc, char* argv[])
{
    double x[4] = {0.0, 1.0, 1.0, 0.0};
    double y[4] = {0.0, 0.0, 1.0, 1.0};

    // Step 4.
    std::ifstream check_file("x_and_y.dat");
    if (check_file.is_open())
    {
        std::cout << "Warning: x_and_y.dat already exists!\n";
        check_file.close();
        
        char choice;
        std::cout << "Do you want to (e)rase or (a)ppend to the existing file?\n";
        std::cin >> choice;

        if (choice == 'a' || choice == 'A')
        {
            // Append mode - open file for appending
            std::ofstream write_output("x_and_y.dat", std::ios::app);
            assert(write_output.is_open());
            write_output.setf(std::ios::scientific);
            write_output.setf(std::ios::showpos);
            write_output.precision(10);
            for (int i = 0; i < 4; i++)
            {
                write_output << x[i] << "  " << y[i] << "\n";
                write_output.flush(); // Data safety: If the program crashes, all previously written lines are safely stored in the file
            }
            write_output.close();
            std::cout << "Data appended to existing file.\n";
            return 0;
        }
        // If user chooses 'erase', continue normally
    }

    // Normal writing (erase existing file or create new one)
    std::ofstream write_output("x_and_y.dat");
    assert(write_output.is_open());
    write_output.setf(std::ios::scientific);
    write_output.setf(std::ios::showpos);
    write_output.precision(10);

    for (int i = 0; i < 4; i++)
    {
        write_output << x[i] << "  " << y[i] << "\n";
        write_output.flush(); // Data safety: If the program crashes, all previously written lines are safely stored in the file
    }
    write_output.close();
    std::cout << "New file created (or existing file overwritten).\n";

    // 2. Extend the code so that the output stream is flushed immediately after each line
    // of the file is written
    // Added write_output.flush() after each write operation in the loop

    // 3. Extend the code so that the precision is set to 10 significant figures, the output is
    // in scientific notation, and plus signs are shown for positive numbers.

    // 4. Amend the program so that it does not automatically create a fresh file
    // x_and_y.dat every time it is run. Have the program first attempt to open
    // the file x_and_y.dat as an ifstream for reading. If the file can be successfully opened then, after closing the ifstream, warn the user. Have the program
    // prompt the user as to whether it should erase the existing file or append to the
    // existing file

    return 0;
}