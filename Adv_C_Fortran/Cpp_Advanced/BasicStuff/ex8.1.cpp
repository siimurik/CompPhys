#include <iostream>
#include <cassert>
#include <iomanip>

// Template class for a specialized array with bounds checking
template <class T>
class ProbabilityArray
{
private:
    T* mData;           // Pointer to array data
    int mSize;          // Size of the array
    const T mTolerance; // Tolerance for numerical errors (10^-6)

public:
    // Constructor
    ProbabilityArray(int size) : mSize(size), mTolerance(1e-6)
    {
        assert(size > 0);
        mData = new T[size];
        // Initialize all values to 0
        for (int i = 0; i < mSize; i++)
        {
            mData[i] = 0.0;
        }
    }

    // Destructor
    ~ProbabilityArray()
    {
        delete[] mData;
    }

    // Copy constructor
    ProbabilityArray(const ProbabilityArray& otherArray)
        : mSize(otherArray.mSize), mTolerance(otherArray.mTolerance)
    {
        mData = new T[mSize];
        for (int i = 0; i < mSize; i++)
        {
            mData[i] = otherArray.mData[i];
        }
    }

    // Assignment operator
    ProbabilityArray& operator=(const ProbabilityArray& otherArray)
    {
        if (this != &otherArray)
        {
            delete[] mData;
            mSize = otherArray.mSize;
            mData = new T[mSize];
            for (int i = 0; i < mSize; i++)
            {
                mData[i] = otherArray.mData[i];
            }
        }
        return *this;
    }

    // Get size of array
    int GetSize() const
    {
        return mSize;
    }

    // Overload [] operator for writing (setting values)
    T& operator[](int index)
    {
        assert(index >= 0 && index < mSize);
        return mData[index];
    }

    // Overload [] operator for reading (getting values with correction)
    T operator[](int index) const
    {
        assert(index >= 0 && index < mSize);
        
        T value = mData[index];
        
        // Case 1: Value is between 0 and 1 inclusive - return as is
        if (value >= 0.0 && value <= 1.0)
        {
            return value;
        }
        // Case 2: Value is between -10^-6 and 0 (exclusive) - return 0
        else if (value >= -mTolerance && value < 0.0)
        {
            return 0.0;
        }
        // Case 3: Value is between 1 (exclusive) and 1 + 10^-6 - return 1
        else if (value > 1.0 && value <= 1.0 + mTolerance)
        {
            return 1.0;
        }
        // Case 4: Value is outside acceptable range - trigger assertion
        else
        {
            std::cerr << "Error: Value at index " << index 
                      << " is " << value 
                      << ", which is outside the acceptable range ["
                      << -mTolerance << ", " << (1.0 + mTolerance) << "]" << std::endl;
            assert(false && "Probability value outside acceptable range!");
            return 0.0; // This line will never be reached
        }
    }

    // Alternative method: Read with explicit correction
    T Read(int index) const
    {
        return (*this)[index];  // Uses the const [] operator
    }

    // Method to set a value directly
    void Write(int index, T value)
    {
        assert(index >= 0 && index < mSize);
        mData[index] = value;
    }

    // Print the array (showing both raw and corrected values)
    void Print() const
    {
        std::cout << std::fixed << std::setprecision(10);
        std::cout << "Index | Raw Value          | Corrected Value" << std::endl;
        std::cout << "------|--------------------|-----------------" << std::endl;
        for (int i = 0; i < mSize; i++)
        {
            std::cout << std::setw(5) << i << " | "
                      << std::setw(18) << mData[i] << " | "
                      << std::setw(15) << (*this)[i] << std::endl;
        }
    }
};

// Main function to demonstrate the ProbabilityArray class
int main()
{
    std::cout << "=== Probability Array Demonstration ===" << std::endl << std::endl;

    // Create a probability array for 10 days
    const int N = 10;
    ProbabilityArray<double> rainProbabilities(N);

    std::cout << "Setting probability values (some with numerical errors):" << std::endl << std::endl;

    // Set various test values
    rainProbabilities[0] = 0.5;              // Valid: between 0 and 1
    rainProbabilities[1] = 0.0;              // Valid: exactly 0
    rainProbabilities[2] = 1.0;              // Valid: exactly 1
    rainProbabilities[3] = -0.0000005;       // Needs correction: between -10^-6 and 0
    rainProbabilities[4] = 1.0000005;        // Needs correction: between 1 and 1+10^-6
    rainProbabilities[5] = 0.25;             // Valid
    rainProbabilities[6] = 0.999999;         // Valid
    rainProbabilities[7] = 0.000001;         // Valid
    rainProbabilities[8] = -0.0000001;       // Needs correction (within tolerance)
    rainProbabilities[9] = 1.0000001;        // Needs correction (within tolerance)

    // Display the array with corrections
    rainProbabilities.Print();

    std::cout << std::endl << "=== Testing Individual Access ===" << std::endl;
    std::cout << "Day 0: " << rainProbabilities.Read(0) << std::endl;
    std::cout << "Day 3 (corrected from negative): " << rainProbabilities.Read(3) << std::endl;
    std::cout << "Day 4 (corrected from >1): " << rainProbabilities.Read(4) << std::endl;

    // Test with a value outside acceptable range (this will trigger assertion)
    std::cout << std::endl << "=== Testing Assertion (Uncomment to test) ===" << std::endl;
    std::cout << "The following line is commented out because it would terminate the program:" << std::endl;
    std::cout << "rainProbabilities[9] = -0.01;  // Too negative!" << std::endl;
    std::cout << "rainProbabilities[9] = 1.01;   // Too large!" << std::endl;
    
    // Uncomment ONE of these lines to test assertion:
    // rainProbabilities.Write(9, -0.01);      // This will fail assertion
    // std::cout << rainProbabilities.Read(9) << std::endl;
    
    // OR:
    // rainProbabilities.Write(9, 1.01);       // This will fail assertion
    // std::cout << rainProbabilities.Read(9) << std::endl;

    std::cout << std::endl << "=== Testing with Float Template ===" << std::endl;
    ProbabilityArray<float> floatProbabilities(3);
    floatProbabilities[0] = 0.3f;
    floatProbabilities[1] = -0.0000005f;
    floatProbabilities[2] = 1.0000005f;
    
    std::cout << "Float array values:" << std::endl;
    for (int i = 0; i < 3; i++)
    {
        std::cout << "Index " << i << ": " << floatProbabilities[i] << std::endl;
    }

    std::cout << std::endl << "Program completed successfully!" << std::endl;

    return 0;
}