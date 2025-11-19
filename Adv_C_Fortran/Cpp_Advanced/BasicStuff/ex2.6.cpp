#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

// Function f(x) = e^x + x^3 - 5
double f(double x) {
    return std::exp(x) + std::pow(x, 3) - 5.0;
}

// Derivative f'(x) = e^x + 3x^2
double f_prime(double x) {
    return std::exp(x) + 3.0 * std::pow(x, 2);
}

int main() {
    std::cout << std::fixed << std::setprecision(10);
    
    // Part 1: Newton-Raphson iteration for f(x) = e^x + x^3 - 5
    // x_i = x_{i-1} - f(x_{i-1}) / f'(x_{i-1})
    //     = x_{i-1} - (e^{x_{i-1}} + x_{i-1}^3 - 5) / (e^{x_{i-1}} + 3x_{i-1}^2)
    
    // Part 2: Using for loop with array storage
    std::cout << "=== Part 2: Newton-Raphson with array storage ===\n";
    std::vector<double> x(101); // Array for iterates x_0 to x_100
    x[0] = 0.0; // Initial guess x_0 = 0
    
    for (int i = 1; i <= 100; i++) {
        x[i] = x[i-1] - f(x[i-1]) / f_prime(x[i-1]);
        std::cout << "Iteration " << i << ": x_" << i << " = " << x[i] 
                  << ", f(x_" << i << ") = " << f(x[i]) << "\n";
    }
    
    // Part 3: Check - verify that f(x_final) is close to zero
    std::cout << "\n=== Part 3: Verification ===\n";
    double final_x = x[100];
    double final_f = f(final_x);
    std::cout << "Final solution: x = " << final_x << "\n";
    std::cout << "f(x) at final solution: " << final_f << "\n";
    std::cout << "|f(x)| should be close to 0: " << std::abs(final_f) << "\n";
    
    if (std::abs(final_f) < 1e-10) {
        std::cout << "✓ Check passed: f(x) is very close to zero!\n";
    } else {
        std::cout << "✗ Check failed: f(x) is not close enough to zero.\n";
    }
    
    // Part 4: Using scalar variables instead of array
    std::cout << "\n=== Part 4: Newton-Raphson with scalar variables ===\n";
    double x_prev = 0.0; // x_0
    double x_next;
    
    std::cout << "Iteration 0: x_prev = " << x_prev << ", f(x_prev) = " << f(x_prev) << "\n";
    
    for (int i = 1; i <= 10; i++) {
        x_next = x_prev - f(x_prev) / f_prime(x_prev);
        std::cout << "Iteration " << i << ": x_next = " << x_next 
                  << ", f(x_next) = " << f(x_next) << "\n";
        x_prev = x_next;
    }
    
    // Part 5: Using while loop with epsilon tolerance
    std::cout << "\n=== Part 5: Newton-Raphson with while loop and epsilon ===\n";
    
    // Test different epsilon values
    std::vector<double> epsilons = {1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12};
    
    for (double epsilon : epsilons) {
        std::cout << "\n--- Testing with epsilon = " << epsilon << " ---\n";
        
        x_prev = 0.0; // Reset initial guess
        int iteration = 0;
        double difference;
        
        std::cout << "Iteration " << iteration << ": x = " << x_prev 
                  << ", f(x) = " << f(x_prev) << "\n";
        
        do {
            x_next = x_prev - f(x_prev) / f_prime(x_prev);
            iteration++;
            difference = std::abs(x_next - x_prev);
            
            std::cout << "Iteration " << iteration << ": x = " << x_next 
                      << ", f(x) = " << f(x_next) 
                      << ", |Δx| = " << difference << "\n";
            
            x_prev = x_next;
            
        } while (difference > epsilon && iteration < 1000);
        
        std::cout << "Converged after " << iteration << " iterations\n";
        std::cout << "Final solution: x = " << x_next << ", f(x) = " << f(x_next) << "\n";
        
        // Additional check: verify f(x) is close to zero
        if (std::abs(f(x_next)) < 1e-10) {
            std::cout << "✓ Root verification: f(x) ≈ 0\n";
        } else {
            std::cout << "✗ Root verification failed: f(x) = " << f(x_next) << "\n";
        }
    }
    
    return 0;
}