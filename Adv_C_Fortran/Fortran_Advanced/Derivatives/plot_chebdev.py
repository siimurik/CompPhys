import numpy as np
import matplotlib.pyplot as plt

# Read Chebyshev derivative values from the data file
data = np.loadtxt('chebyshev_derivative.dat')

# Extract the x values and corresponding derivative values
x_values = np.linspace(-2.0, 2.0, len(data))
derivative_values = data

print(x_values)

## Plot the result
#plt.plot(x_values, derivative_values, label='Chebyshev Derivative')
#plt.title('Chebyshev Derivative of exp(-x^2)')
#plt.xlabel('x')
#plt.ylabel('Derivative Value')
#plt.legend()
#plt.grid()
#plt.show()
