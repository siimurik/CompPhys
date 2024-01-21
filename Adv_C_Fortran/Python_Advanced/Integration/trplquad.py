from scipy import integrate
import numpy as np

# Define the integrand function
def integrand(x, y, z):
    return x**2 + y**2 - z**2

# Define the limits of integration
x_lower = -1
x_upper = 1
y_lower = lambda x: -np.sqrt(1 - x**2)
y_upper = lambda x: np.sqrt(1 - x**2)
z_lower = lambda x, y: -np.sqrt(1 - x**2 - y**2)
z_upper = lambda x, y: np.sqrt(1 - x**2 - y**2)

# Perform the triple integration
result, error = integrate.tplquad(integrand, x_lower, x_upper, y_lower, y_upper, z_lower, z_upper)

# Display the result
print("Result of triple integration:", result)
val = 0.83776623474001310#0.837766290
print("Difference:", abs(result-val))