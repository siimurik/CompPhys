import numpy as np
from scipy.integrate import tplquad

# Define the integrand function
def func(z, y, x):
    return x**2 + y**2 - z**2

# Define the limits for z
def z1(x, y):
    return -np.sqrt(1.0 - x**2 - y**2)

def z2(x, y):
    return np.sqrt(1.0 - x**2 - y**2)

# Define the limits for y
def y1(x):
    return -np.sqrt(1.0 - x**2)

def y2(x):
    return np.sqrt(1.0 - x**2)

# Define the limits for x
x1 = -1.0
x2 = 1.0

# Perform the triple integration
result, error = tplquad(func, x1, x2, y1, y2, z1, z2, epsabs=1.0E-14, epsrel=1.0E-14)

print(f"Result: {result}")
print(f"Estimated error: {error}")
