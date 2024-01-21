import numpy as np
from scipy.integrate import quad
def integrand(x):
    return np.exp(-x**2)
    #return np.cos(1/x)/np.sqrt(x)

a = -2.0
b =  2.0
I = quad(integrand, a, b)

print("Result of the integration: ", I)
fort_val = 1.7641627815248431 
# qgaus    1.7641627815248431 
# romb     1.7641627815248437
# trapzd   1.7641627815247811
# quad     1.7641627815247982E+00
# trap     1.7641627815246661 
print("Difference", abs(fort_val-I[0]))