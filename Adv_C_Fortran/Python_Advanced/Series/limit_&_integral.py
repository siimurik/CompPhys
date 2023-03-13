import sympy as sy
import numpy as np
import scipy.integrate as integrate

x = sy.symbols("x")

r = sy.limit(3*x/((x**2-2)*sy.log(x)), x, sy.oo)

if (r == sy.oo):
    print("The given sequence does not converge to any number")
else:
    print("The given sequence converges to a number " + str(r))

print("Limit of the function")
print(r)

result = integrate.quad(lambda x: 3*x/((x**2-2)*np.log(x)), 2, np.inf)
print("Integral of the function")
print(result)