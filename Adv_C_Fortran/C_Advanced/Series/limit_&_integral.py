import sympy as sy
import numpy as np
import scipy.integrate as integrate
from scipy.special import gamma, factorial
########################################################################
x, n = sy.symbols("x n")
r = sy.limit( ( sy.Abs(3*n/((n**2-2)*sy.log(n))) )**(1/n), n, sy.oo)
#r = 1/((x+1)*(3*x-2))
#x = -4
#r = sy.limit( (x+3)**(n+1)/((n+1)+2)/((x+3)**n/(n+2)), n, sy.oo)
if (r == sy.oo):
    print("The given sequence does not converge to any number")
else:
    print("The given sequence converges to a number " + str(r))

print("Limit of the function")
print(r)
########################################################################
result = integrate.quad(lambda x: 3*x/((x**2-2)*np.log(x)), 2, np.inf, epsabs=1e-10, epsrel=1e-10, limit=100)
#result = integrate.quad(lambda n: (x+3)**n/(n+2), 1, np.inf)
#result = integrate.quad(lambda x: 1/((x+1)*(3*x-2)), 1, np.inf)
#result = integrate.quad(lambda x: np.power(2, x)/gamma(x+1.0), 0, np.inf)
print("Integral of the function")
print(result)
########################################################################
