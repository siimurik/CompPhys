import numpy as np
from scipy.optimize import fsolve

# Define the system of equations
def system(vars):
    x1, x2 = vars
    f1 = x1**2 + x2**2 - 4
    f2 = x1 * x2 - 1
    return [f1, f2]

# Initial guess
x0 = [2.0, 0.5]

# Solve the system
solution, info, ier, msg = fsolve(system, x0, full_output=True)

# Evaluate residuals
residuals = system(solution)

# Output the results
print("----------------------------------------")
print(" NONLINEAR SYSTEM SOLVER RESULTS")
print("----------------------------------------")
print(f" Solution:        x1 = {solution[0]:12.6f}, x2 = {solution[1]:12.6f}")
print(f" Residuals:       f1 = {residuals[0]:12.4e}, f2 = {residuals[1]:12.4e}")
print(f" Function evals:  {info['nfev']} calls")
print(f" Solver status:   {ier} - {msg}")
print("----------------------------------------")

if np.max(np.abs(residuals)) < 1e-6:
    print(" ✅ Verification: solution satisfies tolerance!")
else:
    print(" ⚠️  Warning: residuals above tolerance level.")
