import sympy as sp
import numpy as np
from scipy.special import i0, i1
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt


print("\n=== PART 1: EXACT SYMBOLIC DERIVATION (SymPy) ===")

x = sp.Symbol('x', real=True, positive=True)

# Define numerator and denominator separately to avoid the NotImplementedError
num = x * (1 - x**2) * sp.besseli(1, x)
den = sp.besseli(0, x)

# 1. Expand each to a high order polynomial (O(x^10))
terms = 13
num_poly = num.series(x, 0, terms).removeO()
den_poly = den.series(x, 0, terms).removeO()

# 2. Use the series expansion of the ratio manually
# This is mathematically: (numerator) * (1/denominator)
F_series = (num_poly / den_poly).series(x, 0, terms).removeO()

print(f"\n1. Correct Taylor Series of F(x) (up to x^{terms-1} terms):")
print(sp.simplify(sp.simplify(F_series)))

# 3. Take the symbolic derivative
dF_series = sp.diff(F_series, x)
print("\n2. Derivative of the Series:")
print(dF_series)

# 4. Find the roots
roots = sp.solve(dF_series, x)
valid_roots = [r.evalf() for r in roots if r.is_real and r > 0]
print(f"\nValid roots: {valid_roots}")

# second derivative test
d2F_series = sp.diff(dF_series, x)
print("\nSecond derivative of the Series:")
print(d2F_series)
for r in valid_roots:
    d2F_val = d2F_series.subs(x, r).evalf()
    print(f"At root {r:.4f}, d2F={d2F_val:.4f}")
    if d2F_val < 0:
        print(f"Root {r:.4f} is a local maximum (d2F={d2F_val:.4f})")

# Filter to physically meaningful domain: 0 < x < 1
valid_roots = [
    r.evalf() for r in roots
    if r.is_real and 0 < float(r.evalf()) < 1
]

if valid_roots:
    x_max_sympy = min(valid_roots, key=lambda r: abs(float(r)))
    print(f"\n3. Peak of the Taylor approximation (x_max): {float(x_max_sympy):.4f}")

print("\n\n=== PART 2: EXACT NUMERICAL SOLUTION (SciPy) ===\n")

def exact_F(x_val):
    # Negated for minimization
    return -(x_val * (1 - x_val**2) * i1(x_val) / i0(x_val))

# Using the bounded method to find the peak between 0 and 1
result = minimize_scalar(exact_F, bounds=(0.01, 0.99), method='bounded')

x_exact_max = result.x
F_exact_max = -result.fun

print(f"1. True Maximum Wavenumber (x_max): {x_exact_max:.6f}")
print(f"2. True Maximum Growth Rate (F_max): {F_exact_max:.6f}")

# Final physical spacing
r_0 = 2.0  # mm
lambda_max = (2 * np.pi * r_0) / x_exact_max
print(f"\nExact drop spacing for r0=2mm: {lambda_max:.2f} mm")

print("\n\n=== PART 3: EXACT NUMERICAL SOLUTION (Secant Method) ===\n")

F = lambda x: x * (1 - x**2) * i1(x) / i0(x)

# F(x) = x*(1-x²)*I₁(x) / I₀(x)
# Using quotient rule: d/dx [u/v] = (u'v - uv') / v²
# where u = x*(1-x²)*I₁(x),  v = I₀(x)
#
# u' = (1-3x²)*I₁(x) + x*(1-x²)*I₁'(x)
# I₁'(x) = I₀(x) - I₁(x)/x   (Bessel recurrence)
# v' = I₀'(x) = I₁(x)
# So:
# dF/dx = [((1-3x²)*I₁ + x*(1-x²)*(I₀ - I₁/x)) * I₀ - (x*(1-x²)*I₁) * I₁] / I₀²
def dF_exact(x_val):
    """Correct analytical derivative of F(x) = x(1-x²)I₁(x)/I₀(x)"""
    i0x = i0(x_val)
    i1x = i1(x_val)
    i1p = i0x - i1x / x_val          # I₁'(x) = I₀(x) - I₁(x)/x
    u   = x_val * (1 - x_val**2) * i1x
    up  = (1 - 3*x_val**2) * i1x + x_val * (1 - x_val**2) * i1p
    vp  = i1x                         # I₀'(x) = I₁(x)
    return (up * i0x - u * vp) / i0x**2

def secant(f, x0, x1, tol=1e-10, max_iter=100):
    """
    Secant method to find root of f(x) = 0.
    Apply to f = dF_exact to find the maximum of F.

    Parameters
    ----------
    f           : function - the function whose root we seek (pass dF_exact here)
    x0, x1      : floats   - two initial guesses where f has opposite signs
    tol         : float    - convergence tolerance on |xm - x_new|
    max_iter    : int      - maximum number of iterations

    Returns
    -------
    x_new : float - root of f, i.e. location of maximum of F
    """
    if f(x0) * f(x1) >= 0:
        print("Secant method fails: f(x0) and f(x1) must have opposite signs.")
        return None

    print(f"{'Iter':>4}  {'x0':>12}  {'x1':>12}  {'x_new':>12}  {'F(x_new)':>12}  {'F`(x_new)':>12}")
    print("-" * 72)

    n     = 0
    x_new = 0

    while True:
        # Secant step
        x_new  = (x0 * f(x1) - x1 * f(x0)) / (f(x1) - f(x0))
        fp_new = f(x_new)
        c      = f(x0) * fp_new

        print(f"{n+1:>4}  {x0:>12.8f}  {x1:>12.8f}  {x_new:>12.8f}  {F(x_new):>12.8f}  {fp_new:>12.2e}")

        # convergence check: how much did the estimate move this step?
        if abs(x_new - x1) < tol:
            x0 = x1
            x1 = x_new
            n += 1
            break

        # update interval
        x0 = x1
        x1 = x_new
        n += 1

        # exact root found
        if c == 0:
            print(f"\nFound exact solution at iteration {n}.")
            break

        if n >= max_iter:
            print(f"Reached maximum iterations ({max_iter}).")
            break

    print("-" * 72)
    print(f"\nConverged in {n} iterations.")
    return x_new  # <-- was missing in original; caused NoneType error downstream


# Bracket must straddle the maximum at x≈0.697 without including x=0.
# (dF also vanishes at x=0, so a wide bracket like [0.01, 0.99] can
#  overshoot and converge to that trivial root instead.)
x_max = secant(dF_exact, x0=0.5, x1=0.7, tol=1e-10, max_iter=100)

if x_max is not None:
    F_max = F(x_max)
    print(f"\n1. True Maximum Wavenumber (x_max): {x_max:.6f}")
    print(f"2. True Maximum Growth Rate (F_max): {F_max:.6f}")

# Plot of F(x) and dF/dx to visually confirm the maximum location
x_vals = np.linspace(0.01, 0.99, 500)
F_vals = F(x_vals)
dF_vals = dF_exact(x_vals)
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(x_vals, F_vals, label='F(x)')
plt.axvline(x_max, color='r', linestyle='--', label=f'Max at x={x_max:.4f}')
plt.title('Growth Rate F(x)')
plt.xlabel('x')
plt.ylabel('F(x)')
plt.legend()
plt.subplot(1, 2, 2)
plt.plot(x_vals, dF_vals, label="dF/dx")
plt.axhline(0, color='k', linestyle='--')
plt.axvline(x_max, color='r', linestyle='--', label=f'Root at x={x_max:.4f}')
plt.title('Derivative dF/dx')
plt.xlabel('x')
plt.ylabel('dF/dx')
plt.legend()
plt.tight_layout()
plt.show()