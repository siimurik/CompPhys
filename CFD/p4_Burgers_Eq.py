"""
Step 4: Burgers' Equation

More on the topic: https://en.wikipedia.org/wiki/Burgers'_equation

Burgers' equation in one spatial dimension looks like this:

            ∂u/∂t + u ∂u/∂x = ν ∂²u/∂x²

As you can see, it is a combination of non-linear convection
and diffusion. It is surprising how much you learn from this 
neat little equation!

We can discretize it using the methods we've already detailed 
in Steps 1 to 3. Using forward difference for time, backward 
difference for space and our 2nd-order method for the second 
derivatives yields:

(u^n+1_i - u^n_i)/Δt + u^n_i(u^n_i-u^n_i-1)/Δx = ν (u^n_i+1-2u^n_i+u^n_i-1)/Δx²

As before, once we have an initial condition, the only unknown 
is u^n+1_i. We will step in time as follows:

u^n+1_i = u^n_i - u^n_i Δt/Δx(u^n_i-u^n_i-1) + ν Δt/Δx²(u^n_i+1-2u^n_i+u^n_i-1)

Initial and Boundary Conditions

To examine some interesting properties of Burgers' equation, 
it is helpful to use different initial and boundary conditions 
than we've been using for previous steps.

Our initial condition for this problem is going to be:

    u = -2ν/ϕ ∂ϕ/∂x + 4                 (1)
    ϕ = exp(-x²/4ν) + exp(-(x-2π)²/4ν)  (2)

This has an analytical solution, given by:
    u = -2ν/ϕ ∂ϕ/∂x + 4                                             (3)
    ϕ = exp( -(x-4t)²/(4ν(t+1)) ) + exp( -(x-4t-2π)²/(4ν(t+1)) )    (4)

Our boundary condition will be:

            u(0) = u(2π)

This is called a periodic boundary condition. 
Pay attention! This will cause you a bit of 
headache if you don't tread carefully.

Saving Time with SymPy

The initial condition we're using for Burgers' Equation
can be a bit of a pain to evaluate by hand. The derivative
∂ϕ/∂x isn't too terribly difficult, but it would be easy 
to drop a sign or forget a factor of x somewhere, so we're 
going to use SymPy to help us out.

SymPy is the symbolic math library for Python. It has a lot 
of the same symbolic math functionality as Mathematica with 
the added benefit that we can easily translate its results 
back into our Python calculations (it is also free and open source).

Start by loading the SymPy library, 
together with our favorite library, NumPy.
"""
import numpy as np
import sympy as sy
import matplotlib.pyplot as plt
from sympy import init_printing
init_printing(use_latex=True)

x, nu, t = sy.symbols('x nu t')
phi = (sy.exp(-(x - 4 * t)**2 / (4 * nu * (t + 1))) +
       sy.exp(-(x - 4 * t - 2 * sy.pi)**2 / (4 * nu * (t + 1))))
sy.pprint(phi)

"""
It's maybe a little small, but that looks right. Now to evaluate our partial derivative ∂ϕ∂x
is a trivial task."""
phiprime = phi.diff(x)
#sy.pprint(phiprime)

"""
Now that we have the Pythonic version of our derivative,
we can finish writing out the full initial condition 
equation and then translate it into a usable Python 
expression. For this, we'll use the lambdify function, 
which takes a SymPy symbolic equation and turns it 
into a callable function.
"""
u = -2 * nu * (phiprime / phi) + 4
#print(u)

"""
To lambdify this expression into a useable function, 
we tell lambdify which variables to request and the 
function we want to plug them in to.
"""
ufunc = sy.lambdify((t, x, nu), u)
#print(ufunc(1, 4, 3))

"""
Now that we have the initial conditions set up, 
we can proceed and finish setting up the problem. 
We can generate the plot of the initial condition 
using our lambdify-ed function.
"""
# variable declarations
nx = 101
nt = 100
dx = 2*np.pi / (nx - 1)
nu = 0.07
dt = dx * nu

x  = np.linspace(0, 2*np.pi, nx)
un = np.empty(nx)
t  = 0

u = np.asarray([ufunc(t, x0, nu) for x0 in x])
#print(u)

"""fig, ax = plt.subplots(figsize=(11,7))
plt.plot(x, u, marker='o', lw=2)
plt.xlim([0, 2 * np.pi])
plt.ylim([0, 10])"""

for n in range(nt):
    un = u.copy()
    for i in range(1, nx-1):
        u[i] = un[i] - un[i] * dt / dx *(un[i] - un[i-1]) + nu * dt / dx**2 *\
                (un[i+1] - 2 * un[i] + un[i-1])
    u[0] = un[0] - un[0] * dt / dx * (un[0] - un[-2]) + nu * dt / dx**2 *\
                (un[1] - 2 * un[0] + un[-2])
    u[-1] = u[0]
        
u_analytical = np.asarray([ufunc(nt * dt, xi, nu) for xi in x])

fig, ax = plt.subplots(figsize=(11,7))
#plt.figure(figsize=(11, 7), dpi=100)
plt.plot(x,u, marker='o', lw=2, label='Computational')
plt.plot(x, u_analytical, label='Analytical')
plt.xlim([0, 2 * np.pi])
plt.ylim([0, 10])
plt.legend()
plt.grid()
plt.show()

