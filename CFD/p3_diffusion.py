"""
Step 3: Diffusion Equation in 1-D

The one-dimensional diffusion equation is:
            ∂u/∂t = ν ∂²u/∂x²

The first thing you should notice is that —
unlike the previous two simple equations we 
have studied— this equation has a second-order
derivative. We first need to learn what to 
do with it!

Discretizing ∂²u/∂x²

The second-order derivative can be represented
geometrically as the line tangent to the curve 
given by the first derivative. We will discretize
the second-order derivative with a Central Difference
scheme: a combination of Forward Difference and 
Backward Difference of the first derivative. 
Consider the Taylor expansion of u_i+1
and u_i−1 around ui

u_i+1 = u_i + Δx ∂u/∂x∣_i + Δx²/2 ∂²u/∂x²∣_i + Δx³/3! ∂³u/∂x³∣_i + O(Δx⁴)

u_i-1 = u_i - Δx ∂u/∂x∣_i + Δx²/2 ∂²u/∂x²∣_i - Δx³/3! ∂³u/∂x³∣_i + O(Δx⁴)

If we add these two expansions, you can see 
that the odd-numbered derivative terms will 
cancel each other out. If we neglect any 
terms of O(Δx⁴) or higher (and really, those 
are very small), then we can rearrange the 
sum of these two expansions to solve for our
second-derivative.

u_i+1 + u_i-1 = 2u_i + Δx² ∂²u/∂x²∣_i + O(Δx⁴)

Then rearrange to solve for ∂²u/∂x²∣_i and 
the result is:

  ∂²u/∂x² = (u_i+1 - 2u_i + u_i-1)/Δx² + O(Δx2)

Back to Step 3

We can now write the discretized version 
of the diffusion equation in 1D:

(u^n+1_i - u^n_i) / Δt = ν(u^n_i+1 - 2u^n_i + u^n_i-1)/Δx²

As before, we notice that once we have an 
initial condition, the only unknown is u^n+1_i,
so we re-arrange the equation solving for our unknown:

u^n+1_i = u^n_i + νΔt/Δx²(u^n_i+1 - 2u^n_i + u^n_i-1)

The above discrete equation allows us to write 
a program to advance a solution in time. But we 
need an initial condition. Let's continue using 
our favorite: the hat function. So, at t=0, u=2 
in the interval 0.5≤x≤1 and u=1 everywhere else.
We are ready to number-crunch!
"""
import numpy as np
import matplotlib.pyplot as plt

nx = 41
dx = 2 / (nx - 1)
nt = 20     # the number of timesteps we want to calculate
nu = 0.3    # the value of viscosity
sigma = 0.2 # magical parameter
dt = sigma * dx**2

u = np.ones(nx)     # a numpy array with nx elements all equal to 1.
u[int(0.5 / dx):int(1/dx + 1)] = 2  # setting u =2 between 0.5 and 1 as per our I.C.s
                                    # In layman's terms [10:21] = 2.

un = np.ones(nx)    # our placeholder array, un, to advance the solution in time

for n in range(nt): # iterate throgh time
    un = u.copy()   ## copy the existing values of u into un
    for i in range(1, nx-1):
        u[i] = un[i] + nu*dt/dx**2*(un[i+1] - 2*un[i] + un[i-1])

fig, ax = plt.subplots()
plt.plot(np.linspace(0, 2, nx), u)
plt.grid()
plt.show()