"""
Now we're going to implement nonlinear convection 
using the same methods as in step 1. The 1D convection
equation is:

                ∂u/∂t + u ∂u/∂x = 0

Instead of a constant factor c multiplying the second
term, now we have the solution u multiplying it. Thus,
the second term of the equation is now nonlinear.
We're going to use the same discretization as in Step 1 —
forward difference in time and backward difference in
space. Here is the discretized equation.

(u^n+1_i - u^n_i)/Δt + u^n_i(u^n_i - u^n_i-1)/Δx = 0

Solving for the only unknown term, u^n+1_i, yields:

u^n+1_i = u^n_i - u^n_i Δt/Δx (u^n_i - u^n_i-1)

As before, the Python code starts by loading the 
necessary libraries. Then, we declare some variables
that determine the discretization in space and time 
(you should experiment by changing these parameters 
to see what happens). Then, we create the initial 
condition u0 by initializing the array for the 
solution using u=2 @ 0.5≤x≤1 and u=1 everywhere 
else in (0,2) (i.e., a hat function).
"""
import numpy as np
import matplotlib.pyplot as plt

nx = 41
dx = 2 / (nx - 1)
nt = 20    #nt is the number of timesteps we want to calculate
dt = .025  #dt is the amount of time each timestep covers (delta t)

u = np.ones(nx)      #as before, we initialize u with every value equal to 1.
u[int(.5 / dx) : int(1 / dx + 1)] = 2  #then set u = 2 between 0.5 and 1 as per our I.C.s

fig, ax = plt.subplots()
plt.plot(np.linspace(0, 2, nx), u)

un = np.ones(nx) #initialize our placeholder array un, to hold the time-stepped solution

for n in range(nt):  #iterate through time
    un = u.copy() ##copy the existing values of u into un
    for i in range(1, nx):  ##now we'll iterate through the u array
    
     ###This is the line from Step 1, copied exactly.  Edit it for our new equation.
     ###then uncomment it and run the cell to evaluate Step 2   
      
           ###u[i] = un[i] - c * dt / dx * (un[i] - un[i-1]) 
           u[i] = un[i] - un[i] * dt/dx * (un[i] - un[i-1]) 

plt.plot(np.linspace(0, 2, nx), u) ##Plot the results
plt.grid()
plt.show()