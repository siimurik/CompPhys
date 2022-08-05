"""
The 1-D Linear Convection equation is the simplest, 
most basic model that can be used to learn something
about CFD. It is surprising that this little equation 
can teach us so much! 
Here it is:
                ∂u/∂t + c ∂u/∂x = 0

With given initial conditions (understood as a wave),
the equation represents the propagation of that 
initial wave with speed c, without change of shape. 
Let the initial condition be u(x,0)=u0(x). Then the 
exact solution of the equation is 

                u(x,t) = u0(x-ct)

We discretize this equation in both space and time, 
using the Forward Difference scheme for the time 
derivative and the Backward Difference scheme for the
space derivative. Consider discretizing the spatial
coordinate x into points that we index from i=0 to N,
and stepping in discrete time intervals of size Δt.

From the definition of a derivative (and simply 
removing the limit), we know that:

                ∂u/∂x ≈ (u(x+Δx) - u(x)) / Δx

Our discrete equation, then, is:

(u^(n+1)_i - u^n_i)/Δt + c(u^n_i - u^n_(i-1))/Δx = 0

Where n and n+1 are two consecutive steps in time, 
while i-1 and i are two neighboring points of the 
discretized x coordinate. If there are given initial
conditions, then the only unknown in this discretization
is u^(n+1)_i. We can solve for our unknown to get an
equation that allows us to advance in time, as follows:

u^(n+1)_(i) = u^n_i - c Δt/Δx (u^n_i - u^n_i-1)

"""
from matplotlib import axes
import numpy as np
import matplotlib.pyplot as plt
import time, sys 

# Declaration of variables 
nx = 41  # try changing this number from 41 to 81 and Run All ... what happens?
dx = 2 / (nx-1)
nt = 25     #nt is the number of timesteps we want to calculate
dt = 0.025  #dt is the amount of time each timestep covers (delta t)
c = 1       #assume wavespeed of c = 1

# Initial conditions
u = np.ones(nx)      #numpy function ones()
u[int(.5 / dx):int(1 / dx + 1)] = 2  #setting u = 2 between 0.5 and 1 as per our I.C.s
#print(u,'\n', len(u))

# Initial plot
fig, ax = plt.subplots()
plt.plot(np.linspace(0, 2, nx), u, label='Initial condition')

un = np.ones(nx) #initialize a temporary array

for n in range(nt):  #loop for values of n from 0 to nt, so it will run nt times
    un = u.copy() ##copy the existing values of u into un
    for i in range(1, nx): ## you can try commenting this line and...
    #for i in range(nx): ## ... uncommenting this line and see what happens!
        u[i] = un[i] - c * dt / dx * (un[i] - un[i-1])

plt.plot(np.linspace(0, 2, nx), u, label='Progressed wave after nt seconds')
plt.legend(loc='upper right')
plt.grid()
plt.show()



