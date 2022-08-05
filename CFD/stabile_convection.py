"""
What happened?

To answer that question, we have to think a little 
bit about what we're actually implementing in code.

In each iteration of our time loop, we use the 
existing data about our wave to estimate the speed of
the wave in the subsequent time step. Initially, the
increase in the number of grid points returned more 
accurate answers. There was less numerical diffusion
and the square wave looked much more like a square 
wave than it did in our first example.

Each iteration of our time loop covers a time-step of
length Δt, which we have been defining as 0.025.

During this iteration, we evaluate the speed of the 
wave at each of the x points we've created. In the 
last plot, something has clearly gone wrong.

What has happened is that over the time period Δt,
the wave is travelling a distance which is greater 
than dx. The length dx of each grid box is related 
to the number of total points nx, so stability can 
be enforced if the Δt

step size is calculated with respect to the size of dx.

                σ = u Δt / Δx ≤ σ_max

where u is the speed of the wave; σ is called the 
Courant number and the value of σ_max that will 
ensure stability depends on the discretization used.

In a new version of our code, we'll use the CFL 
number to calculate the appropriate time-step dt 
depending on the size of dx.

-----

Notice that as the number of points nx increases, the
wave convects a shorter and shorter distance. The 
number of time iterations we have advanced the 
solution at is held constant at nt = 20, but 
depending on the value of nx and the corresponding 
values of dx and dt, a shorter time window is being 
examined overall.
"""
import numpy as np                  #numpy is a library for array operations akin to MATLAB
import matplotlib.pyplot as plt     #matplotlib is 2D plotting library

def linearconv(nx):
    dx = 2 / (nx - 1)
    nt = 20    #nt is the number of timesteps we want to calculate
    c = 1
    sigma = 0.5
    
    dt = sigma * dx

    u = np.ones(nx)      #defining a numpy array which is nx elements long with every value equal to 1.
    u[int(.5/dx):int(1 / dx + 1)] = 2  #setting u = 2 between 0.5 and 1 as per our I.C.s

    un = np.ones(nx) #initializing our placeholder array, un, to hold the values we calculate for the n+1 timestep

    for n in range(nt):  #iterate through time
        un = u.copy() ##copy the existing values of u into un
        for i in range(1, nx):
            u[i] = un[i] - c * dt / dx * (un[i] - un[i-1])
    

    plt.plot(np.linspace(0, 2, nx), u)

fig, ax = plt.subplots()
linearconv(41)
linearconv(61)
linearconv(81)
linearconv(101)
plt.grid()
plt.show()
        