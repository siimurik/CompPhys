"""
In the following exercises, you will extend the first four 
steps to 2D. To extend the 1D finite-difference formulas to
partial derivatives in 2D or 3D, just apply the definition: 
a partial derivative with respect to x is the variation in 
the x direction at constant y.

In 2D space, a rectangular (uniform) grid is defined by the 
points with coordinates:
                x(i) = x(0) + iΔx
                y(i) = y(0) + iΔy

Now, define ui,j=u(xi,yj) and apply the finite-difference 
formulas on either variable x,y acting separately on the i 
and j indices. All derivatives are based on the 2D Taylor 
expansion of a mesh point value around ui,j.

Hence, for a first-order partial derivative in the x-direction,
a finite-difference formula is:
        ∂u/∂x∣i,j = (ui+1,j−ui,j)/Δx + O(Δx)

and similarly in the y direction. Thus, we can write backward-
difference, forward-difference or central-difference formulas 
for Steps 5 to 12. Let's get started!

Step 5: 2-D Linear Convection

The PDE governing 2-D Linear Convection is written as
                ∂u/∂t + c ∂u/∂x + c ∂u/∂y = 0

This is the exact same form as with 1-D Linear Convection, 
except that we now have two spatial dimensions to account 
for as we step forward in time.

Again, the timestep will be discretized as a forward difference
and both spatial steps will be discretized as backward differences.

With 1-D implementations, we used i subscripts to denote movement
in space (e.g. u^n_i − u^n_i−1). Now that we have two dimensions 
to account for, we need to add a second subscript, j, to account 
for all the information in the regime.

Here, we'll again use i as the index for our x values, and we'll 
add the j subscript to track our y values.

With that in mind, our discretization of the PDE should be relatively
straightforward. 

    (un+1i,j−uni,j)/Δt + c(uni,j−uni−1,j)/Δx + c(uni,j−uni,j−1)/Δy = 0

As before, solve for the only unknown:
u^n+1_i,j = uni,j − cΔt/Δx(uni,j−uni−1,j) − cΔt/Δy(uni,j−uni,j−1)

We will solve this equation with the following initial conditions:
u(x,y)={2 for1 for0.5≤x,y≤1everywhere else

and boundary conditions:

            [
            | x = 0, 2
u = 1 for   {
            | y = 0, 2
            [
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D    ##New Library required for projected 3d pl
#from matplotlib import cm

### variable declarations
nx = 81
ny = 81
nt = 100
c = 1
dx = 2 / (nx-1)
dy = 2 / (ny-1)
sigma = 0.2
dt = sigma * dx

x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)
X, Y = np.meshgrid(x, y)

u  = np.ones((ny, nx)) ## create a 1xn vector of 1's
un = np.ones((ny, nx)) ##
#print('dx =', dx, '\nu =', np.size(u))

### Assign initial condidtions

### set hat function I.C.: u(0.5<=x<=1 && 0.5<=y<=1 ) is 2
u[int(0.5 / dy):int(1/dy + 1), int(0.5 / dx):int(1/dx + 1)] = 2
"""
### Plot Initial Conditions
## the figsize parameter can be used to produce differeent size images
fig = plt.figure(figsize=(11, 7), dpi = 100)
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, u[:], cmap=plt.cm.viridis)
plt.show()
"""
"""
for n in range(nt + 1): ## loop across number of time steps
    un = u.copy()
    row, col = u.shape
    for i in range(1, row):
        for j in range(1, col):
            u[i,j] = (un[i,j] - (c*dt/dx*(un[i,j]-un[i-1,j])) -
                                (c*dt/dy*(un[i,j]-un[i,j-1])))
            u[ 0, :] = 1
            u[-1, :] = 1
            u[:,  0] = 1
            u[:, -1] = 1

fig = plt.figure(figsize=(11, 7), dpi = 100)
ax = fig.add_subplot(111, projection='3d')
surf2 = ax.plot_surface(X, Y, u[:], cmap=plt.cm.viridis)
plt.show()
"""
# Basically same but without the double for loops
for n in range(nt + 1): ## loop across number of time steps
    un = u.copy()
    u[1:,1:] = (un[1:,1:] - (c*dt/dx*(un[1:,1:]-un[:-1,1:])) -
                            (c*dt/dy*(un[1:,1:]-un[1:,:-1])))
    u[ 0, :] = 1
    u[-1, :] = 1
    u[:,  0] = 1
    u[:, -1] = 1

fig = plt.figure(figsize=(11, 7), dpi = 100)
ax = fig.add_subplot(111, projection='3d')
surf2 = ax.plot_surface(X, Y, u[:], cmap=plt.cm.viridis)
plt.show()    