"""
Step 6: 2-D Convection

Now we solve 2D Convection, represented by the pair of 
coupled partial differential equations below:
                ∂u/∂t + u∂/u∂x + v∂u/∂y = 0
                ∂v/∂t + u∂/v∂x + v∂v/∂y = 0

Discretizing these equations using the methods we've 
applied previously yields:

(un+1i,j−uni,j)/Δt + (uni,juni,j−uni−1,j)/Δx + vni,juni,j−uni,j−1)/Δy = 0
(vn+1i,j−vni,j)/Δt + (uni,jvni,j−vni−1,j)/Δx + vni,jvni,j−vni,j−1)/Δy = 0

Rearranging both equations, we solve for un+1i,j
and vn+1i,j, respectively. Note that these equations are also coupled.

un+1i,j=uni,j − ui,jΔt/Δx(uni,j−uni−1,j) − vni,jΔt/Δy(uni,j−uni,j−1)
vn+1i,j=vni,j − ui,jΔt/Δx(vni,j−vni−1,j) − vni,jΔt/Δy(vni,j−vni,j−1)
Initial Conditions

The initial conditions are the same that we used for 1D convection, 
applied in both the x and y directions.

            [
            | 2 for x,y ∈ (0.5,1)×(0.5,1)
u, v =      {
            | 1     everywhere else
            [

Boundary Conditions

The boundary conditions hold u and v equal to 1 along the boundaries
of the grid.
u=1, v=1 for {x=0,2y=0,2

                [
                | x = 0, 2
u=1, v=1 for    {
                | y = 0, 2
                [
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D    ##New Library required for projected 3d pl
#from matplotlib import cm

### variable declarations
nx = 101
ny = 101
nt = 80
c = 1
dx = 2 / (nx-1)
dy = 2 / (ny-1)
sigma = 0.2
dt = sigma * dx

x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)
X, Y = np.meshgrid(x, y)

u  = np.ones((nx, ny)) ## create a 1xn vector of 1's
v  = np.ones((nx, ny))
un = np.ones((nx, ny))
vn = np.ones((nx, ny))
#print('dx =', dx, '\nu =', np.size(u))

### Assign initial condidtions
### set hat function I.C.: u(0.5<=x<=1 && 0.5<=y<=1 ) is 2
u[int(0.5 / dx):int(1/dx + 1), int(0.5 / dy):int(1/dy + 1)] = 2
### set hat function I.C.: v(0.5<=x<=1 && 0.5<=y<=1 ) is 2
v[int(0.5 / dx):int(1/dx + 1), int(0.5 / dy):int(1/dy + 1)] = 2
"""
### Plot Initial Conditions
## the figsize parameter can be used to produce differeent size images
fig = plt.figure(figsize=(11, 7), dpi = 100)
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, u, cmap=plt.cm.viridis, rstride=2, cstride=2)
ax.set_xlabel('$x$')
ax.set_ylabel('$x$')
plt.show()
"""

for n in range(nt + 1): ## loop across number of time steps
    un = u.copy()
    vn = v.copy()
    u[1:,1:] = (un[1:,1:] - 
                (un[1:,1:]*c*dt/dx*(un[1:,1:]-un[:-1,1:])) -
                (vn[1:,1:]*c*dt/dy*(un[1:,1:]-un[1:,:-1])))
    v[1:,1:] = (vn[1:,1:] - 
                (un[1:,1:]*c*dt/dx*(vn[1:,1:]-vn[:-1,1:])) -
                (vn[1:,1:]*c*dt/dy*(vn[1:,1:]-vn[1:,:-1])))
    u[0, :]  = 1
    u[-1, :] = 1
    u[:, 0]  = 1
    u[:, -1] = 1
    
    v[0, :]  = 1
    v[-1, :] = 1
    v[:, 0]  = 1
    v[:, -1] = 1

fig = plt.figure(figsize=(11, 7), dpi = 100)
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, u, cmap=plt.cm.viridis, rstride=2, cstride=2)
ax.set_xlabel('$x$')
ax.set_ylabel('$x$')
plt.show()