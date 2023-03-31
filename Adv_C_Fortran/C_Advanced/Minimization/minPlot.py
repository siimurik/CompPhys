import numpy as np
import matplotlib.pyplot as plt

# Define the paraboloid function
def paraboloid(x, y, p):
    return p[2] * (x - p[0]) ** 2 + p[3] * (y - p[1]) ** 2 + p[4]

# Define the parameters
par = [1.0, 2.0, 10.0, 20.0, 30.0]

# Generate a grid of points to evaluate the function
x = np.linspace(-10, 10, 100)
y = np.linspace(-10, 10, 100)
X, Y = np.meshgrid(x, y)
Z = paraboloid(X, Y, par)

# Find the minimum point using the C code
# 24  9.920e-01  1.997e+00 f() =  30.001 size = 0.008
x_min = 9.920e-01
y_min = 1.997e+00
f_min = 30.001

# Plot the function and the minimum point
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')
ax.scatter(x_min, y_min, f_min, c='r', marker='o', s=100)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('f(x, y)')
plt.show()
