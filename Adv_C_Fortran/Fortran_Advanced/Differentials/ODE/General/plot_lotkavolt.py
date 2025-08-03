import numpy as np
import matplotlib.pyplot as plt
import os

# Get list of files
file_name = 'lv.dat'

# Read data from file
data = np.loadtxt(file_name)

# Extract third and fourth columns
t = data[:, 1]
x = data[:, 2]
y = data[:, 3]

# Plot data with specified color
plt.plot(t, x, label='x')
plt.plot(t, y, label='y')

# Add labels and title
plt.xlabel('time')
plt.ylabel('x/y')
plt.title('Lotka-Volterra predator-prey model')

# Add legend and grid
plt.legend()
plt.grid(True)

# Show plot
plt.show()
