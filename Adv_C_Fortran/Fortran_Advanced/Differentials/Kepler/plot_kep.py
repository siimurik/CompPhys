import numpy as np
import matplotlib.pyplot as plt
import os

# Read data from file
data = np.loadtxt('kep.dat')

# Extract third and fourth columns
x = data[:, 2]
y = data[:, 3]

plt.plot(x, y, label='Trajectory')
# Add legend and grid
plt.legend()
plt.grid(True)

# Show plot
plt.show()