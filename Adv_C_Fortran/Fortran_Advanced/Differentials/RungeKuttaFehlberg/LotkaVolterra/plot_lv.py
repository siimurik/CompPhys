import numpy as np
import matplotlib.pyplot as plt
import os

# Get list of files
file_names = ['lv.dat', 'lv2.dat', 'lv3.dat', 'predator_prey.dat']
colors = ['#2F6690', '#16425B', '#3A7CA5', '#81C3D7']
for file_name, color in zip(file_names, colors):
    # Read data from file
    data = np.loadtxt(file_name)

    # Extract third and fourth columns
    column_3 = data[:, 2]
    column_4 = data[:, 3]

    # Plot data with specified color
    plt.plot(column_3, column_4, label=file_name, color=color)

# Add labels and title
plt.xlabel('prey(x)')
plt.ylabel('predator(y)')
plt.title('Phase portrait with limit cycle of a Lotka-Volterra predator-prey model')

# Add legend and grid
plt.legend()
plt.grid(True)

# Show plot
plt.show()
