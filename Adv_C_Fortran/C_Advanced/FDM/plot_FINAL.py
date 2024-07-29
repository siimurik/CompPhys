import numpy as np
import matplotlib.pyplot as plt

# Load data from CSV files
def load_csv(filename):
    return np.genfromtxt(filename, delimiter=',', skip_header=0)

# Load temperature data
temperature_data = load_csv("temps.csv")

# Load x and y vector data
xy_data = load_csv("xy.csv")

# Check if temperature_data is loaded correctly
print("Temperature data shape:", temperature_data.shape)
print("XY data shape:", xy_data.shape)

# Extract grid dimensions
ny, nx = temperature_data.shape

# Extract x and y coordinates
x_coords = xy_data[:, 0]
y_coords = xy_data[:, 1]

# Create a 2D grid of x and y coordinates
X, Y = np.meshgrid(x_coords, y_coords)

# Temperature data should be reshaped into the grid shape
Z = temperature_data

# Plotting
fig, ax = plt.subplots(constrained_layout=True)

# Plot filled contour
CS = ax.contourf(X, Y, Z, 10, origin='lower')

# Plot contour lines
CS2 = ax.contour(CS, levels=CS.levels[::1], colors='k', origin='lower', linewidths=0.75)

# Add labels to contour lines
ax.clabel(CS2, fmt='%2.1f', colors='k', fontsize=10)

# Set titles and labels
ax.set_title('Temperature Distribution')
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')

# Add colorbar
cbar = fig.colorbar(CS)
cbar.ax.set_ylabel('Temperature (K)')

# Add contour line levels to the colorbar
cbar.add_lines(CS2)

# Show plot
plt.show()
