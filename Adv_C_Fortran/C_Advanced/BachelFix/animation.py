import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# Define the rotation function
def sequentialRotations(x, y, z, alpha, beta, gamma, delta, iota):
    """
    Applies a series of rotations to a 3D vector using specified angles around the XYZ axes.
    """
    # Rotation around z-axis (alpha)
    a1 = np.array([[np.cos(alpha), -np.sin(alpha), 0], [np.sin(alpha), np.cos(alpha), 0], [0, 0, 1]])
    a = a1 @ np.array([[x], [y], [z]])
    # Rotation around y-axis (beta)
    b2 = np.array([[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [-np.sin(beta), 0, np.cos(beta)]])
    b = b2 @ a
    # Rotation around x-axis (gamma)
    g1 = np.array([[1, 0, 0], [0, np.cos(gamma), -np.sin(gamma)], [0, np.sin(gamma), np.cos(gamma)]])
    g = g1 @ b
    # Rotation around y-axis (delta)
    d2 = np.array([[np.cos(delta), 0, np.sin(delta)], [0, 1, 0], [-np.sin(delta), 0, np.cos(delta)]])
    d = d2 @ g
    # Rotation around z-axis (iota)
    i1 = np.array([[np.cos(iota), -np.sin(iota), 0], [np.sin(iota), np.cos(iota), 0], [0, 0, 1]])
    return (i1 @ d).flatten()

# Create an initial 3D vector (random unit vector)
np.random.seed(0)
initial_vector = np.random.rand(3)
initial_vector /= np.linalg.norm(initial_vector)

# Define the rotation angles (random values)
alpha = np.deg2rad(30)
beta = np.deg2rad(45)
gamma = np.deg2rad(60)
delta = np.deg2rad(30)
iota = np.deg2rad(90)

# Prepare the figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Function to update the animation
def update(frame):
    ax.clear()
    angles = [alpha, beta, gamma, delta, iota]
    vector = initial_vector
    for i in range(frame):
        angle = angles[i % len(angles)]
        if i % len(angles) == 0:
            matrix = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        elif i % len(angles) == 1:
            matrix = np.array([[np.cos(angle), 0, np.sin(angle)], [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]])
        elif i % len(angles) == 2:
            matrix = np.array([[1, 0, 0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]])
        elif i % len(angles) == 3:
            matrix = np.array([[np.cos(angle), 0, np.sin(angle)], [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]])
        elif i % len(angles) == 4:
            matrix = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        vector = matrix @ vector
    
    ax.quiver(0, 0, 0, vector[0], vector[1], vector[2], color='b')
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f"Rotation Step: {frame}")

# Create the animation
ani = FuncAnimation(fig, update, frames=30, interval=500, repeat=False)

# Display the animation
plt.show()
