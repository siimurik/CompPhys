import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Read the CSV data using numpy
data = np.genfromtxt('output.csv', delimiter=',')

# Extract data columns
t = data[:, 0]
x = data[:, 1]
y = data[:, 2]

# Create a figure and axis
fig, ax = plt.subplots()
ax.set_xlabel('Position X')
ax.set_ylabel('Position Y')
ax.set_title('Projectile Trajectory')

# Plot the trajectory
ax.plot(x, y, label='Trajectory')

# Create a point (ball) that will move along the trajectory
ball, = ax.plot([], [], 'o', color='red')

# Define the initialization function
def init():
    ball.set_data([], [])
    return ball,

# Define the animation function
def animate(i):
    ball.set_data([x[i]], [y[i]])  # Pass lists with a single element
    return ball,

# Adjust the speed of the animation by changing the interval (milliseconds)
# For example, interval=20 for faster animation, interval=100 for slower animation
animation_speed = 20  # Adjust this value to change the speed

# Create the animation
ani = FuncAnimation(fig, animate, init_func=init, frames=len(t), interval=animation_speed, blit=True)

# Display the plot with animation
plt.legend()
plt.grid()
plt.show()
