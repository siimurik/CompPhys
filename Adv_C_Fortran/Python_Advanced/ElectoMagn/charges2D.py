import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def calculate_force(q1, q2, r1, r2):
    """
    Calculate the Coulomb force between two charges.

    Args:
        q1, q2: Charges of the particles.
        r1, r2: Positions of the particles (numpy arrays).

    Returns:
        Force on q1 due to q2 (numpy array).
    """
    k = 8.99e9  # Coulomb constant (N*m^2/C^2)
    r = r2 - r1
    distance = np.linalg.norm(r)
    if distance == 0:
        return np.array([0.0, 0.0])  # Avoid division by zero
    force_magnitude = k * q1 * q2 / distance**2
    force_direction = r / distance
    return force_magnitude * force_direction

def update_positions(r, v, a, dt):
    """
    Update positions and velocities using simple Euler integration.

    Args:
        r: Position of the particle (numpy array).
        v: Velocity of the particle (numpy array).
        a: Acceleration of the particle (numpy array).
        dt: Time step.

    Returns:
        Updated position and velocity.
    """
    v += a * dt
    r += v * dt
    return r, v

def animate(i, dt, r1, v1, q1, r2, v2, q2, line1, line2, trail1, trail2):
    """
    Update function for the animation.
    """
    global r1_global, v1_global, r2_global, v2_global

    if i == 0:
        # Initialize global positions and velocities with the initial values
        r1_global, v1_global = r1, v1
        r2_global, v2_global = r2, v2

    force_on_1 = calculate_force(q1, q2, r1_global, r2_global)
    force_on_2 = -force_on_1  # Newton's third law

    # Acceleration (F = ma, assume unit mass for simplicity)
    a1 = -force_on_1
    a2 = -force_on_2

    # Update positions and velocities
    r1_global, v1_global = update_positions(r1_global, v1_global, a1, dt)
    r2_global, v2_global = update_positions(r2_global, v2_global, a2, dt)

    # Update the plot data
    line1.set_data([r1_global[0]], [r1_global[1]])
    line2.set_data([r2_global[0]], [r2_global[1]])

    # Update the trails
    trail1.set_data(np.append(trail1.get_xdata(), r1_global[0]), np.append(trail1.get_ydata(), r1_global[1]))
    trail2.set_data(np.append(trail2.get_xdata(), r2_global[0]), np.append(trail2.get_ydata(), r2_global[1]))

    return line1, line2, trail1, trail2

# Initial conditions
q1 = float(input("Enter charge of particle 1 (in C): "))
q2 = float(input("Enter charge of particle 2 (in C): "))
r1_global = np.array([3.5, 5.0])
v1_global = np.array([0.0, -0.1])  # Moving towards negative y direction
r2_global = np.array([4.5, 3.0])
v2_global = np.array([0.0, 0.1])   # Moving towards positive y direction

# Animation setup
dt = 0.1
fig, ax = plt.subplots()
ax.set_xlim(0, 8)
ax.set_ylim(0, 8)
ax.set_aspect('equal')
ax.set_title("Charged Particle Motion")
ax.set_xlabel("x-axis")
ax.set_ylabel("y-axis")
line1, = ax.plot([], [], 'ro', label='Particle 1')
line2, = ax.plot([], [], 'bo', label='Particle 2')
trail1, = ax.plot([], [], 'r-', alpha=0.5)  # Trail for particle 1
trail2, = ax.plot([], [], 'b-', alpha=0.5)  # Trail for particle 2
ax.legend()

# Animate
ani = FuncAnimation(fig, animate, frames=200,
                    fargs=(dt, r1_global, v1_global, q1,
                           r2_global, v2_global, q2,
                           line1, line2,
                           trail1, trail2),
                    interval=50)

plt.grid()
plt.show()