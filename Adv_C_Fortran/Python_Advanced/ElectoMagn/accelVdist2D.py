import numpy as np
import matplotlib.pyplot as plt

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

# Initial conditions
q1 = 1E-6
q2 = 1E-6
r1 = np.array([3.5, 5.0])
v1 = np.array([0.0, -0.2])  # Moving towards negative y direction
r2 = np.array([4.5, 3.0])
v2 = np.array([0.0, 0.3])   # Moving towards positive y direction

# Time step and total time
dt = 0.1
total_time = 20  # seconds
num_steps = int(total_time / dt)

# Lists to store distance and acceleration values
distances = []
accelerations_1 = []
accelerations_2 = []

for _ in range(num_steps):
    force_on_1 = calculate_force(q1, q2, r1, r2)
    force_on_2 = -force_on_1  # Newton's third law

    # Acceleration (F = ma, assume unit mass for simplicity)
    a1 = -force_on_1
    a2 = -force_on_2

    # Update positions and velocities
    r1, v1 = update_positions(r1, v1, a1, dt)
    r2, v2 = update_positions(r2, v2, a2, dt)

    # Calculate and store distance and acceleration
    distance = np.linalg.norm(r1 - r2)
    acceleration_1 = np.linalg.norm(a1)
    acceleration_2 = np.linalg.norm(a2)
    
    distances.append(distance)
    accelerations_1.append(acceleration_1)
    accelerations_2.append(acceleration_2)

# Plotting the acceleration vs. distance
plt.figure(figsize=(8, 6))
plt.plot(distances, accelerations_1, 'r*-', label='Particle 1')
#plt.plot(distances, accelerations_2, 'bo-', label='Particle 2')

# Marking the start and end positions
plt.scatter(distances[0], accelerations_1[0], color='green', s=100, label='Start')
plt.scatter(distances[-1], accelerations_1[-1], color='purple', s=100, label='End')

plt.title("Acceleration vs. Distance")
plt.xlabel("Distance (m)")
plt.ylabel("Acceleration (m/s^2)")
plt.legend()
plt.grid(True)
plt.show()
