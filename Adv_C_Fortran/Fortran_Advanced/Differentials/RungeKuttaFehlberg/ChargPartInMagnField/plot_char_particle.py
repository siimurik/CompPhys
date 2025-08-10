import numpy as np
import matplotlib.pyplot as plt

def read_trajectory(filename):
    """
    Reads the trajectory file of the particle.
    The file should have columns: t, x, y, z, vx, vy, vz
    Returns arrays: t, x, y, z, vx, vy, vz
    """
    data = np.loadtxt(filename, skiprows=1)
    t = data[:,0]
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]
    vx = data[:,4]
    vy = data[:,5]
    vz = data[:,6]
    return t, x, y, z, vx, vy, vz

def plot_earth(ax, radius=6.378137e6, color='lightblue', alpha=0.5):
    """
    Plots a sphere representing the Earth on the provided Axes3D object.
    """
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    xs = radius * np.outer(np.cos(u), np.sin(v))
    ys = radius * np.outer(np.sin(u), np.sin(v))
    zs = radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(xs, ys, zs, color=color, alpha=alpha, linewidth=0)

def plot_trajectory(filename, earth_radius=6.378137e6):
    """
    Plots the Earth and the particle's trajectory in 3D.
    """
    t, x, y, z, vx, vy, vz = read_trajectory(filename)
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot Earth's surface
    plot_earth(ax, radius=earth_radius)

    # Plot trajectory
    ax.plot3D(x, y, z, color='magenta', label='Particle trajectory')

    # Optionally: plot starting point
    ax.scatter([x[0]], [y[0]], [z[0]], color='black', s=50, label='Start')

    # Format plot
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')
    ax.set_title('Charged Particle Trajectory around Earth')
    ax.legend()
    ax.set_box_aspect([1,1,1])

    # Set limits to make Earth visible and trajectory clear
    lim = np.max(np.abs([x, y, z]))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python plot_trajectory.py timesteps_output.txt")
        sys.exit(1)
    plot_trajectory(sys.argv[1])