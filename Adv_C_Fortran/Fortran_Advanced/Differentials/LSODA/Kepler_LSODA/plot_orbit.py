#!/usr/bin/env python3
"""
Orbital Trajectory Plotter
Reads orbit data from DLSODA simulation and plots trajectory with Earth
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.animation as animation

def plot_static_orbit(filename='orbit.dat'):
    """Plot the complete orbital trajectory as a static plot"""
    
    try:
        # Read data from file
        data = np.loadtxt(filename)
        t = data[:, 0]    # time (s)
        x = data[:, 1]    # x position (km)
        y = data[:, 2]    # y position (km)
        vx = data[:, 3]   # x velocity (km/s)
        vy = data[:, 4]   # y velocity (km/s)
        r = data[:, 5]    # distance from Earth center (km)
        
        print(f"Loaded {len(t)} data points")
        print(f"Time range: {t[0]:.1f} to {t[-1]:.1f} seconds ({t[-1]/3600:.2f} hours)")
        print(f"Distance range: {r.min():.1f} to {r.max():.1f} km")
        print(f"Altitude range: {r.min()-6371:.1f} to {r.max()-6371:.1f} km")
        
    except FileNotFoundError:
        print(f"Error: Cannot find file '{filename}'")
        print("Make sure to run the DLSODA orbital simulation first!")
        return
    except Exception as e:
        print(f"Error reading data: {e}")
        return
    
    # Create figure with subplots
    fig = plt.figure(figsize=(15, 10))
    
    # Main orbital plot
    ax1 = plt.subplot(2, 2, (1, 3))
    
    # Draw Earth
    earth_radius = 6371.0  # km
    earth = Circle((0, 0), earth_radius, color='lightblue', alpha=0.7, label='Earth')
    ax1.add_patch(earth)
    
    # Draw Earth's surface
    earth_surface = Circle((0, 0), earth_radius, fill=False, color='blue', linewidth=2)
    ax1.add_patch(earth_surface)
    
    # Plot trajectory
    ax1.plot(x, y, 'r-', linewidth=1.5, alpha=0.8, label='Orbital trajectory')
    
    # Mark start and end points
    ax1.plot(x[0], y[0], 'go', markersize=8, label='Start')
    ax1.plot(x[-1], y[-1], 'ro', markersize=8, label='End')
    
    # Add some trajectory markers every 1000 points
    step = max(1, len(x) // 20)
    ax1.plot(x[::step], y[::step], 'k.', markersize=3, alpha=0.5)
    
    ax1.set_xlabel('X Position (km)')
    ax1.set_ylabel('Y Position (km)')
    ax1.set_title('Orbital Trajectory Around Earth')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.axis('equal')
    
    # Set axis limits to show orbit clearly
    max_range = max(abs(x).max(), abs(y).max()) * 1.1
    ax1.set_xlim(-max_range, max_range)
    ax1.set_ylim(-max_range, max_range)
    
    # Distance vs time plot
    ax2 = plt.subplot(2, 2, 2)
    ax2.plot(t/3600, r, 'b-', linewidth=1.5)
    ax2.axhline(y=earth_radius, color='r', linestyle='--', alpha=0.7, label='Earth surface')
    ax2.set_xlabel('Time (hours)')
    ax2.set_ylabel('Distance from Earth center (km)')
    ax2.set_title('Orbital Distance vs Time')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Velocity vs time plot
    ax3 = plt.subplot(2, 2, 4)
    v_magnitude = np.sqrt(vx**2 + vy**2)
    ax3.plot(t/3600, v_magnitude, 'g-', linewidth=1.5, label='Speed')
    ax3.plot(t/3600, vx, 'r-', alpha=0.7, label='Vx')
    ax3.plot(t/3600, vy, 'b-', alpha=0.7, label='Vy')
    ax3.set_xlabel('Time (hours)')
    ax3.set_ylabel('Velocity (km/s)')
    ax3.set_title('Orbital Velocity vs Time')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    plt.tight_layout()
    plt.show()
    
    # Print orbital statistics
    print("\nOrbital Statistics:")
    print(f"Perigee (closest): {r.min():.1f} km ({r.min()-6371:.1f} km altitude)")
    print(f"Apogee (farthest): {r.max():.1f} km ({r.max()-6371:.1f} km altitude)")
    print(f"Speed range: {v_magnitude.min():.3f} to {v_magnitude.max():.3f} km/s")
    print(f"Orbital period estimate: {t[-1]/3600:.2f} hours (if completing full orbit)")

def create_animation(filename='orbit.dat', interval=50, save_gif=False):
    """Create an animated plot of the orbital trajectory"""
    
    try:
        # Read data
        data = np.loadtxt(filename)
        x = data[:, 1]
        y = data[:, 2]
        
        # Reduce data points for smoother animation
        step = max(1, len(x) // 1000)  # Max 1000 frames
        x = x[::step]
        y = y[::step]
        
    except FileNotFoundError:
        print(f"Error: Cannot find file '{filename}'")
        return
    
    # Set up the figure and axis
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Draw Earth
    earth_radius = 6371.0
    earth = Circle((0, 0), earth_radius, color='lightblue', alpha=0.7)
    ax.add_patch(earth)
    
    # Initialize empty trajectory line and current position
    trajectory_line, = ax.plot([], [], 'r-', linewidth=1.5, alpha=0.7)
    current_pos, = ax.plot([], [], 'ro', markersize=8)
    
    # Set axis properties
    max_range = max(abs(x).max(), abs(y).max()) * 1.1
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.set_xlabel('X Position (km)')
    ax.set_ylabel('Y Position (km)')
    ax.set_title('Orbital Trajectory Animation')
    ax.grid(True, alpha=0.3)
    ax.axis('equal')
    
    def animate(frame):
        # Update trajectory up to current frame
        trajectory_line.set_data(x[:frame+1], y[:frame+1])
        # Update current position
        if frame < len(x):
            current_pos.set_data([x[frame]], [y[frame]])
        return trajectory_line, current_pos
    
    # Create animation
    anim = animation.FuncAnimation(fig, animate, frames=len(x), 
                                 interval=interval, blit=True, repeat=True)
    
    if save_gif:
        print("Saving animation as orbit_animation.gif...")
        anim.save('orbit_animation.gif', writer='pillow', fps=20)
        print("Animation saved!")
    
    plt.show()
    return anim

def plot_3d_trajectory(filename='orbit.dat'):
    """Create a 3D visualization of the orbit (assuming z=0)"""
    
    try:
        from mpl_toolkits.mplot3d import Axes3D
        
        # Read data
        data = np.loadtxt(filename)
        x = data[:, 1]
        y = data[:, 2]
        z = np.zeros_like(x)  # Assume planar orbit
        
    except ImportError:
        print("3D plotting requires mpl_toolkits.mplot3d")
        return
    except FileNotFoundError:
        print(f"Error: Cannot find file '{filename}'")
        return
    
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Draw Earth as a sphere
    earth_radius = 6371.0
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    earth_x = earth_radius * np.outer(np.cos(u), np.sin(v))
    earth_y = earth_radius * np.outer(np.sin(u), np.sin(v))
    earth_z = earth_radius * np.outer(np.ones(np.size(u)), np.cos(v))
    
    ax.plot_surface(earth_x, earth_y, earth_z, alpha=0.3, color='lightblue')
    
    # Plot trajectory
    ax.plot(x, y, z, 'r-', linewidth=2, label='Orbital trajectory')
    ax.scatter(x[0], y[0], z[0], color='green', s=100, label='Start')
    ax.scatter(x[-1], y[-1], z[-1], color='red', s=100, label='End')
    
    ax.set_xlabel('X Position (km)')
    ax.set_ylabel('Y Position (km)')
    ax.set_zlabel('Z Position (km)')
    ax.set_title('3D Orbital Trajectory')
    ax.legend()
    
    plt.show()

if __name__ == "__main__":
    print("Orbital Trajectory Plotter")
    print("=" * 30)
    
    # Check if data file exists
    import os
    if not os.path.exists('orbit.dat'):
        print("No orbit.dat file found!")
        print("Please run the DLSODA orbital simulation first to generate data.")
        exit(1)
    
    print("Choose plotting option:")
    print("1. Static orbital plot (recommended)")
    print("2. Animated trajectory")
    print("3. 3D visualization")
    print("4. All plots")
    
    try:
        choice = input("Enter choice (1-4): ").strip()
        
        if choice == '1' or choice == '4':
            print("\nGenerating static orbital plot...")
            plot_static_orbit()
        
        if choice == '2' or choice == '4':
            print("\nGenerating animated trajectory...")
            anim = create_animation()
        
        if choice == '3' or choice == '4':
            print("\nGenerating 3D visualization...")
            plot_3d_trajectory()
            
    except KeyboardInterrupt:
        print("\nExiting...")
    except Exception as e:
        print(f"Error: {e}")
        # Default to static plot
        print("Showing static plot...")
        plot_static_orbit()