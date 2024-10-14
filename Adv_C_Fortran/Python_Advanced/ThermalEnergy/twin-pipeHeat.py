import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from numba import njit

# Parameters
Lx, Ly = 10.0, 10.0  # Domain size (meters)
Nx, Ny = 100, 100    # Grid points for better resolution
dx, dy = Lx / Nx, Ly / Ny  # Spatial step sizes
alpha = 1e-3  # Increased thermal diffusivity for faster heat spread (m^2/s)
dt = 0.1      # Larger time step size (s) to make the changes more noticeable //  Courant-Friedrichs-Lewy (CFL) condition for stability dt <= dx**2/(4*alpha) 
Nt = 20000     # Number of time steps   

# The time to reach a steady state is approximately governed by the thermal 
# diffusion equation. A rough estimate can be made by the equation:
#t_steady = Lx**2/alpha
#NumOFTimeSteps = t_steady/dt
#print("Number of time steps is: {:d}".format(int(NumOFTimeSteps)))

# Pipe properties
pipe_radius             = 0.35  # Radius of the pipes (m)
insulation_thickness    = 0.5 # Insulation thickness (m)
pipe_temperature_hot    = 80.0  # Temperature of the hot pipe (°C)
pipe_temperature_cold   = 45.0  # Temperature of the cold pipe (°C)
pipe_pos_hot            = (int(Nx / 4), int(Ny / 2))  # Position of the hot pipe
pipe_pos_cold           = (int(3 * Nx / 4), int(Ny / 2))  # Position of the cold pipe
ground_conductivity     = 1.5  # W/m·K, for moist soil
insulation_conductivity = 0.035  # W/m·K, for fiberglass insulation

# Initial and boundary conditions
T_ground = 10.0 # (°C)
T        = np.full((Nx, Ny), T_ground)  # Initial soil temperature (°C)

# Function to apply pipe temperatures at the given positions
@njit
def apply_pipe_temperatures(T, pipe_pos_hot, pipe_pos_cold, pipe_radius, dx, dy, pipe_temperature_hot, pipe_temperature_cold):
    for i in range(Nx):
        for j in range(Ny):
            # Apply hot pipe temperature
            distance_to_hot_pipe = np.sqrt(((i - pipe_pos_hot[0]) * dx)**2 + ((j - pipe_pos_hot[1]) * dy)**2)
            if distance_to_hot_pipe <= pipe_radius:
                T[i, j] = pipe_temperature_hot
            # Apply cold pipe temperature
            distance_to_cold_pipe = np.sqrt(((i - pipe_pos_cold[0]) * dx)**2 + ((j - pipe_pos_cold[1]) * dy)**2)
            if distance_to_cold_pipe <= pipe_radius:
                T[i, j] = pipe_temperature_cold
    return T


# Function to update the temperature distribution (finite difference method)
@njit
def update_temperature(T, alpha, dx, dy, dt, Nx, Ny, T_ground):
    T_new = T.copy()
    
    # Finite difference update for interior points
    for i in range(1, Nx-1):                                    # If you try to compute the temperature at the boundary  
        for j in range(1, Ny-1):                                # points (i.e., i=0 or i=Nx-1), you will run into a  
            T_new[i, j] = T[i, j] + alpha * dt * (              # problem: the finite difference approximation requires 
                (T[i+1, j] - 2*T[i, j] + T[i-1, j]) / dx**2 +   # access to points outside the  grid (e.g., i=-1 or i=Nx,
                (T[i, j+1] - 2*T[i, j] + T[i, j-1]) / dy**2     # which do not exist).
            )
    
    # Apply boundary conditions
    # Dirichlet (fixed temperature at left, right and bottom)
    T_new[:, 0] = T_ground      # Left boundary
    T_new[:, -1] = T_ground     # Right boundary
    T_new[-1, :] = T_ground     # Bottom boundary

    # Neumann (no heat flux at top)
    T_new[0, :] = T_new[1, :]   # Top boundary
    #T_new[-1, :] = T_new[-2, :] # Bottom boundary
    
    return T_new

@njit
def apply_pipe_temperatures_with_insulation(T, pipe_pos_hot, pipe_pos_cold, pipe_radius, insulation_thickness, dx, dy, pipe_temperature_hot, pipe_temperature_cold):
    for i in range(Nx):
        for j in range(Ny):
            # Calculate the distance to each pipe
            distance_to_hot_pipe = np.sqrt(((i - pipe_pos_hot[0]) * dx)**2 + ((j - pipe_pos_hot[1]) * dy)**2)
            distance_to_cold_pipe = np.sqrt(((i - pipe_pos_cold[0]) * dx)**2 + ((j - pipe_pos_cold[1]) * dy)**2)
            
            # Apply hot pipe temperature within the pipe radius
            if distance_to_hot_pipe <= pipe_radius:
                T[i, j] = pipe_temperature_hot
            # Insulation layer for hot pipe (apply different heat transfer in this region)
            elif pipe_radius < distance_to_hot_pipe <= pipe_radius + insulation_thickness:
                T[i, j] = T[i, j]  # Temperature is updated differently in this region, handled by insulation_conductivity
            
            # Apply cold pipe temperature within the pipe radius
            if distance_to_cold_pipe <= pipe_radius:
                T[i, j] = pipe_temperature_cold
            # Insulation layer for cold pipe
            elif pipe_radius < distance_to_cold_pipe <= pipe_radius + insulation_thickness:
                T[i, j] = T[i, j]  # Temperature is updated differently in this region, handled by insulation_conductivity
    return T

@njit
def update_temperature_with_insulation(T, alpha, dx, dy, dt, Nx, Ny, pipe_pos_hot, pipe_pos_cold, pipe_radius, insulation_thickness, ground_conductivity, insulation_conductivity, T_ground):
    T_new = T.copy()
    
    for i in range(1, Nx-1):
        for j in range(1, Ny-1):
            # Calculate distances to hot and cold pipes
            distance_to_hot_pipe = np.sqrt(((i - pipe_pos_hot[0]) * dx)**2 + ((j - pipe_pos_hot[1]) * dy)**2)
            distance_to_cold_pipe = np.sqrt(((i - pipe_pos_cold[0]) * dx)**2 + ((j - pipe_pos_cold[1]) * dy)**2)
            
            # Determine conductivity based on distance (whether in insulation or ground)
            if pipe_radius < distance_to_hot_pipe <= pipe_radius + insulation_thickness or \
               pipe_radius < distance_to_cold_pipe <= pipe_radius + insulation_thickness:
                conductivity = insulation_conductivity
            else:
                conductivity = ground_conductivity
            
            # Finite difference update with modified conductivity
            T_new[i, j] = T[i, j] + conductivity * alpha * dt * (
                (T[i+1, j] - 2*T[i, j] + T[i-1, j]) / dx**2 +
                (T[i, j+1] - 2*T[i, j] + T[i, j-1]) / dy**2
            )
    
    # Apply boundary conditions (Dirichlet for ground, Neumann for top/bottom)
    T_new[:, 0]  = T_ground     # Left boundary (Dirichlet)
    T_new[:, -1] = T_ground     # Right boundary (Dirichlet)
    T_new[0, :]  = T_new[1, :]  # Top boundary (Neumann)
    T_new[-1, :] = T_ground     # Bottom boundary (Dirichlet)
    
    return T_new

# Main simulation loop
plt.figure(figsize=(6, 6))  # Larger figure for better visibility
for n in range(1, Nt+1):
    # Apply pipe boundary conditions
    T = apply_pipe_temperatures(                T, pipe_pos_hot, pipe_pos_cold, pipe_radius,                       dx, dy, pipe_temperature_hot, pipe_temperature_cold)
    #T = apply_pipe_temperatures_with_insulation(T, pipe_pos_hot, pipe_pos_cold, pipe_radius, insulation_thickness, dx, dy, pipe_temperature_hot, pipe_temperature_cold)
    
    # Update the temperature distribution
    T = update_temperature(                T, alpha, dx, dy, dt, Nx, Ny, T_ground)
    #T = update_temperature_with_insulation(T, alpha, dx, dy, dt, Nx, Ny, pipe_pos_hot, pipe_pos_cold, pipe_radius, insulation_thickness, ground_conductivity, insulation_conductivity, T_ground)
    
    # Plot the heat map every 100 steps without blocking
    if n % 100 == 0 or n == Nt:  # Plot at regular intervals and at the final step
        plt.clf()  # Clear the current plot
        plt.imshow(T.T, origin='lower', extent=[0, Lx, 0, Ly], cmap=cm.plasma)
        plt.colorbar(label="Temperature (°C)")
        plt.title(f"Time step {n}")
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.pause(0.001)  # Pause to update the plot

# Keep the final plot
plt.show()

# Create contour plot after the main simulation
# NOTE: This part only works if the 'for' loop 
# for 'T' has already been done. 
x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(x, y)

fig1, ax1 = plt.subplots(constrained_layout=True)
CS = ax1.contourf(X, Y, T.T, 10, origin='lower', cmap=cm.plasma)  # Create filled contour plot
CS2 = ax1.contour(CS, levels=CS.levels[::1], colors='k', origin='lower', linewidths=0.75)  # Add contour lines
ax1.clabel(CS2, fmt='%2.1f', colors='k', fontsize=10)  # Label contour lines

ax1.set_title('Temperature Distribution Contour')
ax1.set_xlabel('X (m)')
ax1.set_ylabel('Y (m)')

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel('Temperature (°C)')
# Add the contour line levels to the colorbar
cbar.add_lines(CS2)

plt.show()
