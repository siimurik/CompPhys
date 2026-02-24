import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters
x_val = 1.0
u_range = np.linspace(-2.5, 2.5, 100)
v_range = np.linspace(-2.5, 2.5, 100)
U, V = np.meshgrid(u_range, v_range)

# Real part of the exponent Re(psi(t)) = v^3/3 - u^2*v - x*v
# This determines the height (magnitude) of the integrand exp(psi(t))
Z = (V**3 / 3) - (U**2 * V) - (x_val * V)

# Saddle points at (0, sqrt(x)) and (0, -sqrt(x))
s1_u, s1_v = 0, np.sqrt(x_val)
s1_z = (s1_v**3 / 3) - (s1_u**2 * s1_v) - (x_val * s1_v)

s2_u, s2_v = 0, -np.sqrt(x_val)
s2_z = (s2_v**3 / 3) - (s2_u**2 * s2_v) - (x_val * s2_v)

# Plotting
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Surface plot
# Limit Z values for better visualization (avoid huge peaks overshadowing everything)
Z_clipped = np.clip(Z, -10, 10)
surf = ax.plot_surface(U, V, Z_clipped, cmap='RdBu_r', alpha=0.8, edgecolor='none')

# Saddle points
ax.scatter([s1_u], [s1_v], [s1_z], color='black', s=100, label=f'Saddle t1 (ix^{1/2})', zorder=10)
ax.scatter([s2_u], [s2_v], [s2_z], color='purple', s=100, marker='x', label=f'Saddle t2 (-ix^{1/2})', zorder=10)

# Labels
ax.set_title(f'3D Surface of Re(psi(t)) for Airy Function (x={x_val})', fontsize=14)
ax.set_xlabel('Real(t) [u]', fontsize=10)
ax.set_ylabel('Imag(t) [v]', fontsize=10)
ax.set_zlabel('Re(psi(t)) [Height]', fontsize=10)

# Indicate the original path (v=0)
ax.plot(u_range, np.zeros_like(u_range), np.zeros_like(u_range), color='black', linewidth=3, linestyle='--', label='Original Real Axis Path')

# Indicate the steepest descent path (horizontal through t1)
path_u = np.linspace(-1.5, 1.5, 50)
path_v = np.full_like(path_u, s1_v)
path_z = (path_v**3 / 3) - (path_u**2 * path_v) - (x_val * path_v)
ax.plot(path_u, path_v, path_z, color='green', linewidth=4, label='Steepest Descent Path')

ax.legend()
ax.view_init(elev=30, azim=45)

plt.savefig('airy_3d_topography.png')
plt.show()