import numpy as np
import matplotlib.pyplot as plt

# Parameters
x_val = 1.0
limit = 2.5
res = 400

# Create grid
u = np.linspace(-limit, limit, res)
v = np.linspace(-limit, limit, res)
U, V = np.meshgrid(u, v)

# 1. The "Height" of the landscape (Real part of exponent)
# h(u,v) = Re(psi) = v^3/3 - u^2*v - x*v
Z = (V**3 / 3) - (U**2 * V) - (x_val * V)

# 2. The "Flow" direction (Negative Gradient of Height)
# This vector points in the direction of steepest descent.
# Vector F = - grad(h) = ( -dh/du , -dh/dv )
# dh/du = -2uv  =>  -dh/du = 2uv
# dh/dv = v^2 - u^2 - x  =>  -dh/dv = u^2 + x - v^2
FlowU = 2 * U * V
FlowV = U**2 + x_val - V**2

# Plotting
fig, ax = plt.subplots(figsize=(10, 8))

# A. Contour lines of the height (Topography)
# Use a diverging colormap to show valleys (blue) and ridges (red)
cp = ax.contourf(U, V, Z, levels=np.linspace(-8, 8, 40), cmap='RdBu_r', alpha=0.5)
cbar = fig.colorbar(cp, ax=ax, label='Re($\psi(t)$) - Height (Red=Ridge, Blue=Valley)')
# Add thin lines for better definition
ax.contour(U, V, Z, levels=np.linspace(-8, 8, 40), colors='gray', linewidths=0.5, alpha=0.5)


# B. Streamplot to show the flow with arrows
# We use a lower resolution grid for the streamplot to avoid overcrowding arrows
res_stream = 50
u_s = np.linspace(-limit, limit, res_stream)
v_s = np.linspace(-limit, limit, res_stream)
U_s, V_s = np.meshgrid(u_s, v_s)
FlowU_s = 2 * U_s * V_s
FlowV_s = U_s**2 + x_val - V_s**2

st = ax.streamplot(U_s, V_s, FlowU_s, FlowV_s, color='black', linewidth=1, arrowsize=1.2, density=1.5)


# C. Mark key features
# Saddle points
s1 = (0, np.sqrt(x_val))
s2 = (0, -np.sqrt(x_val))
ax.plot(s1[0], s1[1], 'go', markersize=10, label='Contributing Saddle $t_1$', zorder=5)
ax.plot(s2[0], s2[1], 'rx', markersize=10, markeredgewidth=2, label='Non-contributing Saddle $t_2$', zorder=5)

# Original integration path
ax.axhline(0, color='blue', linestyle='--', linewidth=2, label='Original Path (Real Axis)')

# Title and labels
ax.set_title(f'Steepest Descent "Flow" Map for Airy Function (x={x_val})', fontsize=14)
ax.set_xlabel('Real(t) [u]', fontsize=12)
ax.set_ylabel('Imag(t) [v]', fontsize=12)
ax.legend(loc='upper right', framealpha=0.9)
ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)

plt.show()