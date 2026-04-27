import numpy as np
import matplotlib.pyplot as plt

# Define the exponent function psi(t) = i(t^3 / 3 + x*t)
def psi(t, x):
    return 1j * (t**3 / 3 + x * t)

# Set parameters
x_val = 1.0
limit = 2.5
res = 400

# Create grid
u = np.linspace(-limit, limit, res)
v = np.linspace(-limit, limit, res)
U, V = np.meshgrid(u, v)
T = U + 1j * V

# Calculate the Real part of the exponent (this determines the magnitude)
Z = np.real(psi(T, x_val))

# Plotting
plt.figure(figsize=(10, 8))
# Use a diverging color map: blue for valleys (negative), red for ridges (positive)
cp = plt.contourf(U, V, Z, levels=np.linspace(-5, 5, 50), cmap='RdBu_r', extend='both')
plt.colorbar(cp, label='Re(psi(t)) - Magnitude level')

# Mark the saddle points
saddle1 = 1j * np.sqrt(x_val)
saddle2 = -1j * np.sqrt(x_val)
plt.plot(saddle1.real, saddle1.imag, 'ko', markersize=8, label='Saddle Point t1 (Contributes)')
plt.plot(saddle2.real, saddle2.imag, 'kx', markersize=8, label='Saddle Point t2 (Ignored)')

# Draw the original path (Real axis)
plt.axhline(0, color='black', linestyle='--', alpha=0.5, label='Original Path (Real Axis)')

# Draw the deformed path (Steepest Descent through t1)
# For Ai(x), the steepest descent path through t1 is roughly a horizontal line 
# that bends down into the valleys.
path_u = np.linspace(-2, 2, 100)
path_v = np.sqrt(x_val) + 0.1 * path_u**2 # Simple approximation of the curve
# Actually, near the saddle it is purely horizontal. 
plt.plot(path_u, np.full_like(path_u, np.sqrt(x_val)), 'g-', linewidth=2, label='Deformed Path (Steepest Descent)')

# Annotate sectors
plt.text(1.5, 1.0, 'Valley 1', fontsize=12, fontweight='bold')
plt.text(-2.0, 1.0, 'Valley 2', fontsize=12, fontweight='bold')
plt.text(-0.5, -2.0, 'Valley 3', fontsize=12, fontweight='bold')
plt.text(0, 2.0, 'RIDGE', color='white', ha='center', fontweight='bold')

plt.title(f'Topography of the Airy Function Exponent (x={x_val})')
plt.xlabel('Real(t)')
plt.ylabel('Imag(t)')
plt.legend(loc='lower right')
plt.grid(alpha=0.3)
plt.show()