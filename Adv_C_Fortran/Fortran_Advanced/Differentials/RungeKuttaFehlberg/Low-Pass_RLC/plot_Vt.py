import numpy as np
import matplotlib.pyplot as plt

# Read data using genfromtxt, which handles missing/extra whitespace well
data = np.genfromtxt('data1.dat', skip_header=0)

# Extract second and third columns (indexing starts at 0)
t = data[:, 0]  # Time
V = data[:, 1]  # Voltage

plt.plot(t, V, label='Voltage')
plt.legend()
plt.grid(True)
plt.show()
