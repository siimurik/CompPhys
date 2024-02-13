import numpy as np
import matplotlib.pyplot as plt

# Define a function to read the data from the file
def read_data(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        data = [[float(x) for x in line.split()] for line in lines]
    return data

# List of file names
file_names = ['rk_50.dat', 'rk_500.dat', 'rk_5000.dat', 'rk_50000.dat']

# Plotting
for filename in file_names:
    data = read_data(filename)
    data = list(zip(*data))  # Transpose the data
    t = data[0]
    x = data[1]

    plt.plot(t, x, label=filename)

#Analytical solution to Simple Harmonic oscillator
t = np.linspace(0, 6, 50000)
f = lambda x: 0.2*np.cos(np.sqrt(10)*x) #+ v0/omega*sin(omega*x)

plt.plot(t, f(t), label = "Anal. sol." )
plt.xlabel('Time (t)')
plt.ylabel('Position (x)')
plt.title('Position vs Time')
plt.legend()
plt.grid(True)
plt.show()
