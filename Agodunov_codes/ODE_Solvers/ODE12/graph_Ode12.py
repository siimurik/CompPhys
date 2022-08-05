import pandas as pd
import matplotlib.pyplot as plt

# Setting background for graph
plt.style.use('dark_background')
#img = plt.imread("background.png")

# Reading in the general results about PE, KE, Te...
df = pd.read_table('result12.dat', skiprows=3, delim_whitespace=True,
                         names=['t', 'x(t)', 'v(t)', 'energy']      )

# Printing the values from dataframe
print(df.head())

# Assigning names to datacolumns
t = df['t']
x = df['x(t)']
v = df['v(t)']
E = df['energy']

# Plotting x(t), v(t) and E
fig = plt.figure()
ax = fig.add_subplot(111)
#ax.imshow(img, extent = [-0.35, 7.35, -3.35, 4.8])
ax.plot(t, x, color='red',     label="Position")
ax.plot(t, v, color='blue',    label="Velocity")
ax.plot(t, E, color='magenta', label="Energy")
ax.set_xlabel('Time')
ax.set_ylabel('Position, velocity, energy')
#ax.set_xlim(-0.35, 7.35)
#ax.set_ylim(-3.35, 4.80)
plt.title("Harmonic oscillator")
plt.grid()
plt.legend()
plt.show()