import pandas as pd
import matplotlib.pyplot as plt

# Setting background for graph
plt.style.use('dark_background')
#img = plt.imread("background.png")

# Reading in the general results about PE, KE, Te...
df = pd.read_table('shoot2.dat', skiprows=1, 
        delim_whitespace=True, names=['x', 'y', 'dy'] )

# Printing the values from dataframe
print(df.head())

# Assigning names to datacolumns
x = df['x']
y = df['y']
dy = df['dy']

# Plotting x(t), v(t) and E
fig = plt.figure()
ax = fig.add_subplot(111)
#ax.imshow(img, extent = [-0.35, 7.35, -3.35, 4.8])
ax.plot(x, y,   label="Position")
#ax.plot(y, dy,  label="Velocity")
#ax.plot(x, dy, label="Energy")
ax.set_xlabel('Position')
ax.set_ylabel('Position')
#ax.set_xlim(-0.35, 7.35)
#ax.set_ylim(-3.35, 4.80)
plt.title("Shooting method")
plt.grid()
plt.legend()
plt.show()