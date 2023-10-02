import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Set the backend (to move the windows to the desired location on the screen)
matplotlib.use("TkAgg")

# Input parameters
nTrials = int(1E4)
radius = 1

# Counter for the pins inside the circle
nInside = 0
# Counter for the pins dropped
nDrops = 0

# Generate points in a square of side 2 units, from -1 to 1.
XrandCoords = np.random.default_rng().uniform(-1, 1, (nTrials,))
YrandCoords = np.random.default_rng().uniform(-1, 1, (nTrials,))

# Create a single figure for both plots side by side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Set limits and labels for the first plot
ax1.set_xlim(-1, 1)
ax1.set_ylim(-1, 1)
ax1.set_title('Simulation')

# Set limits and labels for the second plot
ax2.set_ylim(2.8, 4.3)
ax2.set_xlim(0, nTrials)
ax2.set_title('π estimate vs no. of pin drops')
ax2.axhline(y=np.pi, c='darkseagreen', linewidth=2, alpha=0.5)
ax2.annotate('π', [0, np.pi], fontsize=20)
ax2.grid()

# Some arrays to store the pi value vs the number of pins dropped
piValueI = []
nDrops_arr = []

# Some arrays to plot the points
insideX = []
outsideX = []
insideY = []
outsideY = []

# Begin Monte Carlo
for i in range(nTrials):
    x = XrandCoords[i]
    y = YrandCoords[i]
    # Increment the counter for the number of total pins dropped
    nDrops = nDrops + 1

    # Check if the points are inside the circle or not
    if x ** 2 + y ** 2 <= radius ** 2:
        nInside = nInside + 1
        insideX.append(x)
        insideY.append(y)
    else:
        outsideX.append(x)
        outsideY.append(y)

    # plot only at some values
    if i % 100 == 0:
        # Draw on the first window
        if i == 0:
            ax1.scatter(insideX, insideY, c='pink', s=50, label='Drop inside')
            ax1.scatter(outsideX, outsideY, c='orange', s=50, label='Drop outside')
            ax1.legend(loc=(0.75, 0.9))
        else:
            ax1.scatter(insideX, insideY, c='pink', s=50)
            ax1.scatter(outsideX, outsideY, c='orange', s=50)

        area = 4 * nInside / nDrops
        ax1.set_title('No. of pin drops = ' + str(nDrops) +
                      '; No. inside circle = ' + str(nInside) +
                      r'; π ≈ $4\frac{N_\mathrm{inside}}{N_\mathrm{total}}=$ ' + str(np.round(area, 6)))
        piValueI.append(area)
        nDrops_arr.append(nDrops)

        # Update the plot on the second window
        ax2.plot(nDrops_arr, piValueI, c='mediumorchid')
        plt.pause(0.01)

area = 4 * nInside / nTrials
print("Final estimated value of Pi:", area)

plt.show()
