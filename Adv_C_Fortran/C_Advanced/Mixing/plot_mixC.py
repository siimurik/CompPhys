import numpy as np
import matplotlib.pyplot as plt

# Replace 'your_file.csv' with the actual path to your CSV file
file_path = 'colsRevamp.csv'
#file_path = 'coloursF.csv'


# Use numpy.genfromtxt to load CSV values into a NumPy matrix
dataRaw = np.genfromtxt(file_path, delimiter=',', filling_values=np.nan)    # shape 512x513
data = dataRaw[:, :len(dataRaw)]  # Removing last column bc it is full of 'nan' values; 
print(data.shape)
print(data)
N = data.shape[0]

# Some magic numbers for the plot to display correctly
px = 1 / 72  # pixel in inches
# Create a figure with one row and two columns
fig, axs = plt.subplots(1, 2, figsize=(1320 * px, 678 * px))  # width, height for each plot

#############################################################################################
# Plot the first image on the left
axs[0].imshow(data, interpolation='none', extent=[0, N, 0, N], origin='lower', aspect='auto')
axs[0].axis('off')
axs[0].set_title('Filled Contour Plot')
axs[0].set_xlabel('x, [0,1]')
axs[0].set_ylabel('y, [0,1]')
#############################################################################################
# Plot the second image on the right
fourier = abs(np.fft.fftshift(np.fft.fft2(data))) ** 2
#print(data)
axs[1].imshow(np.log(fourier), cmap='gray')
#axs[1].imshow(F, cmap='gray')
axs[1].set_title('Plot of $\log{ ( | fftshift( fft2(z) ) |^2 } )$')
axs[1].set_xlabel('x, [0,1]')
axs[1].set_ylabel('y, [0,1]')
axs[1].axis('off')
#############################################################################################
fig.tight_layout()
plt.show()
