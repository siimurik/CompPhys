import numpy as np
import matplotlib.pyplot as plt

# Replace 'your_file.csv' with the actual path to your CSV file
#file_path = 'colsRevamp.csv'
file_path = 'c1_values.csv'

# Use numpy.genfromtxt to load CSV values into a NumPy matrix
dataRaw = np.genfromtxt(file_path, delimiter=',', filling_values=np.nan)    # shape 512x513
data = dataRaw[:, :len(dataRaw)]  # Removing last column bc it is full of 'nan' values; 
print(data.shape)
print(data)
N = data.shape[0]

px = 1/72# pixel in inches
plt.subplots(figsize=(660*px, 678*px)) #width, height. These numbers give exactly 512x512 pixel image, therefore containing no aliasing(?) errors, painstakingly adjusted through much trial and error

#fig = plt.figure(figsize=(5,5), dpi = 128) #lim dpi = ca 156 kanti
#for now it is in grayscale, can change
colbar = plt.imshow(data, interpolation='none', extent=[0, N, 0, N], origin='lower', aspect='auto')
#plt.clim(-1, 1) #fixes limits of color bar
#cmap variants: Greys, hot, terrain, rainbow, cividis, viridis, gray. last one most similar to lecture example
#fig.colorbar(colbar)
plt.axis('off')

plt.savefig('my_fig.png', dpi=72, bbox_inches='tight') #should save the image in working folder, with small borders that need to be cropped out manually

plt.show()


