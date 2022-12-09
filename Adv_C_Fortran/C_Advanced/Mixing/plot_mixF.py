import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Reading in the data
data = pd.read_csv("colors.csv",header = None)
fdata = pd.read_csv("fourier.csv",header = None)
#print('\n Initial data\n',fdata.head())
origin = 'lower'
dim = len(data)
T = data.drop(data.columns[[dim]], axis=1)
#print('\n Initial data\n',T.head())
#print(x_feat.head())
Z = T.values
fourier = abs(np.fft.fftshift(np.fft.fft2(Z)))**2
#print(np.log(fourier))
F = fdata.values
#print(Z)
N = len(T)
#====================================================================
# Plotting
# Creating 2-D grid of features
feature_x = np.linspace(0, 1, N)
feature_y = np.linspace(0, 1, N)
#print(feature_x)

# Creating 2-D grid of features
[X, Y] = np.meshgrid(feature_x, feature_y)

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
axes[0].contourf(X, Y, Z)
axes[0].set_title('Filled Contour Plot')
axes[0].set_xlabel('x, [0,1]')
axes[0].set_ylabel('y, [0,1]')
axes[1].imshow(F, cmap='gray')
axes[1].set_title('Plot of $\log{ ( | fftshift( fft2(z) ) |^2 } )$')
axes[1].set_xlabel('x, [0,1]')
axes[1].set_ylabel('y, [0,1]')
fig.tight_layout()

plt.show()
#====================================================================

#print(nx)
"""
px = 1/72# pixel in inches
plt.subplots(figsize=(660*px, 678*px)) #width, height. These numbers give exactly 512x512 pixel image, therefore containing no aliasing(?) errors, painstakingly adjusted through much trial and error

#fig = plt.figure(figsize=(5,5), dpi = 128) #lim dpi = ca 156 kanti
#for now it is in grayscale, can change
colbar = plt.imshow(F, interpolation='none', extent=[0, N, 0, N], origin='lower', aspect='auto')
#plt.clim(-1, 1) #fixes limits of color bar
#cmap variants: Greys, hot, terrain, rainbow, cividis, viridis, gray. last one most similar to lecture example
#fig.colorbar(colbar)
plt.axis('off')

plt.savefig('my_fig.png', dpi=72, bbox_inches='tight') #should save the image in working folder, with small borders that need to be cropped out manually

plt.show()
"""

#fig, ax = plt.subplots(5, 5)
#plt.subplot(1, 2, 1)
# plots filled contour plot
#plt.contourf(X, Y, Z)
#colbar = plt.imshow(T, interpolation='none', extent=[0, N, 0, N], origin='lower', aspect='auto')

#plt.set_title('Filled Contour Plot')
#plt.set_xlabel('feature_x')
#plt.set_ylabel('feature_y')

#plt.subplot(2, 2, 2)
#plt.imshow(F, cmap='gray')

#plt.show()
#px = 1/72# pixel in inches
#plt.subplots(figsize=(660*px, 678*px))
#plt.subplots(figsize=(660*px, 678*px)) #width, height. These numbers give exactly 512x512 pixel image, therefore containing no aliasing errors, painstakingly adjusted through much trial and error
#colbar = plt.imshow(np.log(fourier), cmap='gray')
#colbar = plt.imshow(F, cmap='gray')

#fig.colorbar(colbar)
#print(fourier)

#plt.axis('off')
#plt.savefig('fourier.png', dpi=72, bbox_inches='tight')
