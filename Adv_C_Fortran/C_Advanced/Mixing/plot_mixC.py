import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("colors.csv",header = None)
#print('\n Initial data\n',data.head())
origin = 'lower'
dim = len(data)
T = data.drop(data.columns[[dim]], axis=1)
#print('\n Initial data\n',T.head())
#print(x_feat.head())
#print(x_feat.tail())
Z = T.values
#print(Z)
N = len(T)
#print(nx)
"""
px = 1/72# pixel in inches
plt.subplots(figsize=(660*px, 678*px)) #width, height. These numbers give exactly 512x512 pixel image, therefore containing no aliasing(?) errors, painstakingly adjusted through much trial and error

#fig = plt.figure(figsize=(5,5), dpi = 128) #lim dpi = ca 156 kanti
#for now it is in grayscale, can change
colbar = plt.imshow(Z, cmap='gray', interpolation='none', extent=[0, N, 0, N], origin='lower', aspect='auto')
#plt.clim(-1, 1) #fixes limits of color bar
#cmap variants: Greys, hot, terrain, rainbow, cividis, viridis, gray. last one most similar to lecture example
#fig.colorbar(colbar)
plt.axis('off')

plt.savefig('my_fig.png', dpi=72, bbox_inches='tight') #should save the image in working folder, with small borders that need to be cropped out manually

plt.show()
"""

# Creating 2-D grid of features
feature_x = np.linspace(0, 1, N)
feature_y = np.linspace(0, 1, N)
#print(feature_x)

# Creating 2-D grid of features
[X, Y] = np.meshgrid(feature_x, feature_y)

fig, ax = plt.subplots(1, 1)
# plots filled contour plot
ax.contourf(X, Y, Z)
  
ax.set_title('Filled Contour Plot')
ax.set_xlabel('feature_x')
ax.set_ylabel('feature_y')
  
plt.show()
