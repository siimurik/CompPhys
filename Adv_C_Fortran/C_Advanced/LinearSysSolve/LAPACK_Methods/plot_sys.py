import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt("out_2.csv")

x = a[:,0]
u = a[:,1]

#x = np.linspace(0, 2, len(u))
print(len(u))
print(u[0:5])
print(x[0:5])
#u = a[:,1]


plt.plot(x, u)
plt.grid()
plt.show()
