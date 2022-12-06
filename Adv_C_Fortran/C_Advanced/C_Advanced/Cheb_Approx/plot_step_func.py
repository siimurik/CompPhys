import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data.dat")
fiq=plt.figure(figsize=(10,10))

p1=plt.subplot()
#plt.xlim(2.4,4.0)
#plt.ylim(0.0,1.0)
p1.grid()
#p1.set_aspect("equal")

x = data[:,0]
f = data[:,1]
f10 = data[:,2]
f40 = data[:,3]

plt.plot(x, f, color='k', label = 'Orig func', linewidth = 2)
plt.plot(x, f10, color='k', linestyle='--', label = '10th Order', linewidth = 1)
plt.plot(x, f40, color='k', linestyle='-.', label = '40th Order', linewidth = 0.5)
plt.xlabel('x')
plt.ylabel('f')
plt.title('Chebyshev approximations to a step function')
plt.legend()
plt.show()

