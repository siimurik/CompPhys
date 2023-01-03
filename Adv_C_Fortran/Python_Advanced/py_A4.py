import time
import numpy as np

# Solve the problem
# -u''(x) + (4 - x)u(x) = x + 5 , x ∈ (0, 2) ,
# u'(0) = 1, u'(2) = 2,
# by means of the Galerkin finite element method
# where xi = ih, h = 0.01

n = 5000
h = 2/n
x = np.arange(0, 2+h, h)
A = np.zeros((n+1,n+1))
y = np.zeros(n+1)
A[0,0] = 1/h + h/2*(4-x[0])
A[0,1] = -1/h
y[0] = h/2*(x[0]+5) - 1

for i in range(1, n):
    A[i,i-1] = -1/h
    A[i,i]   =  2/h + h*(4-x[i])
    A[i,i+1] = -1/h
    y[i]     = h*(x[i]+5)

A[n,n-1] = -1/h
A[n,n]   = 1/h + h/2*(4-x[n])
y[n]     = h/2*(x[n]+5) - 1

start_time = time.time()

u = np.linalg.solve(A, y)

elapsed_time = time.time() - start_time
print('Elapsed time: ', np.round(elapsed_time,6),'seconds.')

# Print out the values of A, y, and u
print("A:")
print(A)
print(len(A), len(A[0]))
print("y:")
print(y[0:5],'\n',y[-5:])
print("u:")
print(u[0:5],'\n',u[-5:])

import matplotlib.pyplot as plt
plt.plot(x, u)
plt.xlabel('x')
plt.ylabel('u')
plt.grid()
plt.show()
