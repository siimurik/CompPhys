import numpy as np

A = np.zeros([1000000])
B = np.zeros([1000000])
for i in range(1000):
    for j in range(1000):
            A[i * 1000 + j] = i*j
            B[i * 1000 + j] = i + j

C = np.matmul(A,B)

print(  A[0:5],
        A[1000:1005],
        A[2000:2005],
        A[3000:4005])
#for i in range(5):
#    print(A[int(1000*i):int(1005*i)])