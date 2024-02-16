import matplotlib.pyplot as plt

# Read the data
with open('rk2.dat', 'r') as file:
    data = file.readlines()

# Extracting columns
column1 = [float(row.split()[0]) for row in data]
column2 = [float(row.split()[1]) for row in data]
column3 = [float(row.split()[2]) for row in data]
column4 = [float(row.split()[3]) for row in data]
column5 = [float(row.split()[4]) for row in data]

# Plotting
plt.figure(figsize=(6, 4))

plt.subplot(221)
plt.plot(column1, column2)
plt.title('x(t)')

plt.subplot(222)
plt.plot(column1, column3)
plt.title('y(t)')

plt.subplot(223)
plt.plot(column1, column4)
plt.title('$v_{x}(t)$')

plt.subplot(224)
plt.plot(column1, column5)
plt.title('$v_{y}(t)$')

plt.tight_layout()
#plt.grid()
plt.show()
