import numpy as np
import matplotlib.pylab as plt
plt.style.use('dark_background')

# Define wavelength in meters
#lam = np.arange(1000.0*1e-10, 20000.*1e-10, 20e-10)
lam = np.arange(100.*1e-09, 2500*10**-9, 1e-09)

h = 6.626068e-34 # m^2 kg / s
c = 2.99792458e8 # m / s
k = 1.380658e-23 # J / K Boltzmann constant

def planck(landa, T): 
    return (2*np.pi*h*c**2)/landa**5 * 1/(np.exp(h*c/(landa*k*T)) - 1.)

plt.xlabel("Lainepikkus [$m$]",fontsize=20)
plt.ylabel("Kiirgustihedus [$\\frac{W}{m^2 \cdot sr \cdot \\frac{1}{s}}$]",fontsize=20)
plt.title("Plancki jaotus",fontsize=20)
plt.plot(lam, planck(lam, 7000.), label='7000 K')
plt.plot(lam, planck(lam, 6000.), label='6000 K')
plt.plot(lam, planck(lam, 5000.), label='5000 K')
plt.plot(lam, planck(lam, 4000.), label='4000 K')
plt.legend(fontsize=20)
plt.grid()
plt.show()
