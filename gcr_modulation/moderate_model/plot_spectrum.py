# Import libraries
import matplotlib.pyplot as plt
import numpy as np

# Load data from file
spectrum1 = np.loadtxt("au1.dat")
spectrum2 = np.loadtxt("au10.dat")
spectrum3 = np.loadtxt("au40.dat")
spectrum4 = np.loadtxt("au70.dat")
spectrum5 = np.loadtxt("au80.dat")
# Prepare arrays
Energy1 = spectrum1[:,0]
ModulatedSpectrum1 = spectrum1[:,1]
UnmodulatedSpectrum1 = spectrum1[:,2]

Energy2 = spectrum2[:,0]
ModulatedSpectrum2 = spectrum2[:,1]
UnmodulatedSpectrum2 = spectrum2[:,2]


Energy3 = spectrum3[:,0]
ModulatedSpectrum3 = spectrum3[:,1]
Energy4 = spectrum4[:,0]
ModulatedSpectrum4 = spectrum4[:,1]
Energy5 = spectrum5[:,0]
ModulatedSpectrum5 = spectrum5[:,1]
# Plot data
ax = plt.figure().add_subplot()
ax.loglog(Energy1, ModulatedSpectrum1, label='au_1')
ax.loglog(Energy4, ModulatedSpectrum4, label='au_70')
ax.loglog(Energy2, ModulatedSpectrum2, label='au_10')
ax.loglog(Energy3, ModulatedSpectrum3, label='au_40')
ax.loglog(Energy5, ModulatedSpectrum5, label='au_80')

# Set axes labels and legend
ax.set_xlabel("Energy")
ax.set_ylabel("Spectrum")
# ax.set_aspect('equal')
ax.legend()

# Show plot
plt.show()