# Import libraries
import matplotlib.pyplot as plt
import numpy as np

# Load data from file
trajectory1 = np.loadtxt("mirrir_30.txt")


p1 = np.array(trajectory1[:,4])
p2 = np.array(trajectory1[:,5])
p3 = np.array(trajectory1[:,6])

energy1 = 0.5 * (p3**2)  # Kinetic energy formula
energy2 = 0.5*(p1**2 + p2**2+p3**2)  # Kinetic energy formula
plt.figure(figsize=(10, 6))
plt.plot(energy1, label='E parallel', color='blue')
plt.plot(energy2, label='E_perpendicular', color='orange')
plt.title('Kinetic Energy vs Time')
plt.xlabel('Time Steps')
plt.ylabel('Kinetic Energy')
plt.grid()
plt.legend()
#plt.savefig('kinetic_energy_plot.png')  # Save the plot as a PNG file
plt.show()  # Display the plot