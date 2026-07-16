# Import libraries
import matplotlib.pyplot as plt
import numpy as np

# Load data from file
trajectory1 = np.loadtxt("sdrift.txt")
#trajectory2 = np.loadtxt("cdrift_mom.txt")
#trajectory3 = np.loadtxt("rk4.txt")

x1 = trajectory1[:,1]
y1 = trajectory1[:,2]
z1 = trajectory1[:,3]
#x2 = trajectory2[:,1]
#y2 = trajectory2[:,2]
#z2 = trajectory2[:,3]
#x3 = trajectory3[:,1]
#y3 = trajectory3[:,2]
#z3 = trajectory3[:,3]
energy = (np.array(trajectory1[:,-1]))**2 + (np.array(trajectory1[:,-2]))**2 + (np.array(trajectory1[:,-3]))**2

# Plot data
ax = plt.figure().add_subplot(projection='3d')
ax.plot(x1, y1, z1, label='shock_drift')
#ax.plot(x2, y2, z2, label='p = 4')
#ax.plot(x3, y3, z3, label='rk4')


# Set axes labels and legend
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
# ax.set_aspect('equal')
ax.legend()
#plt.savefig("both.png")
# Show plot
plt.show()