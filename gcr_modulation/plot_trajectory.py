# Import libraries
import matplotlib.pyplot as plt
import numpy as np

# Load data from file
trajectory1 = np.loadtxt("trajectory_init.dat")
trajectory2 = np.loadtxt("trajectory.dat")


x1 = trajectory1[:,0]
y1 = trajectory1[:,1]
z1 = trajectory1[:,2]
x2 = trajectory2[:,0]
y2 = trajectory2[:,1]
z2 = trajectory2[:,2]


# Plot data
ax = plt.figure().add_subplot(projection='3d')

ax.plot(x1, y1, z1, label='eta = 0.1')
ax.plot(x2, y2, z2, label='eta = 1')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
#ax.set_aspect('equal')
ax.legend()
#plt.savefig("both.png")
# Show plot
plt.show()