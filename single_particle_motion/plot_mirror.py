# Import libraries
import matplotlib.pyplot as plt
import numpy as np

# --- 1. Global Dark Mode & Presentation Settings ---
plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 14,
    'figure.facecolor': 'black',
    'axes.facecolor': 'black',
    'text.color': 'white',
    'axes.labelcolor': 'white',
    'xtick.color': 'white',
    'ytick.color': 'white',
    'figure.dpi': 150
})

# --- 2. Load Data ---
trajectory1 = np.loadtxt("mirrir_30.txt")

x1 = trajectory1[:300000, 1]
y1 = trajectory1[:300000, 2]
z1 = trajectory1[:300000, 3]

# --- 3. Setup Figure and Axes properly ---
fig = plt.figure(figsize=(10, 8))
# Assign the subplot to 'ax' (111 means a 1x1 grid, 1st subplot)
ax = fig.add_subplot(111, projection='3d')

# --- 4. Plot Data ---
# Using a bright neon cyan for high contrast on black
ax.plot(x1, y1, z1, label='Mirror Trajectory', color='#00E5FF', linewidth=1.5)

# Optional: Add a "shadow" projection below the plot for better 3D depth perception
#ax.plot(x1, y1, zs=np.min(z1), zdir='z', color='gray', alpha=0.3, label='z-projection')

# --- 5. Clean up 3D Background ---
# Make the background panes transparent
ax.xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
ax.yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
ax.zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))

# Make the grid lines a subtle dark grey
ax.xaxis._axinfo["grid"].update({"color": "#444444", "linestyle": "--"})
ax.yaxis._axinfo["grid"].update({"color": "#444444", "linestyle": "--"})
ax.zaxis._axinfo["grid"].update({"color": "#444444", "linestyle": "--"})

# --- 6. Set Labels and Legend ---
ax.set_title("Mirror Trajectory", fontsize=18, pad=20)
ax.set_xlabel("x", labelpad=10)
ax.set_ylabel("y", labelpad=10)
ax.set_zlabel("z", labelpad=10)
#ax.set_ylim(-3,3)

# Set legend with a black background to match
#ax.legend(loc="upper left", frameon=True, facecolor='black', edgecolor='none')

# Final Polish to prevent cut-off labels
fig.tight_layout()

# Save image (keeps the black background)
plt.savefig("mirror_trajectory.png", dpi=300, facecolor=fig.get_facecolor())

# Show plot
plt.show()