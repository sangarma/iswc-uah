# Import libraries
import matplotlib.pyplot as plt
import numpy as np

# --- 1. Global Dark Mode Presentation Settings ---
plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 18,
    'axes.labelsize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 14,
    'figure.dpi': 150,
    
    # Dark Mode Colors
    'figure.facecolor': 'black',    # Figure background
    'axes.facecolor': 'black',      # Axes background
    'text.color': 'white',          # Title and general text color
    'axes.labelcolor': 'white',     # X/Y/Z label colors
    'xtick.color': 'white',         # Tick marks color
    'ytick.color': 'white'          # Tick marks color
})

# --- 2. Load Data ---
trajectory1 = np.loadtxt("sdrift.txt")

x1 = trajectory1[:10000, 1]
y1 = trajectory1[:10000, 2]
z1 = trajectory1[:10000, 3]
p1 = trajectory1[:10000, 4]
p2 = trajectory1[:10000, 5]
p3 = trajectory1[:10000, 6]
energy = 0.5 * (p1**2 + p2**2 + p3**2)

# Surface data
y_range = np.linspace(-5, 35, 1000)
z_range = np.linspace(-0.04, 0.04, 1000)
Y, Z = np.meshgrid(y_range, z_range)
X = np.zeros_like(Y)

# --- 3. Setup Figure ---
fig = plt.figure(figsize=(14, 6))

# --- 4. 3D Trajectory Plot ---
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
ax1.set_title("Particle Trajectory")

# Bright red/pink for high contrast on black
ax1.plot(x1, y1, z1, label='Trajectory', color='#FF4D4D', linewidth=2) 

# Surface plotted in a bright cyan
ax1.plot_surface(X, Y, Z, color='#4CC9F0', alpha=0.3, rstride=100, cstride=100, linewidth=0)

# Make 3D background panes transparent so the black figure background shows through
ax1.xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
ax1.yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
ax1.zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))

# Make 3D grid lines a dark grey so they don't distract
ax1.xaxis._axinfo["grid"].update({"color": "#444444", "linestyle": "--"})
ax1.yaxis._axinfo["grid"].update({"color": "#444444", "linestyle": "--"})
ax1.zaxis._axinfo["grid"].update({"color": "#444444", "linestyle": "--"})

# Axes labels & Legend
ax1.set_xlabel("x", labelpad=10)
ax1.set_ylabel("y", labelpad=10)
ax1.set_zlabel("z", labelpad=10)

# Set legend frame edge color to none and facecolor to black
#ax1.legend(loc="upper right", frameon=True, facecolor='black', edgecolor='none')

# --- 5. 2D Energy Plot ---
ax2 = fig.add_subplot(1, 2, 2)
ax2.set_title("Energy")

# Bright neon cyan for the 2D plot line
ax2.plot(energy, label='Total Energy', color='#00E5FF', linewidth=2.5)

ax2.set_xlabel("Time step")
ax2.set_ylabel("Energy")

# 2D Grid and Spines tailored for dark mode
ax2.grid(True, linestyle='--', color='#444444', alpha=0.8) # Dark grey grid
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_color('white')
ax2.spines['left'].set_color('white')

#ax2.legend(loc="upper right", frameon=True, facecolor='black', edgecolor='none')

# --- 6. Final Polish ---
fig.tight_layout(pad=3.0)

# Save the figure. facecolor='black' ensures the saved PNG keeps the black background.
plt.savefig("shock_dark.png", dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor(), transparent=True)

plt.show()