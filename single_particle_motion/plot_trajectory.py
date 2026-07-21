import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# 3D Trajectory Plot (Euler vs Heun vs RK4)
# Transparent Background & White Elements for Dark Presentation
# ============================================================

# Ensure default style so previous dark themes don't conflict
plt.style.use("default")

# Load data from file
trajectory1 = np.loadtxt("euler.txt")
trajectory2 = np.loadtxt("heun.txt")
trajectory3 = np.loadtxt("rk4.txt")

x1, y1, z1 = trajectory1[:,1], trajectory1[:,2], trajectory1[:,3]
x2, y2, z2 = trajectory2[:,1], trajectory2[:,2], trajectory2[:,3]
x3, y3, z3 = trajectory3[:,1], trajectory3[:,2], trajectory3[:,3]

# ============================================================
# Figure Setup
# ============================================================
fig = plt.figure(figsize=(10, 8), dpi=200)

# Make the outer figure background transparent
fig.patch.set_alpha(0)

# Create 3D axes and make the inner background transparent
ax = fig.add_subplot(111, projection='3d')
ax.set_facecolor((0, 0, 0, 0))

# Make the 3D panes (the grey walls) transparent
transparent = (1, 1, 1, 0)
ax.xaxis.set_pane_color(transparent)
ax.yaxis.set_pane_color(transparent)
ax.zaxis.set_pane_color(transparent)

# Add a faint white grid for depth reference
ax.grid(color='white', alpha=0.2)

# ============================================================
# Plot Data
# ============================================================
# Using high-contrast neon colors for dark backgrounds
ax.plot(x1, y1, z1, color='deepskyblue', linewidth=2.0, label='Euler')
ax.plot(x2, y2, z2, color='magenta', linewidth=2.0, label='Heun')
ax.plot(x3, y3, z3, color='springgreen', linewidth=2.0, label='RK4')

# ============================================================
# White Labels, Axes, & Formatting
# ============================================================
ax.set_xlabel("x", color='white', fontsize=12, labelpad=10)
ax.set_ylabel("y", color='white', fontsize=12, labelpad=10)
ax.set_zlabel("z", color='white', fontsize=12, labelpad=10)

# Make tick marks and tick labels white
ax.tick_params(colors='white', labelsize=10)

# Make the physical axis lines (the spine borders) white
ax.xaxis.line.set_color("white")
ax.yaxis.line.set_color("white")
ax.zaxis.line.set_color("white")

# Transparent legend with white text and a white border
leg = ax.legend(loc='best', fontsize=12, facecolor='none', edgecolor='white')
for text in leg.get_texts():
    text.set_color('white')

plt.tight_layout()

# ============================================================
# Save transparent PNG
# ============================================================
plt.savefig(
    "numerical_integration_comparison.png", 
    dpi=300, 
    transparent=True, 
    bbox_inches="tight"
)

# Show plot
plt.show()