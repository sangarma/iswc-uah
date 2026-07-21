import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

# --- 1. Global Presentation Settings ---
plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 18,
    'axes.labelsize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 14,
    'figure.dpi': 150
})

# --- 2. Load Data ---
trajectory1 = np.loadtxt("sdrift.txt")

# Limit to 10000 points
n_points = min(10000, len(trajectory1))
x1 = trajectory1[:n_points, 1]
y1 = trajectory1[:n_points, 2]
z1 = trajectory1[:n_points, 3]
p1 = trajectory1[:n_points, 4]
p2 = trajectory1[:n_points, 5]
p3 = trajectory1[:n_points, 6]

# Calculate energy
energy = 0.5 * (p1**2 + p2**2 + p3**2)
time_steps = np.arange(n_points)

# Surface data
y_range = np.linspace(-5, 35, 1000)
z_range = np.linspace(-0.04, 0.04, 1000)
Y, Z = np.meshgrid(y_range, z_range)
X = np.zeros_like(Y)

# --- 3. Setup Figure and Axes ---
fig = plt.figure(figsize=(14, 6))

ax1 = fig.add_subplot(1, 2, 1, projection='3d')
ax2 = fig.add_subplot(1, 2, 2)

ax1.set_title("Particle Trajectory")
ax2.set_title("Energy Conservation")

# Draw the static surface plane once
ax1.plot_surface(X, Y, Z, color='#457B9D', alpha=0.3, rstride=100, cstride=100, linewidth=0)

# Clean up 3D background panes
ax1.xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax1.yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax1.zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

# Make 3D grid lines lighter and dashed
ax1.xaxis._axinfo["grid"].update({"color": "lightgray", "linestyle": "--"})
ax1.yaxis._axinfo["grid"].update({"color": "lightgray", "linestyle": "--"})
ax1.zaxis._axinfo["grid"].update({"color": "lightgray", "linestyle": "--"})

# Clean up the 2D plot (spines and grid)
ax2.grid(True, linestyle='--', color='gray', alpha=0.3)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Set labels with padding
ax1.set_xlabel("x", labelpad=10)
ax1.set_ylabel("y", labelpad=10)
ax1.set_zlabel("z", labelpad=10)
ax2.set_xlabel("Time step")
ax2.set_ylabel("Energy")

# --- 4. Set Limits (Required for animation) ---
ax1.set_xlim(x1.min(), x1.max())
ax1.set_ylim(y1.min(), y1.max())
ax1.set_zlim(z1.min(), z1.max())

ax2.set_xlim(0, n_points)
e_margin = (energy.max() - energy.min()) * 0.1
if e_margin == 0: e_margin = 0.1
ax2.set_ylim(energy.min() - e_margin, energy.max() + e_margin)

# --- 5. Initialize Empty Line Objects ---
# Using presentation styling (colors and linewidths) for the animated lines
line1, = ax1.plot([], [], [], label='Trajectory', color='#E63946', linewidth=2)
line2, = ax2.plot([], [], label='Total Energy', color='#1D3557', linewidth=2.5)

ax1.legend(loc="upper right", frameon=False)
ax2.legend(loc="upper right", frameon=False)

# --- 6. Define the Update Function ---
step_size = 2  # Plot 20 new points per frame

def update(frame):
    # Calculate the current index based on the frame and step_size
    idx = min(frame * step_size, n_points)
    
    # Update 3D trajectory
    line1.set_data(x1[:idx], y1[:idx])
    line1.set_3d_properties(z1[:idx])
    
    # Update Energy plot
    line2.set_data(time_steps[:idx], energy[:idx])
    
    return line1, line2

# --- 7. Create the Animation ---
total_frames = (n_points // step_size) + 1

ani = FuncAnimation(
    fig, 
    update, 
    frames=total_frames, 
    interval=100,
    blit=False
)

# Optional: Save for a presentation slide (e.g., PowerPoint allows inserting MP4s)
# ani.save("presentation_animation.mp4", writer='ffmpeg', fps=30, dpi=150)

fig.tight_layout(pad=3.0)
plt.show()