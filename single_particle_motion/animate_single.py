# Import libraries
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

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

n_points = min(50000, len(trajectory1))

# Position data
x1 = trajectory1[:n_points, 1]
y1 = trajectory1[:n_points, 2]
z1 = trajectory1[:n_points, 3]

# --- 3. Setup Figure and Axes ---
# Adjusted figure size for a single plot
fig = plt.figure(figsize=(10, 8)) 

# Single 3D Plot
ax = fig.add_subplot(111, projection='3d')

# --- 4. Initialize Empty Line for Animation ---
line, = ax.plot([], [], [], label='Mirror Trajectory', color='#00E5FF', linewidth=1.5)

# --- 5. Clean up Backgrounds & Style ---
# 3D Plot Background
ax.xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
ax.yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
ax.zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
ax.xaxis._axinfo["grid"].update({"color": "#444444", "linestyle": "--"})
ax.yaxis._axinfo["grid"].update({"color": "#444444", "linestyle": "--"})
ax.zaxis._axinfo["grid"].update({"color": "#444444", "linestyle": "--"})

# --- 6. Set Limits, Labels and Titles ---
ax.set_title("Mirror Trajectory", fontsize=18, pad=10)
ax.set_xlabel("x", labelpad=10)
ax.set_ylabel("y", labelpad=10)
ax.set_zlabel("z", labelpad=10)

# CRITICAL FOR ANIMATION: Lock the axes limits
ax.set_xlim(x1.min(), x1.max())
ax.set_ylim(-5, 5)
ax.set_zlim(z1.min(), z1.max())

# --- 7. Define Animation Update Function ---
step_size = 50  # Number of points to draw per frame

def update(frame):
    idx = min(frame * step_size, n_points)
    
    # Update 3D trajectory
    line.set_data(x1[:idx], y1[:idx])
    line.set_3d_properties(z1[:idx])
    
    return line,

# --- 8. Create and Run Animation ---
total_frames = (n_points // step_size) + 1

ani = FuncAnimation(
    fig, 
    update, 
    frames=total_frames, 
    interval=30,      
    blit=False        
)

fig.tight_layout()

# Optional: Save as an mp4 for your presentation
#ani.save("mirror_trajectory_animation.mp4", writer='ffmpeg', fps=30, savefig_kwargs={'facecolor': fig.get_facecolor()})
ani.save("mirror_trajectory_animation.gif", writer='pillow', fps=30, savefig_kwargs={'facecolor': 'black'})
plt.show()