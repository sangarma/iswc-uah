import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ============================================================
# Electron Current Sheet Drift (Meandering Orbit) - 3D Plot
# Transparent Background & White Elements for Dark Presentation
# ============================================================

plt.style.use("default")

# ----------------------------
# Physical Parameters
# ----------------------------
B0 = 1.0            # Asymptotic magnetic field strength
L = 1.0             # Current sheet thickness parameter
m_e = 1.0           # Scaled electron mass
q_e = -1.0          # Electron charge

# Initial states [x, y, z, vx, vy, vz]
# Electron 1: Base velocity
state1 = [0.0, 0.0, 0.0, 0.0, 2.0, 0.0]

# Electron 2: Higher momentum/velocity
state2 = [0.0, 0.0, 0.0, 0.0, 4.0, 0.0]

# Time span for the numerical integration
t_span = (0, 40)
t_eval = np.linspace(t_span[0], t_span[1], 3000)

# ============================================================
# Numerical Solver for Non-Uniform Magnetic Field
# ============================================================
def lorentz_force(t, state, q, m):
    x, y, z, vx, vy, vz = state
    
    # Harris Current Sheet Magnetic Field: Reverses across y=0
    Bz = B0 * np.tanh(y / L)
    
    # F = q(v x B)
    ax = (q / m) * (vy * Bz)
    ay = (q / m) * (-vx * Bz)
    az = 0.0
    
    return [vx, vy, vz, ax, ay, az]

# Calculate Electron 1 Path (Base Velocity)
sol_e1 = solve_ivp(lorentz_force, t_span, state1, args=(q_e, m_e), t_eval=t_eval, method='RK45')
x_e1, y_e1, z_e1 = sol_e1.y[0], sol_e1.y[1], sol_e1.y[2]
vx_e1, vy_e1, vz_e1 = sol_e1.y[3], sol_e1.y[4], sol_e1.y[5]

# Calculate Electron 2 Path (Higher Velocity)
sol_e2 = solve_ivp(lorentz_force, t_span, state2, args=(q_e, m_e), t_eval=t_eval, method='RK45')
x_e2, y_e2, z_e2 = sol_e2.y[0], sol_e2.y[1], sol_e2.y[2]
vx_e2, vy_e2, vz_e2 = sol_e2.y[3], sol_e2.y[4], sol_e2.y[5]

# ============================================================
# Figure Setup
# ============================================================
fig = plt.figure(figsize=(14, 8), dpi=200)
fig.patch.set_alpha(0)

ax = fig.add_subplot(111, projection='3d')
ax.set_facecolor((0, 0, 0, 0))

# Transparent panes (removes the default grey box walls)
transparent = (1, 1, 1, 0)
ax.xaxis.set_pane_color(transparent)
ax.yaxis.set_pane_color(transparent)
ax.zaxis.set_pane_color(transparent)

# White grid lines (semi-transparent)
ax.grid(color='white', alpha=0.2)

# ============================================================
# Plot the Neutral Sheet (Current Sheet Plane)
# ============================================================
# Draw a faint plane at Y=0 to clearly show the magnetic reversal boundary
min_x = min(x_e1.min(), x_e2.min())
max_x = max(x_e1.max(), x_e2.max())
max_z = max(z_e1.max(), z_e2.max(), 1.0) # Ensure z is at least 1 for the plane

xx, zz = np.meshgrid(np.linspace(min_x, max_x, 10),
                     np.linspace(0, max_z, 10))
yy = np.zeros_like(xx)
ax.plot_surface(xx, yy, zz, color='white', alpha=0.1, edgecolor='none')

# ============================================================
# Plot Trajectories
# ============================================================
# Base velocity electron
ax.plot(x_e1, y_e1, z_e1, color='deepskyblue', linewidth=2.5, label=r'Vy = 2')

# Higher velocity electron
ax.plot(x_e2, y_e2, z_e2, color='springgreen', linewidth=2.5, label=r'Vy = 4')

# ============================================================
# Velocity Arrows at the end of the paths
# ============================================================
arrow_length = 3

# Arrow for Electron 1

# ============================================================
# White Labels, Axes, & Formatting
# ============================================================
ax.set_xlabel('X', fontsize=12, color='white', labelpad=10)
ax.set_ylabel('Y', fontsize=12, color='white', labelpad=10)
ax.set_zlabel('Z', fontsize=12, color='white', labelpad=10)

# Make tick marks and tick labels white
ax.tick_params(colors='white', labelsize=10)

# Make the physical axis lines white
ax.xaxis.line.set_color("white")
ax.yaxis.line.set_color("white")
ax.zaxis.line.set_color("white")
ax.set_xlim(-100, 0)
# Transparent legend with white text
leg = ax.legend(loc='upper left', fontsize=12, facecolor='none', edgecolor='white')
for text in leg.get_texts():
    text.set_color('white')

# ============================================================
# Viewing Angle and Limits
# ============================================================
# Angled to see the 3D shape: snaking across the Y=0 plane while drifting in X and Z
ax.view_init(elev=25, azim=-55)

# Axis limits as requested
ax.set_zlim(0, 10)
ax.set_xlim(-100, 0)

plt.tight_layout()

# ============================================================
# Save transparent PNG
# ============================================================
plt.savefig(
    "electron_current_sheet_drift_2_particles_3D.png",
    dpi=300,
    transparent=True,
    bbox_inches="tight"
)

plt.show()