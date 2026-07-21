import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# ExB Drift of Charged Particles
# Transparent Background & White Elements for Dark Presentation
# ============================================================

plt.style.use("default")

# ----------------------------
# Physical Parameters
# ----------------------------
B_field = 1.0       # B-field in +Z direction
E_field = 0.5       # E-field in +Y direction

m_e = 1.0
m_p = 5.0

q_e = -1.0
q_p = +1.0

# Initial velocity [vx, vy, vz]
# Starting from rest in XY plane creates beautiful cycloids
v0 = np.array([0.0, 0.0, 0.15])

# Time array
t = np.linspace(0, 100, 3000)

# ============================================================
# Exact analytical solution for crossed E and B fields
# ============================================================
def calc_exb_kinematics(m, q, v0, B, E, t):
    
    omega = q * B / m
    vx0, vy0, vz0 = v0
    
    # Drift velocity magnitude (E x B drift is in +X direction)
    v_drift = E / B
    
    # Transform initial velocity to the guiding center drift frame
    v_prime_x0 = vx0 - v_drift
    v_prime_y0 = vy0
    
    wt = omega * t
    
    # Velocity (transform back to lab frame)
    vx = v_prime_x0 * np.cos(wt) + v_prime_y0 * np.sin(wt) + v_drift
    vy = -v_prime_x0 * np.sin(wt) + v_prime_y0 * np.cos(wt)
    vz = np.full_like(t, vz0)
    
    # Position
    x = (v_prime_x0 / omega) * np.sin(wt) - (v_prime_y0 / omega) * (np.cos(wt) - 1) + (v_drift * t)
    y = (v_prime_x0 / omega) * (np.cos(wt) - 1) + (v_prime_y0 / omega) * np.sin(wt)
    z = vz0 * t
    
    return x, y, z, vx, vy, vz

# Electron Path
x_e, y_e, z_e, vx_e, vy_e, vz_e = calc_exb_kinematics(
    m_e, q_e, v0, B_field, E_field, t
)

# Proton Path
x_p, y_p, z_p, vx_p, vy_p, vz_p = calc_exb_kinematics(
    m_p, q_p, v0, B_field, E_field, t
)

# ============================================================
# Figure Setup
# ============================================================
fig = plt.figure(figsize=(12, 8), dpi=200)

# Transparent figure background
fig.patch.set_alpha(0)

ax = fig.add_subplot(111, projection='3d')
ax.set_facecolor((0, 0, 0, 0))

# Transparent panes (removes the grey walls)
transparent = (1, 1, 1, 0)
ax.xaxis.set_pane_color(transparent)
ax.yaxis.set_pane_color(transparent)
ax.zaxis.set_pane_color(transparent)

# White grid lines (semi-transparent)
ax.grid(color='white', alpha=0.2)

# ============================================================
# Plot trajectories
# ============================================================
ax.plot(
    x_e, y_e, z_e,
    color='deepskyblue',
    linewidth=2.5,
    label='Electron ($e^{-}$)'
)

ax.plot(
    x_p, y_p, z_p,
    color='magenta',
    linewidth=2.5,
    label='Proton ($p^{+}$)'
)

# ============================================================
# Velocity arrows at the end of the path
# ============================================================
arrow_length = 3

ax.quiver(
    x_e[-1], y_e[-1], z_e[-1],
    vx_e[-1], vy_e[-1], vz_e[-1],
    color='deepskyblue', normalize=True,
    length=arrow_length, arrow_length_ratio=0.4, linewidth=2
)

ax.quiver(
    x_p[-1], y_p[-1], z_p[-1],
    vx_p[-1], vy_p[-1], vz_p[-1],
    color='magenta', normalize=True,
    length=arrow_length, arrow_length_ratio=0.4, linewidth=2
)

# ============================================================
# Field Vectors (E and B)
# ============================================================
mid_x = 25 
zmax = max(z_e.max(), z_p.max())

# Magnetic Field (Gold, +Z direction)
ax.quiver(
    mid_x, 0, 0, 
    0, 0, zmax + 5,
    color='gold', linewidth=3.0,
    arrow_length_ratio=0.08, label=r'Magnetic Field $\mathbf{B}$'
)

# Electric Field (LimeGreen, +Y direction)
ax.quiver(
    mid_x, 0, 0, 
    0, 5, 0,  
    color='limegreen', linewidth=3.0,
    arrow_length_ratio=0.2, label=r'Electric Field $\mathbf{E}$'
)

# ============================================================
# White Labels & Axis Lines
# ============================================================
ax.set_title(
    r'$\mathbf{E} \times \mathbf{B}$ Drift of Charged Particles',
    fontsize=18,
    pad=20,
    fontweight='bold',
    color='white'
)

ax.set_xlabel('X', fontsize=12, color='white')
ax.set_ylabel('Y', fontsize=12, color='white')
ax.set_zlabel('Z', fontsize=12, color='white')

# Make tick marks and tick labels white
ax.tick_params(colors='white', labelsize=10)

# Make the physical axis lines white
ax.xaxis.line.set_color("white")
ax.yaxis.line.set_color("white")
ax.zaxis.line.set_color("white")

# Adjust axes limits
ax.set_ylim(-3, 6)

# Transparent legend with white text
leg = ax.legend(loc='upper left', fontsize=11, facecolor='none', edgecolor='white')
for text in leg.get_texts():
    text.set_color('white')

# ============================================================
# Viewing angle
# ============================================================
ax.view_init(elev=20, azim=180-45)

plt.tight_layout()

# ============================================================
# Save transparent PNG
# ============================================================
plt.savefig(
    "exb_drift_transparent_white_axes.png",
    dpi=300,
    transparent=True,
    bbox_inches="tight"
)

plt.show()