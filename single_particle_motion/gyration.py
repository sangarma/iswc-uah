import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Charged Particle Gyration in a Uniform Magnetic Field
# B = B * z_hat
# Exact analytical solution of:
#     m dv/dt = q (v x B)
# ============================================================

plt.style.use('dark_background')

# ----------------------------
# Physical Parameters
# ----------------------------
B_field = 1.0

m_e = 1.0
m_p = 5.0

q_e = -1.0
q_p = +1.0

# Initial velocity [vx, vy, vz]
v0 = np.array([2.0, 0.0, 0.5])

# Time array
t = np.linspace(0, 150, 3000)


# ============================================================
# Analytical solution
# ============================================================
def calc_kinematics(m, q, v0, B, t):

    omega = q * B / m

    vx0, vy0, vz0 = v0

    wt = omega * t

    # -------------------------
    # Velocity
    # -------------------------
    vx = vx0 * np.cos(wt) + vy0 * np.sin(wt)

    vy = -vx0 * np.sin(wt) + vy0 * np.cos(wt)

    vz = np.full_like(t, vz0)

    # -------------------------
    # Position
    # (Assuming x=y=z=0 at t=0)
    # -------------------------
    x = (vx0 / omega) * np.sin(wt) \
        + (vy0 / omega) * (np.cos(wt) - 1)

    y = (vx0 / omega) * (np.cos(wt) - 1) \
        + (vy0 / omega) * np.sin(wt)

    z = vz0 * t

    return x, y, z, vx, vy, vz


# Electron
x_e, y_e, z_e, vx_e, vy_e, vz_e = calc_kinematics(
    m_e, q_e, v0, B_field, t
)

# Proton
x_p, y_p, z_p, vx_p, vy_p, vz_p = calc_kinematics(
    m_p, q_p, v0, B_field, t
)

# ============================================================
# Plot
# ============================================================
fig = plt.figure(figsize=(10, 8), dpi=150)

ax = fig.add_subplot(111, projection='3d')

# Black panes
ax.xaxis.set_pane_color((0, 0, 0, 1))
ax.yaxis.set_pane_color((0, 0, 0, 1))
ax.zaxis.set_pane_color((0, 0, 0, 1))

ax.grid(color='#444444', linestyle='--', linewidth=0.5)

# ============================================================
# Trajectories
# ============================================================
ax.plot(
    x_e, y_e, z_e,
    color='cyan',
    linewidth=2,
    label='Electron ($e^{-}$)'
)

ax.plot(
    x_p, y_p, z_p,
    color='magenta',
    linewidth=2,
    label='Proton ($p^{+}$)'
)

# ============================================================
# Arrow showing instantaneous velocity
# ============================================================
arrow_length = 5

ax.quiver(
    x_e[-1], y_e[-1], z_e[-1],
    vx_e[-1], vy_e[-1], vz_e[-1],
    color='cyan',
    normalize=True,
    length=arrow_length,
    arrow_length_ratio=0.4,
    linewidth=2
)

ax.quiver(
    x_p[-1], y_p[-1], z_p[-1],
    vx_p[-1], vy_p[-1], vz_p[-1],
    color='magenta',
    normalize=True,
    length=arrow_length,
    arrow_length_ratio=0.4,
    linewidth=2
)

# ============================================================
# Magnetic field arrow
# ============================================================
zmax = max(z_e.max(), z_p.max())

ax.quiver(
    0, 0, 0,
    0, 0, zmax + 5,
    color='yellow',
    linewidth=2,
    arrow_length_ratio=0.08,
    label='Magnetic Field ($\\mathbf{B}$)'
)

# ============================================================
# Labels
# ============================================================
ax.set_title(
    "Electron and Proton Gyration in a Uniform Magnetic Field",
    fontsize=16,
    pad=20
)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z (Along Magnetic Field)")

ax.legend()

# Nice viewing angle
ax.view_init(elev=20, azim=180)
ax.grid(False)
plt.tight_layout()
plt.show()