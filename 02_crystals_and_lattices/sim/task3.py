hbar=1
a=0.5
m=1
P = [0.1,1,5,20]
q_vals = np.linspace(0.01, 25.0, 4000)
for idx,p in enumerate(P):
    reduced_k = []
    reduced_E = []
    extended_k = []
    extended_E = []
    for q in q_vals:
        val = np.cos(q*a)+ ( p/ (q*a) ) * np.sin(q*a)  # cos(ka)

        if -1<=  val <= 1:
            k = np.arccos(val)/a
            E = hbar * q**2 / (2*m)

            # Reduced zone scheme (-pi/a to pi/a)
            for sign in [-1, 1]:
                reduced_k.append(sign * k)
                reduced_E.append(E)

           # Extended zone scheme (-3pi/a to 3pi/a)
            for n in [-2, -1, 0, 1, 2]:
                extended_k.append(sign * k + n * (2 * np.pi / a))
                extended_E.append(E)


fig, axes = plt.subplots(1, 2, figsize=(14, 5))
# 1. Reduced Zone Scheme
axes[0].scatter(
    np.array(reduced_k) * a / np.pi,
    reduced_E,
    s=2,
    color='royalblue',
    alpha=0.7,
)
axes[0].set_title(f'Reduced Zone Scheme ($P = {P}$)', fontsize=12)
axes[0].set_xlabel(r'Crystal Momentum $k \cdot (a / \pi)$', fontsize=11)
axes[0].set_ylabel('Energy $E$ (eV)', fontsize=11)
axes[0].set_xlim(-1, 1)
axes[0].set_ylim(0, 100)
axes[0].axvline(0, color='gray', linestyle=':', alpha=0.7)
axes[0].grid(True, linestyle='--', alpha=0.5)

# 2. Extended Zone Scheme
axes[1].scatter(
    np.array(extended_k) * a / np.pi,
    extended_E,
    s=1,
    color='darkorange',
    alpha=0.5,
)
axes[1].set_title(f'Extended Zone Scheme ($P = {P}$)', fontsize=12)
axes[1].set_xlabel(r'Crystal Momentum $k \cdot (a / \pi)$', fontsize=11)
axes[1].set_ylabel('Energy $E$ (eV)', fontsize=11)
axes[1].set_xlim(-3, 3)
axes[1].set_ylim(0, 100)
axes[1].axvline(0, color='gray', linestyle=':', alpha=0.7)
axes[1].grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.show()




####### attempting it in a combined way ###
P_kp = 1.0
q_vals = np.linspace(0.01, 25.0, 4000)
kp_k, kp_E = [], []

for q in q_vals:
  rhs = np.cos(q * a) + (P_kp / (q * a)) * np.sin(q * a)
  if -1.0 <= rhs <= 1.0:
    k = np.arccos(rhs) / a
    E = (hbar**2) * (q**2) / (2 * m)
    for sign in [-1, 1]:
      kp_k.append(sign * k)
      kp_E.append(E)

N_cell = 100
dx = a / N_cell
t_fd = hbar**2 / (2 * m * dx**2)
V0_fd = 1.5  # Equivalent potential strength

x_cell = np.linspace(0, a, N_cell, endpoint=False)
V_x = 2 * V0_fd * np.cos(2 * np.pi * x_cell / a)

fd_energies = []
for k in k_grid:
  H = np.zeros((N_cell, N_cell), dtype=complex)
  for i in range(N_cell):
    H[i, i] = 2 * t_fd + V_x[i]
    if i < N_cell - 1:
      H[i, i + 1] = -t_fd
      H[i + 1, i] = -t_fd
    else:
      H[i, 0] = -t_fd * np.exp(1j * k * a)
      H[0, i] = -t_fd * np.exp(-1j * k * a)
  evals, _ = np.linalg.eigh(H)
  fd_energies.append(evals)
fd_energies = np.array(fd_energies)


# Plot Real-Space Finite-Difference (Lines)

plt.scatter(
    np.array(kp_k) * a / np.pi,
    kp_E,
    s=4,
    color='royalblue',
    alpha=0.6,
    label='Kronig-Penney (Analytical)',
)
for band in range(min(3, N_cell)):
  label_str = (
      'Real-Space FD (Numerical)' if band == 0 else ''
  )  # Label only once
  plt.plot(
      k_grid * a / np.pi,
      fd_energies[:, band],
      color='darkorange',
      linestyle='-',
      linewidth=2,
      label=label_str,
  )
    
plt.title(
    'Comparison: Kronig-Penney vs. Real-Space Finite-Difference', fontsize=12
)
plt.xlabel(r'Crystal Momentum $k \cdot (a / \pi)$', fontsize=11)
plt.ylabel('Energy $E$ (eV)', fontsize=11)
plt.xlim(-3, 3)
plt.ylim(0, 100)
plt.axvline(0, color='gray', linestyle=':', alpha=0.7)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=10)

plt.tight_layout()
plt.show()
