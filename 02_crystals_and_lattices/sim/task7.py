m0 = 9.1093837e-31  # Free electron mass (kg)
hbar = 1.0545718e-34  # Reduced Planck's constant (J s)
eV = 1.60217663e-19  # Joules per eV
meV = 1e-3 * eV

m_star = 0.067 * m0  # GaAs effective mass
L = 10e-9  # Confinement width (10 nm)


# Calculate fundamental confinement energy unit E1 in meV
# E1 = (hbar^2 * pi^2) / (2 * m_star * L^2)
E1_J = (hbar**2 * np.pi**2) / (2 * m_star * L**2)
E1_meV = E1_J / meV
print(f'Ground State Confinement Energy E1 = {E1_meV:.2f} meV')

# Energy grid for DOS evaluation (0 to 400 meV)
E_vals = np.linspace(0.1, 400, 2000)


# A. Quantum Well (2D DOS) -> Step functions
dos_2d = np.zeros_like(E_vals)
for n in range(1, 15):
  En = n**2 * E1_meV
  dos_2d += np.heaviside(E_vals - En, 1.0)

# B. Quantum Wire (1D DOS) -> Inverse square-root spikes
dos_1d = np.zeros_like(E_vals)
for ny in range(1, 10):
  for nz in range(1, 10):
    En = (ny**2 + nz**2) * E1_meV
    diff = E_vals - En
    valid = diff > 0
    contrib = np.zeros_like(E_vals)
    contrib[valid] = 1.0 / np.sqrt(diff[valid])
    dos_1d += np.minimum(contrib, 15.0)  # Cap numerical singularity peaks

# C. Quantum Box / Dot (0D DOS) -> Lorentzian broadened delta peaks
eta = 3.0  # Broadening parameter in meV
dos_0d = np.zeros_like(E_vals)
for nx in range(1, 6):
  for ny in range(1, 6):
    for nz in range(1, 6):
      En = (nx**2 + ny**2 + nz**2) * E1_meV
      if En <= 400:
        dos_0d += (1.0 / (np.pi * eta)) / (1.0 + ((E_vals - En) / eta) ** 2)

fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

# 2D Quantum Well Plot
axes[0].plot(E_vals, dos_2d, color='royalblue', lw=2)
axes[0].set_title(
    'Quantum Well (2D) Density of States', fontsize=12, fontweight='bold'
)
axes[0].set_ylabel('DOS (arb. units)', fontsize=10)
axes[0].grid(True, linestyle='--', alpha=0.5)

# 1D Quantum Wire Plot
axes[1].plot(E_vals, dos_1d, color='darkorange', lw=2)
axes[1].set_title(
    'Quantum Wire (1D) Density of States', fontsize=12, fontweight='bold'
)
axes[1].set_ylabel('DOS (arb. units)', fontsize=10)
axes[1].grid(True, linestyle='--', alpha=0.5)

# 0D Quantum Box Plot
axes[2].plot(E_vals, dos_0d, color='forestgreen', lw=2)
axes[2].set_title(
    'Quantum Box / Dot (0D) Density of States', fontsize=12, fontweight='bold'
)
axes[2].set_xlabel('Energy E (meV)', fontsize=11)
axes[2].set_ylabel('DOS (arb. units)', fontsize=10)
axes[2].grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.show()
