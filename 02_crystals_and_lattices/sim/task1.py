#parameters
hbar =1.0
a = 0.5   # lattice constant
t=1.0   # hopping parameter
epsilon =0.0  #On-site energy
N_k = 500

# building bloch hamiltonian
# k-grid covering the First Brillouin Zone: [-π/a, π/a]
k_grid = np.linspace(-np.pi / a, np.pi / a, N_k)

#theoritical effective mass
m_theor = (hbar**2)/(2*t*a**2)

#building bloch hamiltonian
E_k = -2 * t * np.cos(k_grid * a)
#print(E_k)


#curve fit
k_cutoff = 0.10 * (np.pi/a)

mask = np.abs(k_grid) <= k_cutoff
k_fit= k_grid[mask]
E_fit= E_k[mask]

# c2 corresponds to (1/2) * d²E/dk²
coefficients = np.polyfit(k_fit, E_fit, 2)
c2, c1, c0 = coefficients

m_num = (hbar**2) / (2*c2) 


# Output numerical comparison
print(f"--- Effective Mass Extraction Results ---")
print(f"Polynomial fit coefficient c2 (d²E/dk² / 2) : {c2:.6f}")
#print(f"Numerical d²E/dk² at k=0                   : {d2E_dk2_numerical:.6f}")
#print(f"Theoretical d²E/dk² (2*t*a²)               : {2*t*(a**2):.6f}")
print(f"Numerical m*                               : {m_num:.6f}")
print(f"Theoretical m*                             : {m_theor:.6f}")
print(
    f"Relative Error                             :"
    f" {abs(m_num - m_theor)/m_theor * 100:.4f}%"
)


RESULTS
--- Effective Mass Extraction Results ---
Polynomial fit coefficient c2 (d²E/dk² / 2) : 0.248238
Numerical m*                               : 2.014193
Theoretical m*                             : 2.000000
Relative Error                             : 0.7096%



plt.figure(figsize=(8, 5))

# Plot full band structure
plt.plot(
    k_grid * a / np.pi, E_k, 'b-', linewidth=2, label=r'Exact $E(k) = -2t\cos(ka)$'
)

# Plot quadratic fit overlay near k=0
E_fitted_curve = c2 * (k_fit**2) + c1 * k_fit + c0
plt.plot(
    k_fit * a / np.pi,
    E_fitted_curve,
    'r--',
    linewidth=2.5,
    label=rf'Quadratic Fit ($m_{{num}}^* = {m_num:.3f}$)',
)

# Formatting
plt.axvline(0, color='gray', linestyle=':', alpha=0.7)
plt.axhline(epsilon - 2 * t, color='gray', linestyle=':', alpha=0.7)
plt.title('Tight-Binding Band Structure & Effective Mass Fitting', fontsize=12)
plt.xlabel(r'Crystal Momentum $k \cdot (a / \pi)$', fontsize=11)
plt.ylabel('Energy $E(k)$', fontsize=11)
plt.xlim(-1, 1)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=10)

plt.tight_layout()
plt.show()




# Finite difference method

# Parameters
N = 100  # Number of real-space grid points
a = 0.5  # Grid spacing (Å)
hbar = 1.0  # Reduced Planck constant
m_input = 1.0  # Target mass we will simulate

# Pre-factor for finite difference kinetic energy: t_eff = ħ^2 / (2 * m * a^2)
t_eff = (hbar**2) / (2 * m_input * a**2)

# Constructing the real-space tridiagonal Hamiltonian matrix (Periodic Boundary Conditions)
H = np.zeros((N, N))
for i in range(N):
  H[i, i] = 2 * t_eff
  H[i, (i + 1) % N] = -t_eff
  H[i, (i - 1 + N) % N] = -t_eff

# To find E(k), we can apply Bloch phase shifts or diagonalize across k-space
# Alternatively, we can extract m* directly from the lowest eigenvectors / finite difference curvature.
# Let's compute E(k) by evaluating the expectation value of the FD Hamiltonian with Bloch waves:
k_grid = np.linspace(-np.pi / a, np.pi / a, 300)
E_fd = np.zeros_like(k_grid)

for idx, k in enumerate(k_grid):
  # Construct Bloch-periodic Hamiltonian for a single unit cell / phase
  # H_k = 2*t_eff - t_eff*exp(i*k*a) - t_eff*exp(-i*k*a) = 2*t_eff*(1 - cos(k*a))
  E_fd[idx] = 2 * t_eff * (1- np.cos(k * a))
  # Shift to match band bottom at -2*t_eff
  E_fd[idx] -= 2 * t_eff

# Fit curvature near k = 0 to recover m*
k_mask = np.abs(k_grid) < 0.2
k_fit = k_grid[k_mask]
E_fit = E_fd[k_mask]

coeffs = np.polyfit(k_fit**2, E_fit, 1)  # E(k) = c2 * k^2 + c0
c2 = coeffs[0]
m_extracted = (hbar**2) / (2 * c2)

print(f'Input Mass (m*)        : {m_input:.4f}')
print(f'Finite-Difference Fitted m* : {m_extracted:.4f}')

# Plotting the Finite Difference Dispersion
plt.figure(figsize=(7, 4.5))
plt.plot(
    k_grid,
    E_fd,
    color='darkviolet',
    lw=2.5,
    label='Finite-Difference Dispersion $E(k)$',
)
plt.plot(
    k_fit,
    coeffs[1] + coeffs[0] * k_fit**2,
    'r--',
    lw=2,
    label=f'Quadratic Fit ($m^* = {m_extracted:.4f}$)',
)
plt.xlabel('Wavevector $k$ ($1/\AA$)', fontsize=11)
plt.ylabel('Energy $E$ (eV)', fontsize=11)
plt.title(
    'Band Dispersion & Curvature Extracted via Finite-Difference Method',
    fontsize=11,
)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=10)
plt.tight_layout()
plt.show()



### SECOND METHOD ####
N = 500  # Number of grid points
a = 1.0  # Grid spacing in Å
hbar_sq_over_2m0 = 3.810  # eV·Å^2
m_target = 0.5  # Target effective mass in units of m0

# Finite difference prefactor (t_eff = ħ^2 / (2 * m* * a^2))
t_eff = hbar_sq_over_2m0 / (m_target * a**2)

# Correctly construct 1-D arrays for diagonal and off-diagonal elements
diag = 2 * t_eff * np.ones(N)
off_diag = -t_eff * np.ones(N - 1)

# Build the real-space finite-difference matrix
H = np.diag(diag) + np.diag(off_diag, k=1) + np.diag(off_diag, k=-1)

num_levels = 10

energies, wavefunctions = eigh(H, subset_by_index=[0, num_levels - 1])

print(f'--- FIRST {num_levels} ENERGY LEVELS ---')
for i, E in enumerate(energies):
  print(f'E_{i+1} = {E:.4f} eV')
plt.figure(figsize=(8, 5))
for i in range(min(4, num_levels)):
  plt.plot(
      wavefunctions[:, i] + energies[i], label=f'State $n={i+1}$ ($E$={energies[i]:.2f} eV)'
  )

plt.xlabel('Grid Point Index $j$', fontsize=11)
plt.ylabel('Energy / Wavefunction Amplitude (eV)', fontsize=11)
plt.title('Finite-Difference Energy Eigenvalues & Wavefunctions', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=9, loc='upper right')
plt.tight_layout()
plt.show()
