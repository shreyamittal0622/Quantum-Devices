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
