
hbar_eV_A = 3.81  # ħ in eV·Å (with m0 scaled out: ħ^2 / (2 m0) = 3.81 eV·Å^2)

# --- 1. Synthetic Band Data (Simulated numerical output near k=0) ---
m_true = 0.067  # True effective mass (GaAs-like)
k_pts = np.linspace(-0.04, 0.04, 21)  # k in Å^-1
E_true = (hbar_eV_A / m_true) * k_pts**2  # E(k) in eV
# Add tiny numerical noise
np.random.seed(42)
E_data = E_true + np.random.normal(0, 1e-4, size=len(k_pts))

# --- 2. Fit quadratic polynomial: E(k) = c2 * k^2 + c0 ---
# Fit against k^2 directly
poly_coeffs = np.polyfit(k_pts**2, E_data, 1)  # returns [c2, c0]
c2 = poly_coeffs[0]

# --- 3. Extract m* ---
m_fitted = hbar_eV_A / c2

print(f'Fitted c2 coefficient : {c2:.4f} eV·Å^2')
print(f'Extracted m*          : {m_fitted:.4f} m0')
print(f'True m*               : {m_true:.4f} m0')

# --- 4. Plotting Verification ---
plt.figure(figsize=(7, 4.5))
plt.plot(k_pts, E_data, 'ob', label='Numerical Data Points')
plt.plot(
    k_pts,
    poly_coeffs[0] * k_pts**2 + poly_coeffs[1],
    'r--',
    lw=2,
    label=f'Quadratic Fit ($m^* = {m_fitted:.3f} m_0$)',
)
plt.xlabel('Wavevector $k$ ($\AA^{-1}$)')
plt.ylabel('Energy $E$ (eV)')
plt.title('Curvature Fitting at Band Bottom')
plt.grid(True, alpha=0.4)
plt.legend()
plt.tight_layout()
plt.show()
