a = 1.0
t1 = 1.0
t2 = [0.4, 1.0, 1.6]  # Different t2/t1 ratios
N_k = 300
k_grid = np.linspace(-np.pi / a, np.pi / a, N_k)

plt.figure(figsize=(10,5))

for i,t2 in enumerate(t2):
    ratios = t2/ t1
    E_minus = np.zeros(N_k)
    E_plus = np.zeros(N_k)

    for j ,k in enumerate(k_grid):
        off_diag = t1 +t2*np.exp(-1j*k*a)
        H_k = np.array([[0,off_diag],[np.conj(off_diag),0]])

        evals,_ = np.linalg.eigh(H_k)
        
        E_minus[j]=evals[0]
        E_plus[j]=evals[1]

        
    plt.subplot(1, 2, 1)
    plt.plot(
      k_grid * a / np.pi, E_plus, label=f'$t_2/t_1 = {t2}$ (Conduction)' )
    plt.plot(
      k_grid * a / np.pi, E_minus, linestyle='--', color='gray', alpha=0.7)

plt.subplot(1, 2, 1)
plt.title('SSH Band Structures for Various $t_2/t_1$', fontsize=11)
plt.xlabel(r'Crystal Momentum $k \cdot (a / \pi)$', fontsize=10)
plt.ylabel('Energy $E(k)$', fontsize=10)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=8)


# band gaps#
t2_range = np.linspace(0.1, 2.0, 100) # took for values of t2 of a different range 
band_gaps = 2 * np.abs(t1 - t2_range)

plt.subplot(1, 2, 2)
plt.plot(t2_range / t1, band_gaps, 'r-', linewidth=2.5)
plt.title('Band Gap vs. Bond Alternation ($t_2 / t_1$)', fontsize=11)
plt.xlabel(r'Hopping Ratio $t_2 / t_1$', fontsize=10)
plt.ylabel('Band Gap $\Delta$', fontsize=10)
plt.axvline(1.0, color='black', linestyle=':', label='Gap Closes ($t_1 = t_2$)')
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=9)

plt.tight_layout()
plt.show()


##band gap for t2 given ###
# Parameters (took for values of t2 given to us)
t1 = 1.0
t2_values = [0.4, 1.0, 1.6]

# Calculate band gaps for your discrete list
# For the SSH model, the band gap is analytically given by Delta = 2 * abs(t1 - t2)
band_gaps = [2 * abs(t1 - t2) for t2 in t2_values]
ratios = [t2 / t1 for t2 in t2_values]

# Plotting discrete points
plt.figure(figsize=(6, 4))
plt.plot(
    ratios,
    band_gaps,
    'ro-',
    linewidth=2,
    markersize=8,
    label='Calculated Points',
)

plt.title('Band Gap at Specific $t_2 / t_1$ Values', fontsize=11)
plt.xlabel(r'Hopping Ratio $t_2 / t_1$', fontsize=10)
plt.ylabel('Band Gap $\Delta$', fontsize=10)
plt.axvline(1.0, color='black', linestyle=':', label='Gap Closes ($t_1 = t_2$)')
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=9)

plt.tight_layout()
plt.show()
