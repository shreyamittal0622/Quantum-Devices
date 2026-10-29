t = 2.7  # eV
acc = 0.142  # nm (carbon-carbon distance)


N = 200
k_max = 4 * np.pi / (3 * np.sqrt(3) * acc)
kx = np.linspace(-k_max, k_max, N)
ky = np.linspace(-k_max, k_max, N)
KX, KY = np.meshgrid(kx, ky)

### Squate Lattice Dispersion
a_sq = 0.246  # Effective square lattice spacing
E_square = -2 * t * (np.cos(KX * a_sq) + np.cos(KY * a_sq))

### Honeycomb Dispersion
f = (
    1
    + 4 * np.cos(np.sqrt(3) * KY * acc / 2) * np.cos(3 * KX * acc / 2)
    + 4 * np.cos(3 * KY * acc / 2) ** 2
)
E_honeycomb_pos = t * np.sqrt(3 + f)
E_honeycomb_neg = -t * np.sqrt(3 + f)

fig = plt.figure(figsize=(12, 5))

# Square Lattice Plot
ax1 = fig.add_subplot(121, projection='3d')
ax1.plot_surface(KX, KY, E_square, cmap='viridis', edgecolor='none', alpha=0.9)
ax1.set_title('Square Lattice Band', fontsize=11)
ax1.set_xlabel('$k_x$', fontsize=9)
ax1.set_ylabel('$k_y$', fontsize=9)
ax1.set_zlabel('Energy (eV)', fontsize=9)

# Honeycomb Lattice Plot (Conduction & Valence Bands)
ax2 = fig.add_subplot(122, projection='3d')
ax2.plot_surface(
    KX, KY, E_honeycomb_pos, cmap='plasma', edgecolor='none', alpha=0.9
)
ax2.plot_surface(
    KX, KY, E_honeycomb_neg, cmap='plasma', edgecolor='none', alpha=0.9
)
ax2.set_title(
    'Honeycomb Lattice (Dirac Cones at $K, K^{\\prime}$)', fontsize=11
)
ax2.set_xlabel('$k_x$', fontsize=9)
ax2.set_ylabel('$k_y$', fontsize=9)
ax2.set_zlabel('Energy (eV)', fontsize=9)

plt.tight_layout()
plt.show()
