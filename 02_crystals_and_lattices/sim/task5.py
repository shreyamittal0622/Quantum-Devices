# Parameters
t = 1.0  # Hopping parameter (eV)
eta = 0.05  # Lorentzian broadening factor
hbar = 1
m_star = 1
E_grid1 = np.linspace(-6.2, 6.2, 400)
#E_grid = np.linspace(0.005, 3.0, 400)
# --- 1. 1D Chain DOS ---
N_1d = 5000
k_1d = np.linspace(-np.pi, np.pi, N_1d)
E_1d = -2 * t * np.cos(k_1d)

dos_1d_num = np.zeros_like(E_grid)
for e in E_1d:
  dos_1d_num += (1 / np.pi) * (eta / ((E_grid - e) ** 2 + eta**2))
dos_1d_num /= N_1d

dos_1d_theory = (1 / (np.pi * hbar)) * np.sqrt(2 * m_star / E_grid)


#-----2. 2D chain DOS -----
N_2d = 300
kx_2d = np.linspace(-np.pi, np.pi, N_2d)
ky_2d = np.linspace(-np.pi, np.pi, N_2d)
KX, KY = np.meshgrid(kx_2d, ky_2d)

E_2d = -2 * t * (np.cos(KX) + np.cos(KY))
E_2d_flat = E_2d.flatten()

dos_2d_num = np.zeros_like(E_grid)
for e in E_2d_flat:
  dos_2d_num += (1 / np.pi) * (eta / ((E_grid - e) ** 2 + eta**2))
dos_2d_num /= len(E_2d_flat)

dos_2d_theory = np.full_like(E_grid, (m_star / (np.pi * hbar**2)) / t)

# --- 3. 3D Simple Cubic Lattice DOS ---
N_3d = 60
kx_3d = np.linspace(-np.pi, np.pi, N_3d)
ky_3d = np.linspace(-np.pi, np.pi, N_3d)
kz_3d = np.linspace(-np.pi, np.pi, N_3d)
K_X, K_Y, K_Z = np.meshgrid(kx_3d, ky_3d, kz_3d, indexing='ij')
E_3d = -2 * t * (np.cos(K_X) + np.cos(K_Y) + np.cos(K_Z))
E_3d_flat = E_3d.flatten()

dos_3d_num = np.zeros_like(E_grid)
for e in E_3d_flat:
  dos_3d_num += (1 / np.pi) * (eta / ((E_grid - e) ** 2 + eta**2))
dos_3d_num /= len(E_3d_flat)

dos_3d_theory = ( (1 / (2 * np.pi**2)) * ((2 * m_star / hbar**2) ** 1.5) * np.sqrt(E_grid))


fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

## plots analytical graphs###
# 1D Plot
axes[0].plot(E_grid1, dos_1d_num, color='royalblue', linewidth=2)
axes[0].set_title('1D Chain DOS', fontsize=11)
axes[0].set_xlabel('Energy $E$ (eV)', fontsize=10)
axes[0].set_ylabel('$g(E)$', fontsize=10)
axes[0].grid(True, linestyle='--', alpha=0.5)

# 2D Plot
axes[1].plot(E_grid1, dos_2d_num, color='darkorange', linewidth=2)
axes[1].set_title('2D Square Lattice DOS', fontsize=11)
axes[1].set_xlabel('Energy $E$ (eV)', fontsize=10)
axes[1].grid(True, linestyle='--', alpha=0.5)

# 3D Plot
axes[2].plot(E_grid1, dos_3d_num, color='forestgreen', linewidth=2)
axes[2].set_title('3D Cubic Lattice DOS', fontsize=11)
axes[2].set_xlabel('Energy $E$ (eV)', fontsize=10)
axes[2].grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.show()



## plots a comparison between theoritcal and numerical values##
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# 1D Comparison
axes[0].plot(E_grid,dos_1d_num,color='royalblue',linewidth=2, label='Numerical (1D)',)
axes[0].plot( E_grid, dos_1d_theory,color='crimson',linestyle='--',linewidth=2,label='Eq. (4) [$\\propto E^{-1/2}$]',)
axes[0].set_title('1D Band-Edge DOS', fontsize=12)
axes[0].set_xlabel('Energy from Band Bottom (eV)', fontsize=10)
axes[0].set_ylabel('$g(E)$', fontsize=10)
axes[0].set_ylim(0, 5)
axes[0].grid(True, linestyle='--', alpha=0.5)
axes[0].legend(fontsize=9)

# 2D Comparison
axes[1].plot(E_grid,  dos_2d_num,color='darkorange',linewidth=2,  label='Numerical (2D)',)
axes[1].plot(E_grid,dos_2d_theory,color='crimson',linestyle='--',linewidth=2,label='Eq. (4) [Constant]',)
axes[1].set_title('2D Band-Edge DOS', fontsize=12)
axes[1].set_xlabel('Energy from Band Bottom (eV)', fontsize=10)
axes[1].set_ylim(0, 0.5)
axes[1].grid(True, linestyle='--', alpha=0.5)
axes[1].legend(fontsize=9)

# 3D Comparison
axes[2].plot(E_grid,dos_3d_num,color='forestgreen',linewidth=2,label='Numerical (3D)',)
axes[2].plot(E_grid,dos_3d_theory,color='crimson',linestyle='--',linewidth=2,label='Eq. (4) [$\\propto \\sqrt{E}$]',)
axes[2].set_title('3D Band-Edge DOS', fontsize=12)
axes[2].set_xlabel('Energy from Band Bottom (eV)', fontsize=10)
axes[2].set_ylim(0, 0.3)
axes[2].grid(True, linestyle='--', alpha=0.5)
axes[2].legend(fontsize=9)

plt.tight_layout()
plt.show()
