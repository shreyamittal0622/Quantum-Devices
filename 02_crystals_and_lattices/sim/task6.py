# III-V Semiconductor Data (Vurgaftman et al. 0K values)
materials = ['InSb', 'InAs', 'GaSb', 'InP', 'GaAs']
Eg = np.array([0.235, 0.417, 0.812, 1.424, 1.519])  # Band gap in eV
m_star_measured = np.array([0.014, 0.026, 0.039, 0.079, 0.063])  # Measured m*/m0

# Theoretical Curve: m0 / m* = Ep / Eg (Assuming Ep = 28.8 eV)
Ep = 28.8
Eg_fine = np.linspace(0.2, 1.6, 100)
m_star_theory = 1 / (Ep / Eg_fine)

# Plotting
plt.figure(figsize=(9, 5.5))

# Theoretical Trend Line
plt.plot(Eg_fine,m_star_theory,color='crimson',linestyle='--',linewidth=2,label='Two-Band Theory (m0/m* = Ep/Eg)',)

# Measured Data Points
plt.scatter(Eg,m_star_measured, color='royalblue', s=90, zorder=5, label='Measured (Vurgaftman et al.)',)

# Annotate each material point
for i, mat in enumerate(materials):
  plt.annotate(mat,(Eg[i], m_star_measured[i]),textcoords='offset points',xytext=(0, 10),ha='center',fontsize=10,fontweight='bold',)

plt.title('Effective Mass vs. Band Gap in III-V Semiconductors', fontsize=12)
plt.xlabel('Band Gap $E_g$ (eV)', fontsize=11)
plt.ylabel('Effective Mass $m^* / m_0$', fontsize=11)
plt.xlim(0.1, 1.7)
plt.ylim(0, 0.09)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=10, loc='upper left')

plt.tight_layout()
plt.show()
