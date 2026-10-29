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
