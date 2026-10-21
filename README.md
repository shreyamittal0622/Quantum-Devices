# Quantum-Mechanics-on-a-Grid

## 1. Scales and Units

For an electron of mass $m^* = \mu\, m_0$ it is convenient to work in nm and eV, where

$$\frac{\hbar^2}$$ $${2m_0}$$ = 0.0381\ \text{eV}\cdot\text{nm}^2, \qquad E_n^{\infty} = \frac{n^2\pi^2\hbar^2}{2m^*L^2}. \tag{1}$$

The second expression is the infinite square well of width $L$. For $m^* = m_0$ and $L = 1$ nm it gives $E_1 \approx 0.376$ eV, an atomic scale; the small effective masses of semiconductors ($\mu = 0.067$ in GaAs) and the larger sizes of lithographic dots bring the scale down to meV.

## 2. Finite-Difference Hamiltonian

On a uniform grid $x_j = ja$ the second derivative becomes $\psi_j'' \approx (\psi_{j+1} - 2\psi_j + \psi_{j-1})/a^2$, so that

$$
(\mathcal{H}\psi)_j = -t\,(\psi_{j+1} + \psi_{j-1}) + (2t + V_j)\,\psi_j, \qquad t = \frac{\hbar^2}{2m^*a^2}. \tag{2}
$$

This is a tridiagonal, real, symmetric matrix. The truncation error is $\mathcal{O}(a^2)$, and the lattice dispersion $E(k) = 2t(1 - \cos ka)$ agrees with the free-particle $\hbar^2k^2/2m^*$ only for $ka \ll 1$; this is the origin of the *lattice* that you will meet physically as the tight-binding model in problem 02.

For a position-dependent effective mass, as at a GaAs/AlGaAs interface, current conservation requires continuity of $\psi$ and of $\psi'/m^*$. The Hermitian discretisation (BenDaniel–Duke) uses a hopping on each bond,

$$
t_{j+\frac12} = \frac{\hbar^2}{2a^2}\cdot\frac12\left(\frac{1}{m_j} + \frac{1}{m_{j+1}}\right), \qquad
(\mathcal{H}\psi)_j = -t_{j+\frac12}\psi_{j+1} - t_{j-\frac12}\psi_{j-1} + \left(t_{j+\frac12} + t_{j-\frac12} + V_j\right)\psi_j. \tag{3}
$$

## 3. Analytic Benchmarks

You need at least four exactly solvable systems to validate the code.

- **Infinite well:** Eq. (1).
- **Harmonic oscillator:** $V = \tfrac12 m^*\omega^2x^2$ gives $E_n = \hbar\omega\left(n + \tfrac12\right)$.
- **Finite well** of depth $V_0$ and width $L$: with $z_0 = \frac{L}{2\hbar}\sqrt{2m^*V_0}$ and $z = \frac{L}{2\hbar}\sqrt{2m^*(E + V_0)}$ (for $E$ measured from the top, $E < 0$), the even and odd bound states satisfy

$$
\tan z = \sqrt{\frac{z_0^2}{z^2} - 1}\ \ (\text{even}), \qquad -\cot z = \sqrt{\frac{z_0^2}{z^2} - 1}\ \ (\text{odd}), \tag{4}
$$

  and the number of bound states is $N = \lceil 2z_0/\pi \rceil$.
- **Attractive delta:** $V = -\alpha\,\delta(x)$ has one bound state, $E = -m^*\alpha^2/2\hbar^2$. On the grid, $V_j = -\alpha/a$ at the central site.

## 4. Time Evolution

The explicit Euler update $\psi \to (1 - i\mathcal{H}\Delta t/\hbar)\psi$ is not unitary and blows up. The Crank–Nicolson scheme,

$$\left(1 + \frac{i\mathcal{H}\Delta t}{2\hbar}\right)\psi^{n+1} = \left(1 - \frac{i\mathcal{H}\Delta t}{2\hbar}\right)\psi^{n}, \tag{5}$$

is unitary, second-order accurate in $\Delta t$, and requires only a tridiagonal solve per step. A free Gaussian packet $\psi_0 \propto \exp\left[-(x - x_0)^2/4\sigma^2 + ik_0x\right]$ moves with group velocity $\hbar k_0/m^*$ and spreads as $\sigma(t) = \sigma\sqrt{1 + \left(\hbar t/2m^*\sigma^2\right)^2}$.

## 5. Scattering from a Barrier

A packet with mean energy $E$ incident on a rectangular barrier of height $V_0$ and width $b$ has, for $E < V_0$,

$$
T(E) = \left[1 + \frac{V_0^2\sinh^2(\kappa b)}{4E(V_0 - E)}\right]^{-1}, \qquad \kappa = \frac{\sqrt{2m^*(V_0 - E)}}{\hbar}. \tag{6}
$$

The transmitted probability of a wavepacket is the energy-weighted average of $T(E)$ over its spectrum, which gives you an independent numerical check of Eq. (6).

