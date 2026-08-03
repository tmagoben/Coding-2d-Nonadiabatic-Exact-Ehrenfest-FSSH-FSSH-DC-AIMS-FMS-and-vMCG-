# Rigorous Tutorial: Two-State Nonadiabatic Quantum Dynamics with FFT and Direct Diagonalization

## Purpose

This tutorial builds the simplest complete nonadiabatic quantum-dynamics
workflow in one and two nuclear dimensions. It begins with the static
electronic problem—adiabatic states, derivative couplings, and
adiabatic-to-diabatic transformations—and then propagates a coupled nuclear
wavepacket using two independent algorithms:

1. a Fourier split-operator propagator
2. exact propagation from direct diagonalization of the complete finite-grid
   Hamiltonian.

The second method is computationally expensive but supplies a clean numerical
reference for the first. Every equation uses atomic units, so

$\hbar=e=m_e=4\pi\epsilon_0=1$.

The examples use two electronic states. After reading this derivation, continue to `IMPLEMENTATION_WALKTHROUGH.md`, which maps each equation to the exact Python array operation. All matrix formulas extend to more
states, but the two-state case makes the geometry and numerical implementation
fully visible.

---

# Part I. Molecular separation and the coupled nuclear equation

## 1. Born--Oppenheimer electronic problem

At fixed nuclear coordinates $\(R\)$, the clamped-nuclei electronic Hamiltonian
obeys

$$\hat H_e(r;R)\phi_i(r;R)=E_i(R)\phi_i(r;R)$$.

The electronic coordinates are denoted by $\(r\)$, and $\(R\)$ denotes all nuclear
coordinates. The adiabatic electronic functions are orthonormal at every
geometry:

$\langle\phi_i(R)|\phi_j(R)\rangle_r=\delta_{ij}$

The complete molecular wavefunction is expanded as

$$\Psi(r,R,t)=
\sum_i \chi_i(R,t)\phi_i(r;R).
$$

Substitution into the time-dependent molecular Schrödinger equation produces
coupled nuclear equations. In one nuclear coordinate $\(R\)$, with nuclear mass
$\(M\)$,

$$i\frac{\partial\boldsymbol\chi}{\partial t}=
\left[
-\frac{1}{2M}
\left(
I\frac{\partial^2}{\partial R^2}
+
2\tau(R)\frac{\partial}{\partial R}
+
D(R)
\right)
+
E_{\mathrm{ad}}(R)
\right]\boldsymbol\chi,
$$

where

$$\tau_{ij}(R)=
\left\langle
\phi_i(R)
\middle|
\frac{\partial}{\partial R}
\phi_j(R)
\right\rangle
$$

and

$$
D_{ij}(R)=
\left\langle
\phi_i(R)
\middle|
\frac{\partial^2}{\partial R^2}
\phi_j(R)
\right\rangle.
$$

For a complete electronic basis,

$$D=\frac{\partial\tau}{\partial R}+\tau^2$$.

Consequently, the adiabatic nuclear kinetic operator can be written in a
gauge-covariant form:

$$
\hat T_{\mathrm{ad}}=
-\frac{1}{2M}
\left(
I\frac{\partial}{\partial R}+\tau
\right)^2.
$$

This form is exact but numerically delicate because the derivative couplings
can become sharply peaked near avoided crossings and singular at exact
degeneracies.

---

# Part II. Diabatic representation

## 2. Basis transformation

Let the columns of $\(U(R)\)$ be adiabatic eigenvectors expressed in a fixed
diabatic electronic basis:

$$
|\phi_i(R)\rangle =
\sum_a |\chi_a\rangle U_{ai}(R).
$$

The potential matrices satisfy

$$
U^\mathsf{T}(R)V_{\mathrm{d}}(R)U(R) =
E_{\mathrm{ad}}(R).
$$

If the diabatic basis is independent of nuclear coordinates, its first-order
derivative coupling is zero. The nuclear equation becomes

$$
i\frac{\partial\boldsymbol\psi_{\mathrm{d}}}{\partial t}
\left[
-\frac{1}{2M}\frac{\partial^2}{\partial R^2}I
+
V_{\mathrm{d}}(R)
\right]\boldsymbol\psi_{\mathrm{d}}.
$$

This is the equation propagated in the examples.

The relation between nuclear amplitudes is

$$
\boldsymbol\psi_{\mathrm{d}}(R,t)=
U(R)\boldsymbol\chi_{\mathrm{ad}}(R,t),
$$

and therefore

$$
\boldsymbol\chi_{\mathrm{ad}}(R,t)=
U^\mathsf{T}(R)\boldsymbol\psi_{\mathrm{d}}(R,t)
$$

for a real orthogonal two-state transformation.

The dynamics are representation invariant when every transformation term is
included consistently. We use the diabatic representation because it produces
a simple diagonal kinetic operator and a local matrix-valued potential.

---

# Part III. Static one-dimensional avoided crossing

## 3. Smooth Tully-like single avoided-crossing model

The first model is a fully smooth version of the standard single
avoided-crossing topology:

$$
V_{\mathrm{d}}(x)=\begin{pmatrix}V_{11}(x)&V_{12}(x)\\V_{12}(x)&-V_{11}(x)\end{pmatrix}$$,

where

$$
V_{11}(x)=A\tanh(Bx)
$$

and

$$
V_{12}(x)=Ce^{-Dx^2}.
$$

The hyperbolic tangent gives the same asymptotic scattering structure as the
usual piecewise Tully model but avoids a higher-derivative discontinuity at
the origin. This makes the finite-difference derivative-coupling validation
cleaner without changing the physical lesson.

The adiabatic energies are

$$
E_\pm(x)=
\pm\sqrt{V_{11}^2(x)+V_{12}^2(x)}.
$$

At $\(x=0\)$,

$$
E_+(0)-E_-(0)=2|C|.
$$

Thus $\(C\neq0\)$ converts the diabatic crossing into an adiabatic avoided
crossing.

## 4. Numerical diagonalization and state tracking

At every grid point $\(x_n\)$, we solve

$$
V_{\mathrm{d}}(x_n)U(x_n)=
U(x_n)E_{\mathrm{ad}}(x_n)
$$.

An eigensolver is free to return either $\(\phi_i\)$ or $\(-\phi_i\)$. If the signs
are chosen independently at neighbouring points, the finite-difference
derivative

$$
\frac{\phi_i(x_{n+1})-\phi_i(x_{n-1})}{2\Delta x}
$$

contains artificial sign discontinuities of order $\(1/\Delta x\)$.

The implementation therefore performs:

1. overlap calculation,
   $$
   O_{ij}=
   \left\langle\phi_i(x_{n-1})|\phi_j(x_n)\right\rangle
   $$
2. state assignment maximizing $\(|O_{ij}|\)$
3. sign correction enforcing positive diagonal overlap.

## 5. Derivative couplings

For nondegenerate states,

$$
\tau_{ij}(x)=
\frac{
\langle\phi_i|\partial_xV|\phi_j\rangle
}{
E_j-E_i
},
\qquad i\ne j.
$$

The code independently evaluates

$$
\tau(x)=U^\mathsf{T}(x)\frac{dU(x)}{dx}
$$

with finite differences. Agreement tests both the derivative of the
Hamiltonian and phase tracking.

For a real orthonormal basis,

$$
\tau^\mathsf{T}=-\tau.
$$

The diagonal elements are gauge dependent and are set to zero in the chosen
real parallel-transport gauge.

## 6. Pathwise adiabatic-to-diabatic transformation

Define

$$
|\chi_a(x)\rangle=
\sum_i|\phi_i(x)\rangle A_{ia}(x).
$$

Demanding vanishing derivative coupling in the transformed basis gives

$$
\frac{dA}{dx}=-\tau(x)A(x).
$$

The exact formal solution is a path-ordered exponential:

$$
A(x)=\mathcal P\exp\left[-\int_{x_0}^x\tau(x')\,dx'\right]A(x_0).
$$

Numerically, the interval is divided into steps and the midpoint approximation
is used:

$$
A_{n+1}=
\exp\left[-\frac{\tau_n+\tau_{n+1}}{2}\Delta x\right]A_n$$.

Because $\(\tau\)$ is antisymmetric itself, every exponential is orthogonal. The
recovered diabatic matrix is

$$
V_{\mathrm{d}}^{\mathrm{rec}}(x)=
A^\mathsf{T}(x)
E_{\mathrm{ad}}(x)
A(x).
$$

The initial condition $\(A(x_0)\)$ fixes a constant diabatic gauge. In a real
ab initio calculation, different valid initial choices produce matrices
related by a constant orthogonal rotation.

Run:

```bash
python examples/01_static_1d_derivative_coupling_and_adt.py
```

---

# Part IV. Nuclear grid and Fourier representation

## 7. Uniform periodic grid

Choose \(N\) coordinate points over a box of length $\(L=x_{\max}-x_{\min}\)$:

$$
x_n=x_{\min}+n\Delta x,
\qquad
\Delta x=\frac{L}{N},
\qquad
n=0,\ldots,N-1.
$$

The periodic Fourier momenta returned by the discrete FFT are

$$
k_m=
2\pi\,\mathrm{fftfreq}(N,\Delta x).
$$

The largest representable momentum is approximately the Nyquist value

$$
|k|_{\max}\approx\frac{\pi}{\Delta x}.
$$

A wavepacket whose momentum distribution reaches the Nyquist boundary is
aliased and is not numerically resolved.

The periodic grid also implies that probability leaving one side of the box
re-enters from the opposite side. Every example therefore stops before the
wavepacket reaches the boundaries. A production scattering calculation would
normally add a complex absorbing potential or mask function and then validate
the absorber independently.

## 8. Fourier kinetic operator

The one-dimensional nuclear kinetic energy is

$$
\hat T=-\frac{1}{2M}\frac{\partial^2}{\partial x^2}.
$$

For a plane wave,

$$
\hat T e^{ikx}=
\frac{k^2}{2M}e^{ikx}.
$$

Therefore the kinetic propagator is diagonal in momentum space:

$$
e^{-i\hat T\Delta t}\tilde\psi(k)=
e^{-ik^2\Delta t/(2M)}\tilde\psi(k).
$$

The algorithm is:

1. Fourier transform $\(\psi(x)\rightarrow\tilde\psi(k)\)$
2. multiply by $\(e^{-ik^2\Delta t/(2M)}\)$
3. inverse Fourier transform.

The cost is

$
O(N\log N)
$

per electronic state and time step.

---

# Part V. Exact local potential propagation

## 9. Matrix exponential at every coordinate

The diabatic potential is a \(2\times2\) matrix. It must not be propagated by
independently exponentiating its elements.

Write a real symmetric matrix as

\[
V=v_0I+h_x\sigma_x+h_z\sigma_z,
\]

where

\[
v_0=\frac{V_{11}+V_{22}}{2},
\qquad
h_z=\frac{V_{11}-V_{22}}{2},
\qquad
h_x=V_{12}.
\]

Define

\[
q=\sqrt{h_x^2+h_z^2}.
\]

Using

\[
(h_x\sigma_x+h_z\sigma_z)^2=q^2I,
\]

the exact local exponential is

\[
e^{-iV\Delta t}
=
e^{-iv_0\Delta t}
\left[
\cos(q\Delta t)I
-i\frac{\sin(q\Delta t)}{q}
(h_x\sigma_x+h_z\sigma_z)
\right].
\]

At \(q=0\), the stable limiting expression is

\[
\frac{\sin(q\Delta t)}{q}\rightarrow\Delta t.
\]

The code constructs this unitary matrix once for a time-independent
potential.

---

# Part VI. FFT split-operator dynamics

## 10. Why splitting is required

The full propagator is

\[
e^{-i(\hat T+\hat V)\Delta t}.
\]

Since

\[
[\hat T,\hat V]\ne0,
\]

we cannot write this exactly as

\[
e^{-i\hat T\Delta t}e^{-i\hat V\Delta t}.
\]

The symmetric Strang factorization is

\[
e^{-i(\hat T+\hat V)\Delta t}
=
e^{-i\hat V\Delta t/2}
e^{-i\hat T\Delta t}
e^{-i\hat V\Delta t/2}
+
O(\Delta t^3)
\]

for one time step. Over a fixed final time, the global error is

\[
O(\Delta t^2).
\]

## 11. One complete step

For the two-component diabatic wavefunction

\[
\boldsymbol\psi(x,t)
=
\begin{pmatrix}
\psi_1(x,t)\\
\psi_2(x,t)
\end{pmatrix},
\]

one step is:

1. apply half the local potential evolution,
   $\boldsymbol\psi^{(1)}=e^{-iV(x)\Delta t/2}\boldsymbol\psi(t)$

2. FFT both electronic components,

   $\tilde{\boldsymbol\psi}^{(1)}(k)=\mathcal F[\boldsymbol\psi^{(1)}(x)]$

3. apply the full kinetic phase,

   $\tilde{\boldsymbol\psi}^{(2)}(k)=e^{-ik^2\Delta t/(2M)}\tilde{\boldsymbol\psi}^{(1)}(k)$

4. inverse FFT,
   $\boldsymbol\psi^{(2)}(x)=\mathcal F^{-1}[\tilde{\boldsymbol\psi}^{(2)}(k)]$

5. apply the second half potential step,
   $\boldsymbol\psi(t+\Delta t) = e^{-iV(x)\Delta t/2}\boldsymbol\psi^{(2)}(x)$.

Both subpropagators are unitary, so norm conservation should be close to
machine precision. Norm conservation alone does not prove convergence:
Strang splitting can conserve norm while still having phase and population
errors. Time-step convergence must be checked.

---

# Part VII. Initial wavepacket and observables

## 12. Gaussian nuclear wavepacket

The one-dimensional scalar packet is

$$g(x)=\mathcal N\exp\left[-\frac{(x-x_0)^2}{4\sigma^2}+ip_0(x-x_0)\right]$$.

Its probability density has variance $\(\sigma^2\)$:

$$|g(x)|^2\propto\exp\left[-\frac{(x-x_0)^2}{2\sigma^2}\right]$$.

The discrete normalization is


$$\Delta x\sum_n|g(x_n)|^2=1$$.

To prepare the packet on adiabatic state $\(i\)$, convert it to diabatic
components:

$\psi_a(x,0)=U_{ai}(x)g(x)$.

This construction is valid when the packet is localized in a region where the
chosen adiabatic state is smooth.

## 13. Diabatic populations

The diabatic state populations are

$$P_a^{\mathrm{d}}(t)=\int|\psi_a(x,t)|^2dx$$.

On the grid,
$$
P_a^{\mathrm{d}}(t)\approx\Delta x\sum_n|\psi_a(x_n,t)|^2$$.

## 14. Adiabatic populations

At every coordinate,

$$\boldsymbol\chi_{\mathrm{ad}}(x,t)=U^\mathsf{T}(x)\boldsymbol\psi_{\mathrm{d}}(x,t)$$.

The adiabatic population is

$P_i^{\mathrm{ad}}(t)=\int|\chi_i(x,t)|^2dx$.

These populations reveal transfer between the local Born--Oppenheimer
surfaces. They are not generally identical to diabatic populations.

---

# Part VIII. Direct diagonalization

## 15. Finite-grid Hamiltonian matrix

The FFT kinetic operation defines a finite periodic spectral matrix

\[
$T=
F^\dagger
\operatorname{diag}\left(
\frac{k_m^2}{2M}
\right)
F$,

where $\(F\)$ is the discrete Fourier transform matrix.

For two electronic states, flatten the wavefunction in state-major order:

$$
\boldsymbol\Psi=\begin{pmatrix}\psi_1(x_0)\\
\vdots\\
\psi_1(x_{N-1})\\
\psi_2(x_0)\\
\vdots\\
\psi_2(x_{N-1})
\end{pmatrix}.
$$

The full matrix is


$$H=\begin{pmatrix}
T+\operatorname{diag}(V_{11})&\operatorname{diag}(V_{12})\\
\operatorname{diag}(V_{21})&T+\operatorname{diag}(V_{22})
\end{pmatrix}$$.

Its dimension is $\(2N\times2N\)$.

## 16. Exact finite-basis propagator

Diagonalize

$HW=W\varepsilon$.

Because $\(H\)$ is Hermitian, $\(W\)$ is unitary. The exact propagator within this
finite periodic basis is

$$
e^{-iHt}=
W e^{-i\varepsilon t}W^\dagger$$.

Thus

$$
\boldsymbol\Psi(t)=
W e^{-i\varepsilon t}W^\dagger
\boldsymbol\Psi(0)$$.


There is no time-step error. There is still finite-grid error, finite-box
error, and periodic-boundary error.

Dense diagonalization scales approximately as

$
O((2N)^3)
$

in time and

$
O((2N)^2)
$

in memory. It is therefore a reference method, not a large-grid production
method.

## 17. Comparing the two propagators

The normalized wavefunction fidelity is

$$
\mathcal F(t)=
\frac{
|\langle\Psi_{\mathrm{direct}}(t)|
\Psi_{\mathrm{FFT}}(t)\rangle|^2
}{
\langle\Psi_{\mathrm{direct}}|\Psi_{\mathrm{direct}}\rangle
\langle\Psi_{\mathrm{FFT}}|\Psi_{\mathrm{FFT}}\rangle
}$$.


Because a global phase is physically irrelevant, the code also removes the
best global phase before computing an $\(L^2\)$ error.

Run:

```bash
python examples/02_dynamics_1d_fft_vs_diagonalization.py
```

---

# Part IX. Two-dimensional conical intersection

## 18. Linear vibronic-coupling model

The two-dimensional diabatic potential is

$$
V_{\mathrm{d}}(x,y)=
\frac{k(x^2+y^2)}{2}I+
\begin{pmatrix}
\kappa x&\lambda y\\
\lambda y&-\kappa x
\end{pmatrix}.
$$

The adiabatic surfaces are

$$
E_\pm(x,y)=
\frac{k(x^2+y^2)}{2}
\pm
\sqrt{\kappa^2x^2+\lambda^2y^2}.
$$

They are exactly degenerate at

$x=y=0$.

The coordinates $\(x\)$ and $\(y\)$ form the idealized branching plane.

## 19. Derivative-coupling vector field

Define

$$
\theta(x,y)=
\operatorname{atan2}(\lambda y,\kappa x).
$$

In a local real gauge,

$$
\boldsymbol\tau_{01}=
\frac12\nabla\theta.
$$

Therefore

$$
\tau_{01}^{(x)}=
-\frac{\kappa\lambda y}
{2(\kappa^2x^2+\lambda^2y^2)}
$$

and

$$
\tau_{01}^{(y)}=
\frac{\kappa\lambda x}
{2(\kappa^2x^2+\lambda^2y^2)}$$.

The field diverges at the conical intersection. This singularity is a feature
of the adiabatic representation, not a divergence in the smooth diabatic
potential matrix.

## 20. Berry phase

For a closed loop enclosing the origin,
$$
\oint\boldsymbol\tau_{01}\cdot d\mathbf R=
\pm\pi.
$$

Parallel transport around the loop gives

$$
A_{\mathrm{final}}\approx-A_{\mathrm{initial}}.
$$

The potential matrix is unchanged by this sign, but a nuclear wavefunction
on one adiabatic sheet acquires a geometric phase. A globally smooth,
single-valued real adiabatic basis cannot cover a region enclosing the
intersection without a branch cut or equivalent gauge construction.

Run:

```bash
python examples/03_static_2d_conical_intersection.py
```

---

# Part X. Two-dimensional FFT dynamics

## 21. Two-dimensional kinetic operator

With masses $\(M_x\)$ and $\(M_y\)$,

$$
\hat T=-\frac{1}{2M_x}\frac{\partial^2}{\partial x^2}-\frac{1}{2M_y}\frac{\partial^2}{\partial y^2}.
$$

In Fourier space,


$$
T(k_x,k_y)=
\frac{k_x^2}{2M_x}
+
\frac{k_y^2}{2M_y}.
$$

The kinetic propagator is

$$
e^{-iT(k_x,k_y)\Delta t}.
$$

The two-dimensional split step is exactly the same five-stage sequence as in
one dimension, replacing the FFT by a two-dimensional FFT.

For $\(N_xN_y\)$ spatial points, the split-operator cost scales approximately as

$$
O(N_xN_y\log(N_xN_y)).
$$

## 22. Two-dimensional direct benchmark

The spatial kinetic matrix is a Kronecker sum:


$$T_{2D}=I_y\otimes T_x+T_y\otimes I_x$$.

For two electronic states, the complete dimension is

$$
2N_xN_y.
$$

Dense diagonalization therefore scales as

$$
O((2N_xN_y)^3),
$$

which becomes prohibitive very quickly. The example performs:

1. a meaningful $\(64\times64\)$ FFT propagation;
2. a separate $\(14\times14\)$ direct-diagonalization benchmark.

The tiny direct benchmark verifies the implementation but is not asserted to
be a converged physical grid.

Run:

```bash
python examples/04_dynamics_2d_fft_and_direct_benchmark.py
```

---

# Part XI. Numerical convergence protocol

A scientifically useful calculation must vary each independent numerical
parameter.

## 23. Spatial spacing

Repeat with increasing $\(N\)$ at fixed box length. Check:

$$
\Delta P_i,
\qquad
\Delta\langle x\rangle,
\qquad
1-\mathcal F.
$$

The initial and evolved momentum distributions must remain well inside the
Nyquist interval.

## 24. Box size

Increase $\(x_{\max}-x_{\min}\)$ while maintaining similar spacing. Confirm that
the packet does not interact with its periodic image.

## 25. Time step

For Strang splitting, repeat with

$$
\Delta t,\quad\Delta t/2,\quad\Delta t/4.
$$

An observable $\(O(\Delta t)\)$ should approach

$$
O(\Delta t)=O(0)+C\Delta t^2+O(\Delta t^4).
$$

A useful Richardson estimate is

$$
O(0)
\approx
\frac{4O(\Delta t/2)-O(\Delta t)}{3}.
$$

## 26. Direct-diagonalization comparison

At a grid size where dense diagonalization is possible, compare:

- wavefunction fidelity;
- phase-aligned $\(L^2\)$ error;
- diabatic populations;
- adiabatic populations;
- norm.

This isolates the split-operator time-step error from the finite-grid error.

## 27. Static electronic validation

Before propagating:

1. verify $\(U^\mathsf{T}U=I\)$;
2. verify derivative-coupling antisymmetry;
3. compare Hellmann--Feynman and finite-difference couplings;
4. reconstruct the known diabatic potential from the ADT;
5. test local potential-propagator unitarity.

---

# Part XII. Mapping the tutorial to ab initio data

The model potential can be replaced by electronic-structure data generated at
each nuclear geometry.

An ab initio workflow supplies:

$$
E_i(R),
\qquad
\nabla_R E_i(R),
\qquad
\tau_{ij}(R),
$$

together with consistent state tracking and phase alignment.

For a one-dimensional path:

1. perform state-averaged electronic-structure calculations;
2. order states by overlap rather than energy alone;
3. phase-align neighboring states;
4. calculate derivative couplings;
5. integrate
   $$
   dA/dR=-\tau A;
   $$
6. form
   $$
   V_{\mathrm{d}}=A^\mathsf{T}E_{\mathrm{ad}}A;
   $$
7. interpolate the smooth matrix elements;
8. propagate the nuclear wavepacket with the same FFT code.

For a multidimensional system, a pathwise transformation always exists away
from exact degeneracies, but a global diabatic representation may be
path-dependent. The Berry phase in the 2D example is the simplest explicit
demonstration of this obstruction.

---

# Part XIII. Code map

## `nonadiabatic_dynamics/models.py`

Defines the 1D avoided crossing and 2D conical-intersection potentials and
their analytic derivatives.

## `nonadiabatic_dynamics/adiabatic.py`

Performs adiabatic diagonalization, overlap-based state tracking, phase
correction, Hellmann--Feynman derivative couplings, and finite-difference
couplings.

## `nonadiabatic_dynamics/diabatization.py`

Integrates the path-ordered ADT with midpoint matrix exponentials and
reconstructs the diabatic potential.

## `nonadiabatic_dynamics/grids.py`

Constructs periodic coordinate grids and FFT momentum grids.

## `nonadiabatic_dynamics/wavepackets.py`

Builds normalized Gaussian packets and embeds them on local adiabatic states.

## `nonadiabatic_dynamics/propagators.py`

Constructs exact local \(2\times2\) potential exponentials and performs 1D and
2D Strang split-operator steps.

## `nonadiabatic_dynamics/direct.py`

Builds the complete Fourier-spectral finite-grid Hamiltonian, diagonalizes it,
and applies the exact finite-basis propagator.

## `nonadiabatic_dynamics/observables.py`

Calculates norm, diabatic populations, adiabatic populations, fidelity, and
phase-aligned wavefunction errors.

---

# Part XIV. Recommended sequential execution

```bash
python -m pip install -e .

python examples/01_static_1d_derivative_coupling_and_adt.py
python examples/02_dynamics_1d_fft_vs_diagonalization.py
python examples/03_static_2d_conical_intersection.py
python examples/04_dynamics_2d_fft_and_direct_benchmark.py
python examples/05_convergence_study_1d.py

pytest
```

The first and third scripts establish the electronic representation. The
second and fourth scripts then propagate the nuclear dynamics. The tests use
smaller grids to verify unitarity and agreement between FFT propagation and
direct diagonalization.
