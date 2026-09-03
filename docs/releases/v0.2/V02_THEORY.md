# v0.2 Theory: Multistate Gaussian Dynamics, Adiabatic Couplings, Gauge Continuity, and Dynamic Spawning

Version 0.2 extends the single-surface Gaussian framework to a two-electronic-state
nonadiabatic problem. The central new issue is that the electronic basis itself
depends on nuclear geometry.

Atomic units are used.

---

## 1. Two-state diabatic molecular Hamiltonian

Consider one nuclear coordinate $x$ and two electronic diabatic states:

$$
\boxed{
\hat H_d
=
-\frac{1}{2M}\frac{d^2}{dx^2}I_2
+
V_d(x),
}
$$

with

$$
V_d(x)
=
\begin{pmatrix}
V_{11}(x) & V_{12}(x)\\
V_{12}(x) & V_{22}(x)
\end{pmatrix}.
$$

The exact nuclear wavefunction is a two-component spinor,

$$
\Psi_d(x,t)
=
\begin{pmatrix}
\psi_1(x,t)\\
\psi_2(x,t)
\end{pmatrix}.
$$

The diabatic representation is useful because the kinetic operator is diagonal in
electronic space and nonadiabaticity appears explicitly through $V_{12}$.

---

## 2. Exact two-state split-operator reference

Write

$$
\hat H_d=\hat T I_2+V_d(x).
$$

Second-order Strang splitting gives

$$
\boxed{
e^{-i\hat H_d\Delta t}
=
e^{-iV_d\Delta t/2}
e^{-i\hat T\Delta t}
e^{-iV_d\Delta t/2}
+
\mathcal O(\Delta t^3).
}
$$

At each grid point, $V_d(x_n)$ is a finite Hermitian matrix. Diagonalize

$$
V_d(x_n)U_n=U_n\varepsilon_n
$$

so that

$$
\boxed{
e^{-iV_d(x_n)\Delta t/2}
=
U_n
e^{-i\varepsilon_n\Delta t/2}
U_n^\dagger.
}
$$

The kinetic step is identical for both electronic components and is diagonal after an
FFT.

This calculation is the reference against which all Gaussian approximations in v0.2
are judged.

---

## 3. Adiabatic states

Diagonalize the diabatic potential:

$$
V_d(x)U(x)=U(x)E(x),
$$

where

$$
E(x)=
\begin{pmatrix}
E_0(x)&0\\
0&E_1(x)
\end{pmatrix}.
$$

The columns of $U$ are the adiabatic electronic states expressed in the diabatic
basis.

The adiabatic nuclear amplitudes are

$$
\boxed{
\chi_a(x)=\sum_i U_{ia}^*(x)\psi_i(x),
}
$$

or

$$
\boldsymbol\chi=U^\dagger\boldsymbol\psi_d.
$$

---

## 4. Gauge freedom of adiabatic eigenvectors

An electronic eigenvector is not unique. For a nondegenerate state,

$$
|\phi_a(x)\rangle
\rightarrow
e^{i\theta_a(x)}|\phi_a(x)\rangle
$$

leaves the energy unchanged.

For a real Hamiltonian this reduces locally to arbitrary signs,

$$
|\phi_a\rangle\rightarrow\pm|\phi_a\rangle.
$$

Raw diagonalization at successive geometries can therefore create artificial sign
flips.

For the present real nondegenerate two-state model, v0.2 chooses the sign of each new
eigenvector so that its overlap with the preceding eigenvector is positive:

$$
\boxed{
\langle\phi_a(x_{n-1})|\phi_a(x_n)\rangle>0.
}
$$

This is a simple discrete parallel-transport gauge.

It is **not** sufficient for a genuinely degenerate multidimensional subspace, where
one must align the whole subspace by a unitary transformation.

---

## 5. Derivative couplings

Define

$$
\boxed{
d_{ab}(x)
=
\langle\phi_a(x)|\partial_x\phi_b(x)\rangle.
}
$$

Differentiate the electronic eigenvalue equation. For $a\ne b$,

$$
\boxed{
d_{ab}
=
\frac{
\langle\phi_a|\partial_x V_d|\phi_b\rangle
}{
E_b-E_a
}.
}
$$

For a real orthonormal basis,

$$
d_{ba}=-d_{ab}.
$$

The adiabatic gradients are

$$
\boxed{
E_a'(x)
=
\langle\phi_a|\partial_xV_d|\phi_a\rangle.
}
$$

Both relations are evaluated directly in `adiabatic.py`.

---

## 6. Exact adiabatic nuclear Hamiltonian

The adiabatic representation does not simply replace $V_d$ by diagonal energies.
The kinetic operator differentiates the coordinate-dependent basis.

Write

$$
\Psi=\sum_a\chi_a(x)\phi_a(x).
$$

Projection gives

$$
\boxed{
(H_{\mathrm{ad}}\chi)_a
=
-\frac{1}{2M}
\left[
\chi_a''
+
2\sum_b d_{ab}\chi_b'
+
\sum_b\tau_{ab}\chi_b
\right]
+
E_a\chi_a,
}
$$

where

$$
\tau_{ab}=\langle\phi_a|\partial_x^2\phi_b\rangle.
$$

For a complete electronic subspace,

$$
\boxed{
\tau=d'+d^2.
}
$$

Therefore the kinetic energy can be written compactly as

$$
\boxed{
T_{\mathrm{ad}}
=
-\frac{1}{2M}
(\partial_x I+d)^2.
}
$$

This is a covariant derivative.

The $2d\partial_x$ term is the first-order nonadiabatic kinetic coupling.
The $d'+d^2$ term contains the second-order derivative coupling and diagonal
Born-Oppenheimer correction contributions.

Dropping these terms while calling the calculation "adiabatic propagation" is
incorrect.

---

## 7. Representation equivalence

Let

$$
\Psi_d(x)=U(x)\chi(x).
$$

Then exact operators satisfy

$$
\boxed{
U^\dagger H_d U
=
H_{\mathrm{ad}}
}
$$

when the derivatives acting on $U(x)$ are included.

This identity supplies an important numerical test: apply both operators to the same
smooth wavefunction, transform the diabatic result to the adiabatic basis, and compare
with the covariant-derivative result.

Any discrepancy beyond discretization error indicates a missing derivative-coupling
term or inconsistent gauge.

---

## 8. Adiabatic Gaussian trajectory basis functions

The v0.2 basis state is

$$
\boxed{
|G_{ka}\rangle
=
g_k(x;q_k,p_k,\alpha_k)|\phi_a(x)\rangle.
}
$$

The total molecular wavefunction is approximated by

$$
\boxed{
|\Psi\rangle
=
\sum_{ka}C_{ka}(t)|G_{ka}(t)\rangle.
}
$$

At equal nuclear geometry, different electronic adiabatic states are orthogonal, so

$$
\langle G_{ka}|G_{lb}\rangle
=
\delta_{ab}
\langle g_k|g_l\rangle.
$$

However, the Hamiltonian couples different electronic states through the derivative
coupling terms in $H_{\mathrm{ad}}$.

---

## 9. Projected moving-basis equations

For any time-dependent nonorthogonal basis,

$$
|\Psi\rangle=\sum_\mu C_\mu|G_\mu\rangle.
$$

Projection of the TDSE gives

$$
\boxed{
iS\dot C=(H-i\tau^{\mathrm{basis}})C,
}
$$

where

$$
S_{\mu\nu}=\langle G_\mu|G_\nu\rangle,
$$

$$
H_{\mu\nu}=\langle G_\mu|\hat H|G_\nu\rangle,
$$

and

$$
\tau^{\mathrm{basis}}_{\mu\nu}
=
\langle G_\mu|\dot G_\nu\rangle.
$$

Do not confuse this basis-motion matrix with the electronic second derivative coupling
$\tau_{ab}=\langle\phi_a|\partial_x^2\phi_b\rangle$.

The notation collision is common in the literature; the code names them separately.

---

## 10. Classical motion of adiabatic TBF centers

A trajectory basis function attached to adiabatic state $a$ is guided by

$$
\boxed{
\dot q_k=\frac{p_k}{M},
}
$$

$$
\boxed{
\dot p_k=-E_a'(q_k).
}
$$

This does **not** make the total molecular wavefunction classical. The basis centers
are trajectory guided, but the complex coefficients remain quantum mechanically
coupled and the basis functions interfere.

---

## 11. Dynamic spawning criterion

The purpose of spawning is to enlarge the basis only where the current basis becomes
unable to represent nonadiabatic transfer efficiently.

For a parent on state $a$ and target $b$, define the scalar coupling along the parent
velocity,

$$
\boxed{
\eta_{ab}
=
\left|
\frac{p}{M}d_{ab}(q)
\right|.
}
$$

The v0.2 pedagogical criterion is

$$
\eta_{ab}>\eta_{\mathrm{spawn}}.
$$

A child is not created if a sufficiently overlapping target-state Gaussian already
exists.

This is deliberately simple. Full FMS/AIMS spawning uses more elaborate criteria and
placement optimization.

---

## 12. Child momentum from local energy conservation

For the simple 1D demonstration, a spawned child on state $b$ is placed at the parent
position and assigned momentum from

$$
\frac{p_b^2}{2M}+E_b(q)
=
\frac{p_a^2}{2M}+E_a(q).
$$

Therefore

$$
\boxed{
p_b^2
=
p_a^2
+
2M[E_a(q)-E_b(q)].
}
$$

If the right side is nonnegative,

$$
\boxed{
p_b
=
{\mathrm{sgn}}(p_a)
\sqrt{
p_a^2+2M(E_a-E_b)
}.
}
$$

If it is negative, no real energy-conserving child momentum exists under this
particular placement rule.

This is a pedagogical initialization rule, not a claim to reproduce the full
optimization procedure used by production AIMS implementations.

---

## 13. Why spawning differs from FSSH

In FSSH, a trajectory switches an active label stochastically:

$$
a\rightarrow b.
$$

In spawning, the parent basis function remains and a child basis function is added:

$$
\boxed{
G_{ka}
\rightarrow
G_{ka}+G_{{\mathrm{child}},b}.
}
$$

The child's initial coefficient may be zero; Hamiltonian coupling subsequently moves
amplitude into it.

Thus FMS/AIMS branching is **wavefunction-basis branching**, not trajectory-count
branching.

---

## 14. v0.2 numerical contracts

The implementation is accepted only if:

1. two-state exact split propagation conserves norm;
2. adiabatic energies reproduce eigenvalues of the diabatic potential;
3. derivative couplings are antisymmetric;
4. gauge-aligned adjacent eigenvectors maintain positive overlap;
5. the covariant adiabatic Hamiltonian agrees with the transformed diabatic
   Hamiltonian to finite-difference accuracy;
6. adiabatic Gaussian Hamiltonian matrices are Hermitian within quadrature error;
7. the spawn rule is deterministic;
8. an energy-conserving child satisfies the local energy equation;
9. basis growth leaves the pre-spawn wavefunction unchanged when the new child
   coefficient is initialized to zero.
