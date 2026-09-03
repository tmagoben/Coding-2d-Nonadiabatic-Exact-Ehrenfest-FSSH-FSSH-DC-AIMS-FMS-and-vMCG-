# v0.12 Theory: Representation-Consistent Coherence, Exact LVC Gaussian Propagation, and Born-Huang Projection

Version 0.12 begins with a correction to the interpretation of the v0.10-v0.11
strong-conical-intersection benchmark.

The earlier releases correctly showed that a small center-based Gaussian electronic
representation can disagree strongly with the exact two-dimensional TDSE.

However, v0.12 demonstrates that the discrepancy contained **two different errors**:

1. an **initial electronic-representation error**;
2. a **subsequent dynamics/basis-propagation error**.

Those errors were not cleanly separated before.

The central purpose of v0.12 is therefore

$$
\boxed{
\text{compare like with like before assigning an error to the dynamics}.
}
$$

This release builds three complementary representations:

1. a center-frozen/local-diabatic spinor Gaussian representation;
2. a spinor-complete global-diabatic Gaussian representation;
3. a coordinate-dependent Born-Huang Gaussian representation projected on the exact
   two-dimensional grid.

The first two admit exact analytic LVC Gaussian Hamiltonian matrix elements.
The third reproduces the coordinate dependence of an adiabatic TBF directly and uses
the same FFT kinetic operator as the exact TDSE benchmark.

Atomic units are used throughout.

---

# 1. Exact benchmark state

The exact initial benchmark is

$$
\boxed{
|\Psi_0\rangle
=
g(R;q_0,p_0,A)
|\Phi_{+}(R)\rangle,
}
$$

where

$$
|\Phi_{+}(R)\rangle
$$

is the upper adiabatic electronic eigenstate at the **integration coordinate** $R$.

This distinction matters.

It is not

$$
g(R)
|\Phi_{+}(q_0)\rangle.
$$

The electronic factor varies across the full nuclear wavepacket.

---

# 2. Global diabatic expansion of the adiabatic state

Let

$$
\{|d_1\rangle,|d_2\rangle\}
$$

be the fixed global diabatic electronic basis.

Then

$$
|\Phi_+(R)\rangle
=
\sum_a
U_{a+}(R)|d_a\rangle.
$$

Therefore

$$
|\Psi_0\rangle
=
g(R)
\sum_a U_{a+}(R)|d_a\rangle.
$$

The diabatic nuclear component on electronic basis state $a$ is

$$
\boxed{
\psi_a(R)
=
g(R)U_{a+}(R).
}
$$

---

# 3. Reduced electronic density of a single adiabatic packet

Trace over the nuclear coordinate:

$$
\rho_{ab}
=
\int
\psi_a(R)
\psi_b^*(R)
\,dR.
$$

Substitute:

$$
\boxed{
\rho_{ab}
=
\int
|g(R)|^2
U_{a+}(R)
U_{b+}^*(R)
\,dR.
}
$$

Unless

$$
U_{a+}(R)
$$

is effectively constant over the support of $g$, this matrix is not rank one.

Therefore a wavepacket occupying one adiabatic state can still have

$$
\boxed{
\operatorname{Tr}(\rho_e^2)<1
}
$$

when the electronic reduced density is expressed relative to a fixed global
electronic basis.

This is not numerical decoherence.

It is electron-nuclear correlation arising from the coordinate dependence of the
adiabatic electronic eigenvector.

---

# 4. Center-frozen electronic approximation

The center-frozen approximation replaces

$$
|\Phi_+(R)\rangle
$$

by

$$
|\Phi_+(q_0)\rangle.
$$

Then

$$
|\Psi_{\mathrm{cf}}\rangle
=
g(R)
|\Phi_+(q_0)\rangle.
$$

Its reduced electronic density is

$$
\boxed{
\rho_{\mathrm{cf}}
=
|\Phi_+(q_0)\rangle
\langle\Phi_+(q_0)|.
}
$$

Therefore

$$
\boxed{
\operatorname{Tr}(\rho_{\mathrm{cf}}^2)=1.
}
$$

For the default strong-CI packet this approximation is already visibly different from
the coordinate-dependent exact initial state **before any time propagation occurs**.

That observation changes how the v0.10-v0.11 full-density error should be interpreted.

---

# 5. Correction to the v0.10-v0.11 benchmark interpretation

The v0.10-v0.11 reduced-density helper reconstructed each TBF using an electronic
vector evaluated at the TBF center.

The exact reference used

$$
g(R)\Phi_a(R).
$$

The Gaussian diagnostic used approximately

$$
g(R)\Phi_a(q_i).
$$

These are different ansätze.

Therefore the earlier quantity

$$
\|\rho_{\mathrm{Gaussian}}-\rho_{\mathrm{exact}}\|_F
$$

combined:

- representation error;
- basis incompleteness;
- dynamical error;
- phase/coherence error.

The earlier result remains useful as a diagnostic of the center-frozen model.

It should **not** be interpreted as a clean error of the time propagator alone.

v0.12 explicitly separates these contributions.

---

# 6. Spinor-complete Gaussian ansatz

Introduce a set of nuclear Gaussians

$$
g_k(R)
$$

and attach a complete two-component electronic coefficient vector in the fixed
diabatic basis:

$$
\boxed{
|\Psi_G(R,t)\rangle
=
\sum_k
g_k(R,t)
\sum_{a=1}^2
C_{ka}(t)|d_a\rangle.
}
$$

Equivalently,

$$
\Psi_G(R,t)
=
\sum_k
g_k(R,t)\mathbf C_k(t).
$$

Unlike a one-state-per-TBF basis, every nuclear Gaussian contains the complete
two-state electronic subspace.

This makes electronic representation changes exact within the selected nuclear
Gaussian span.

---

# 7. Spinor-complete overlap matrix

Flatten the compound index

$$
(k,a).
$$

Because the diabatic electronic states are orthonormal and coordinate independent,

$$
\langle d_a|d_b\rangle
=
\delta_{ab}.
$$

Hence

$$
\boxed{
S_{ka,lb}
=
\langle g_k|g_l\rangle
\delta_{ab}.
}
$$

In block form,

$$
\boxed{
S
=
S^{\mathrm{nuc}}\otimes I_2.
}
$$

No derivative coupling or electronic gauge appears in this overlap.

---

# 8. Exact LVC potential matrix element

The repository LVC model is

$$
V_d(x,y)
=
\frac12\omega^2(x^2+y^2)I
+
\kappa x\sigma_z
+
\lambda y\sigma_x.
$$

For Gaussian pair $i,j$, define

$$
S_{ij}
=
\langle g_i|g_j\rangle,
$$

the complex cross centroid

$$
\mu_{ij},
$$

and cross covariance

$$
\Sigma_{ij}
=
(A_i+A_j)^{-1}.
$$

Then

$$
\langle x\rangle_{ij}
=
\mu_x,
$$

$$
\langle y\rangle_{ij}
=
\mu_y,
$$

and

$$
\langle x^2+y^2\rangle_{ij}
=
\mu_x^2+\mu_y^2
+
\Sigma_{xx}+\Sigma_{yy}.
$$

Therefore

$$
\boxed{
\begin{aligned}
\langle g_i|V_d|g_j\rangle
=
S_{ij}
\Big[
&
\frac12\omega^2
(
\mu_x^2+\mu_y^2+
\Sigma_{xx}+\Sigma_{yy}
)I
\\
&+
\kappa\mu_x\sigma_z
+
\lambda\mu_y\sigma_x
\Big].
\end{aligned}
}
$$

For this quadratic/linear LVC Hamiltonian the potential matrix element is exact.

There is no SPA0 or SPA1 truncation in this v0.12 path.

---

# 9. Exact spinor-complete Hamiltonian

The kinetic matrix element is scalar in electronic space:

$$
T_{ka,lb}
=
T_{kl}^{\mathrm{nuc}}\delta_{ab}.
$$

Thus

$$
\boxed{
H_{ka,lb}
=
T_{kl}^{\mathrm{nuc}}\delta_{ab}
+
[V_{kl}]_{ab}.
}
$$

Every interstate coupling is generated directly by the global diabatic potential
matrix.

No explicit derivative coupling is required.

---

# 10. Moving-basis connection in the global diabatic basis

The global diabatic electronic vectors are fixed:

$$
\frac{d}{dt}|d_a\rangle=0.
$$

Therefore

$$
\frac{d}{dt}
[g_j(R,t)|d_b\rangle]
=
\dot g_j(R,t)|d_b\rangle.
$$

Hence

$$
\boxed{
T_{ka,lb}
=
\langle g_k|\dot g_l\rangle
\delta_{ab}.
}
$$

This is substantially cleaner than a center-following adiabatic spinor basis, where an
additional electronic connection appears.

---

# 11. Why an exact Hamiltonian can still give the wrong answer

An exact projected Hamiltonian does not fix an inadequate initial basis.

Suppose the target state is

$$
g(R)\Phi_+(R),
$$

while the spinor-complete Gaussian basis contains only one nuclear Gaussian $g(R)$.

Then the basis can represent only

$$
g(R)
\begin{pmatrix}
c_1\\c_2
\end{pmatrix},
$$

where

$$
c_1,c_2
$$

are constants.

It cannot reproduce the coordinate dependence

$$
U_{a+}(R).
$$

Therefore even exact subsequent propagation begins from a projected initial state that
may differ strongly from the intended exact packet.

This is the central v0.12 diagnosis.

---

# 12. Hilbert-space projection of the exact initial packet

Let the spinor-complete Gaussian basis functions be

$$
|B_{ia}\rangle
=
g_i(R)|d_a\rangle.
$$

We seek

$$
|\Psi_{\mathrm{proj}}\rangle
=
\sum_{ia}C_{ia}|B_{ia}\rangle
$$

that minimizes

$$
\|
\Psi_{\mathrm{target}}
-
\Psi_{\mathrm{proj}}
\|^2.
$$

The normal equations are

$$
\boxed{
SC=b,
}
$$

with

$$
S_{ia,jb}
=
\langle B_{ia}|B_{jb}\rangle
$$

and

$$
\boxed{
b_{ia}
=
\langle B_{ia}|\Psi_{\mathrm{target}}\rangle.
}
$$

v0.12 evaluates $b$ on the same exact two-dimensional grid and solves the linear
least-squares problem.

---

# 13. Projection residual and fidelity

The squared residual is

$$
\boxed{
\epsilon_{\mathrm{proj}}^2
=
\|
\Psi_{\mathrm{target}}
-
\Psi_{\mathrm{proj}}
\|^2.
}
$$

The relative residual is

$$
\boxed{
\epsilon_{\mathrm{rel}}
=
\frac{
\|\Psi_{\mathrm{target}}-\Psi_{\mathrm{proj}}\|^2
}{
\|\Psi_{\mathrm{target}}\|^2
}.
}
$$

The wavefunction fidelity is

$$
\boxed{
F_{\mathrm{proj}}
=
\frac{
|\langle\Psi_{\mathrm{target}}|\Psi_{\mathrm{proj}}\rangle|^2
}{
\langle\Psi_{\mathrm{target}}|\Psi_{\mathrm{target}}\rangle
\langle\Psi_{\mathrm{proj}}|\Psi_{\mathrm{proj}}\rangle
}.
}
$$

These metrics quantify the **initial representation** independently of propagation.

---

# 14. Projected exact dynamics

After forming

$$
\Psi_{\mathrm{proj}}(0),
$$

v0.12 propagates that projected wavefunction with the same exact two-dimensional TDSE
solver:

$$
\boxed{
\Psi_{\mathrm{proj}}^{\mathrm{exact}}(t)
=
e^{-iHt}
\Psi_{\mathrm{proj}}(0).
}
$$

The Gaussian calculation starts from the **same projected state**.

Thus

$$
\boxed{
\epsilon_{\mathrm{dyn}}(t)
=
\|
\rho_G(t)
-
\rho_{\mathrm{proj}}^{\mathrm{exact}}(t)
\|_F
}
$$

measures Gaussian dynamical error without contamination from initial-state projection
error.

This is the most important new validation quantity in v0.12.

---

# 15. Target error

The physically intended exact benchmark remains

$$
\Psi_{\mathrm{target}}^{\mathrm{exact}}(t).
$$

Therefore v0.12 also reports

$$
\boxed{
\epsilon_{\mathrm{target}}(t)
=
\|
\rho_G(t)
-
\rho_{\mathrm{target}}^{\mathrm{exact}}(t)
\|_F.
}
$$

The two errors answer different questions:

- $\epsilon_{\mathrm{dyn}}$: how well does the Gaussian method propagate what it actually
  represents?
- $\epsilon_{\mathrm{target}}$: how well does the complete calculation reproduce the
  intended physical benchmark?

They must not be conflated.

---

# 16. The errors are not additive

It is tempting to write

$$
\epsilon_{\mathrm{target}}
=
\epsilon_{\mathrm{init}}
+
\epsilon_{\mathrm{dyn}}.
$$

That is not generally valid.

Quantum evolution changes the direction of the error vector in Hilbert space, and the
reduced-density map is nonlinear with respect to wavefunction normalization and
partial tracing.

v0.12 therefore reports:

$$
\epsilon_{\mathrm{init}},
\qquad
\epsilon_{\mathrm{dyn}},
\qquad
\epsilon_{\mathrm{target}}
$$

separately.

---

# 17. Projection ladder

The compact release campaign uses three initial nuclear Gaussian banks.

## One Gaussian

One nuclear Gaussian with two electronic components.

This is a useful limiting case but cannot reproduce substantial electronic variation
across the packet.

## Five Gaussians

A center Gaussian plus four shifted Gaussians:

$$
(0,0),
\quad
(\pm0.35,0),
\quad
(0,\pm0.35),
$$

with a narrower coordinate width.

## Nine Gaussians

A $3\times3$ grid of centers with spacing

$$
0.45
$$

and width scale

$$
A_i=3A_0.
$$

These banks are not claimed to be globally optimal.

They form a transparent deterministic convergence ladder.

---

# 18. v0.12 nine-Gaussian result

For the nine-Gaussian reference bank, the initial wavefunction projection gives

$$
\boxed{
F_{\mathrm{proj}}\approx0.832276.
}
$$

The initial reduced-density error is

$$
\boxed{
\epsilon_{\rho,0}\approx0.035455.
}
$$

The Gaussian propagation compared with the exact propagation of the **same projected
initial state** gives

$$
\boxed{
\epsilon_{\mathrm{dyn}}
\approx
2.90\times10^{-4}.
}
$$

This is the key result.

The exact analytic LVC Gaussian propagation is extremely accurate once initial
representation error is removed from the comparison.

---

# 19. Error versus the original exact target

Compared with the original coordinate-dependent adiabatic packet, the nine-Gaussian
reference gives final reduced-density error

$$
\boxed{
\epsilon_{\mathrm{target}}
\approx0.0350003.
}
$$

The diagonal population error is

$$
\boxed{
\epsilon_P
\approx0.0281090.
}
$$

The target purity is

$$
\mathcal P_{\mathrm{exact}}
\approx0.676208,
$$

while the nine-Gaussian result gives

$$
\boxed{
\mathcal P_G
\approx0.662382.
}
$$

The coherence-phase error is only

$$
\boxed{
\Delta\phi_{\rho_{01}}
\approx1.96\times10^{-3}\ {\mathrm{rad}}.
}
$$

Thus v0.12 substantially improves the full electronic density/coherence benchmark.

---

# 20. Relation to v0.11

The v0.11 compact reference gave excellent diagonal populations but failed the
full-density criterion.

Its reported values were approximately

$$
\epsilon_P^{(0.11)}
=
0.01288,
$$

and

$$
\epsilon_\rho^{(0.11)}
=
0.1599.
$$

v0.12 gives a slightly larger population error,

$$
\epsilon_P^{(0.12)}
\approx0.02811,
$$

but reduces the full-density error to

$$
\boxed{
\epsilon_\rho^{(0.12)}
\approx0.03500.
}
$$

The coherence phase is also dramatically more accurate.

This illustrates why diagonal populations alone are not enough to validate
nonadiabatic quantum dynamics.

---

# 21. Coherence observables

For two-state density matrix

$$
\rho
=
\begin{pmatrix}
\rho_{00} & \rho_{01}\\
\rho_{10} & \rho_{11}
\end{pmatrix},
$$

v0.12 reports the complex coherence

$$
\boxed{
c=\rho_{01}.
}
$$

Its magnitude is

$$
|c|.
$$

For candidate and reference coherences,

$$
c_G,
\qquad
c_R,
$$

the phase difference is

$$
\boxed{
\Delta\phi
=
\operatorname{Arg}
(
c_Gc_R^*
).
}
$$

The result is wrapped to

$$
[-\pi,\pi].
$$

A phase is not reported if either coherence magnitude is below a configured floor.

---

# 22. Trace distance

The trace distance between normalized density matrices is

$$
\boxed{
D_{\mathrm{tr}}(\rho,\sigma)
=
\frac12
\|\rho-\sigma\|_1.
}
$$

For Hermitian

$$
\Delta=\rho-\sigma,
$$

this becomes

$$
\boxed{
D_{\mathrm{tr}}
=
\frac12
\sum_k
|\lambda_k(\Delta)|.
}
$$

This provides a basis-independent density-matrix discrepancy measure.

---

# 23. Bloch-vector diagnostic

For a normalized two-state density matrix,

$$
\rho
=
\frac12
(I+\mathbf r\cdot\boldsymbol\sigma).
$$

The implemented Bloch vector is

$$
\boxed{
\mathbf r
=
\begin{pmatrix}
2\operatorname{Re}\rho_{01}\\
-2\operatorname{Im}\rho_{01}\\
\rho_{00}-\rho_{11}
\end{pmatrix}.
}
$$

The $z$ component measures population imbalance.

The $x$ and $y$ components measure electronic coherence.

Therefore Bloch-vector error exposes cases where populations agree but coherence does
not.

---

# 24. Coordinate-dependent Born-Huang TBF benchmark

v0.12 also implements the more literal basis

$$
\boxed{
\Xi_i(R)
=
g_i(R)\Phi_{a_i}(R).
}
$$

This representation reproduces the exact initial localized adiabatic packet with one
TBF.

The overlap is computed directly on the benchmark grid and the Hamiltonian is obtained
by applying the **global diabatic Hamiltonian** to each basis field.

This avoids explicitly inserting divergent first- or second-order adiabatic derivative
couplings.

---

# 25. FFT-projected Hamiltonian

Let

$$
\Xi_i(R)
$$

be represented on the same periodic grid as the exact TDSE.

The exact benchmark Hamiltonian is

$$
\hat H_d
=
-\frac{1}{2M}\nabla^2I
+
V_d(R).
$$

v0.12 evaluates

$$
\boxed{
H_{ij}
=
\langle\Xi_i|\hat H_d|\Xi_j\rangle_{\mathrm{grid}}.
}
$$

The kinetic action is performed by the same FFT momentum lattice used by the exact
split-operator code.

Thus the projection includes:

- the coordinate dependence of the adiabatic eigenvectors;
- gauge/branch structure represented on the grid;
- all kinetic effects generated by that discrete basis field.

This is a benchmark bridge, not a scalable molecular algorithm.

---

# 26. What the Born-Huang path shows

The Born-Huang representation fixes the initial ansatz exactly on the selected grid,
but a small set of classically guided Gaussian centers still does not reproduce the
full exact final branching for the demanding strong-CI passage.

Therefore v0.12 identifies two distinct requirements:

$$
\boxed{
\text{correct electronic representation}
}
$$

and

$$
\boxed{
\text{sufficient nuclear Gaussian span / trajectory motion}.
}
$$

Fixing only one is not enough in every representation.

---

# 27. Why the nine-Gaussian projected global-diabatic path works well

For the release benchmark, the nine-Gaussian spinor-complete bank has enough nuclear
flexibility to represent the electronically varying initial packet to small
reduced-density error.

Once that projection is made, the Hamiltonian is polynomial and integrated exactly
between the Gaussians.

The subsequent projected-state dynamics error is therefore very small.

This does **not** mean nine Gaussians are universally sufficient.

It means that for this specific two-state, two-mode LVC benchmark and propagation
interval, the dominant remaining error is initial representation rather than the exact
Gaussian matrix-element propagation itself.

---

# 28. Dynamic spawning after the nine-Gaussian bank

v0.12 also tests extending the nine-Gaussian bank with additional spawned TBFs.

For this compact benchmark, adding a few spawned functions does not materially reduce
the final target-density error and increases the overlap condition number.

That result is useful:

$$
\boxed{
\text{more adaptive basis functions are not automatically better}.
}
$$

The next function should be added because it resolves a measured missing direction,
not simply because a coupling threshold was crossed.

---

# 29. Implication for PySCF/direct dynamics

A real molecular calculation does not normally provide a global diabatic basis.

The v0.12 lesson is still transferable.

The correct ab initio analogue is to maintain a sufficiently complete electronic
subspace and transport it using many-electron overlaps.

The v0.6-v0.8 infrastructure already provides:

- CASSCF wavefunction snapshots;
- cross-geometry many-electron overlaps;
- state tracking;
- graph gauge transport;
- local subspace alignment.

A future molecular version of the v0.12 representation would combine those overlap
links with a locally complete electronic subspace attached to each Gaussian region.

v0.12 does not claim that this ab initio extension is already implemented.

---

# 30. v0.12 scientific conclusion

The strongest result of v0.12 is not merely a smaller number.

It is a cleaner error decomposition:

$$
\boxed{
\text{initial representation}
\quad\neq\quad
\text{projected-state dynamics}
\quad\neq\quad
\text{target-observable error}.
}
$$

For the nine-Gaussian reference:

$$
\epsilon_{\rho,0}
\approx
3.55\times10^{-2},
$$

$$
\epsilon_{\mathrm{dyn}}
\approx
2.90\times10^{-4},
$$

and

$$
\epsilon_{\mathrm{target}}
\approx
3.50\times10^{-2}.
$$

Therefore the strong-CI benchmark is now limited primarily by how the intended
coordinate-dependent initial electronic state is represented in the finite Gaussian
bank, not by the exact LVC Gaussian coefficient propagation.

That is the correct starting point for v0.13.
