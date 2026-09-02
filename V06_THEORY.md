# v0.6 Theory: Electronic-State Tracking, Gauge Continuity, and Many-Electron Overlap

Version 0.6 addresses a problem that becomes unavoidable once the v0.5 PySCF backend
is used along a trajectory:

> the electronic root returned as "state 0" at one geometry need not represent the
> same physical electronic state as "state 0" at the next geometry.

Near avoided crossings, conical-intersection regions, or root reordering, simple
energy sorting can destroy state identity. Independent eigenvectors can also change
sign even when the state itself is perfectly smooth.

The purpose of v0.6 is to make those issues explicit and provide a reproducible
tracking/gauge layer.

Atomic units are used for dynamics.

---

# 1. Why energy ordering is not state tracking

At geometry $\mathbf R_n$, suppose SA-CASSCF returns

$$
E_0(\mathbf R_n)<E_1(\mathbf R_n).
$$

At the next geometry the physical characters may exchange while the solver again
sorts roots by energy.

Then

$$
\boxed{
\text{raw root index}
\neq
\text{persistent electronic-state identity}.
}
$$

For a dynamics algorithm, this matters because state labels are attached to

$$
E_I,
\qquad
\nabla E_I,
\qquad
\mathbf d_{IJ},
\qquad
\text{TBF electronic labels}.
$$

A root permutation left untreated can turn a smooth surface into a discontinuous
one.

---

# 2. Electronic phase/gauge freedom

A normalized adiabatic electronic state is not unique:

$$
\boxed{
|\Phi_I(\mathbf R)\rangle
\rightarrow
e^{i\theta_I(\mathbf R)}
|\Phi_I(\mathbf R)\rangle.
}
$$

For the real RHF/ROHF SA-CASSCF calculations targeted here,

$$
e^{i\theta_I}
\rightarrow
s_I,
\qquad
s_I=\pm1.
$$

The energy and gradient do not change under this sign transformation.

The derivative coupling does:

$$
\mathbf d_{IJ}
=
\langle\Phi_I|\nabla\Phi_J\rangle.
$$

If

$$
|\Phi_I'\rangle=s_I|\Phi_I\rangle,
$$

then

$$
\boxed{
\mathbf d'_{IJ}
=
s_Is_J\mathbf d_{IJ}.
}
$$

Therefore state tracking must transform the NAC tensor whenever it changes electronic
signs.

---

# 3. Cross-geometry many-electron overlap

The natural continuity diagnostic is

$$
\boxed{
O_{IJ}^{(n,n+1)}
=
\langle
\Psi_I(\mathbf R_n)
|
\Psi_J(\mathbf R_{n+1})
\rangle.
}
$$

Here $\Psi_I$ is the complete restricted CASSCF electronic wavefunction represented
by the doubly occupied core orbitals and the active-space CI expansion.

The state at the new geometry with the largest many-electron overlap with a previous
state is the natural candidate for the same tracked identity.

---

# 4. Nonorthogonal orbital bases at different geometries

Atomic orbitals centered on different nuclear geometries are not the same one-particle
basis.

Define the cross-AO overlap

$$
\boxed{
S_{\mu\nu}^{AB}
=
\langle
\chi_\mu(\mathbf R_A)
|
\chi_\nu(\mathbf R_B)
\rangle.
}
$$

PySCF supplies this object through its cross-molecule one-electron integral interface.

If the core+active MO coefficient matrices are

$$
C_A,
\qquad
C_B,
$$

the cross-MO overlap is

$$
\boxed{
S_{\mathrm{MO}}^{AB}
=
C_A^\dagger
S_{\mathrm{AO}}^{AB}
C_B.
}
$$

This matrix is generally not the identity.

---

# 5. CASSCF wavefunction in core + active orbitals

For restricted CASSCF, the core spatial orbitals are doubly occupied in every
configuration.

Write schematically

$$
|\Psi_I\rangle
=
\sum_{KL}
C^{(I)}_{KL}
|
K_\alpha L_\beta;
\text{core doubly occupied}
\rangle.
$$

The active-space CI coefficient tensor alone is not enough for an exact
cross-geometry overlap if the previous core subspace overlaps with the current active
subspace.

A simple factorization such as

$$
(\det S_{\rm core})^2
\times
\langle\Psi_{\rm CAS}^A|\Psi_{\rm CAS}^B\rangle
$$

assumes away those core-active cross blocks.

v0.6 does **not** make that approximation.

---

# 6. Exact core+active CI embedding

Let

$$
n_c
$$

be the number of doubly occupied core orbitals and

$$
n_a
$$

the number of active orbitals.

The correlated one-particle space used for overlap tracking is ordered as

$$
\boxed{
[
\text{core}_1,\ldots,\text{core}_{n_c},
\text{active}_1,\ldots,\text{active}_{n_a}
].
}
$$

Every active alpha determinant is extended by occupying all core alpha orbitals.
Every active beta determinant is extended in the same way.

Thus each active-space configuration is embedded into a determinant of the full
core+active orbital space.

The embedded CI vector has zero amplitude for every determinant that violates the
CASSCF core-occupation constraint.

This produces the exact restricted-CASSCF wavefunction within the occupied
core+active orbital space.

---

# 7. Determinant overlap in a nonorthogonal orbital basis

For two Slater determinants built from nonorthogonal spin orbitals,

$$
|\mathcal D_A\rangle,
\qquad
|\mathcal D_B\rangle,
$$

the overlap is the determinant of the occupied-orbital overlap matrix:

$$
\boxed{
\langle\mathcal D_A|\mathcal D_B\rangle
=
\det S_{\rm occ}^{AB}.
}
$$

For a spin-separated alpha/beta representation,

$$
\boxed{
\langle
D_{\alpha A}D_{\beta A}
|
D_{\alpha B}D_{\beta B}
\rangle
=
\det S_\alpha
\det S_\beta.
}
$$

Therefore the CASSCF state overlap is

$$
\boxed{
O_{IJ}
=
\sum_{KL,MN}
(C^{A,I}_{KL})^*
C^{B,J}_{MN}
\det S_{\alpha;KM}
\det S_{\beta;LN}.
}
$$

v0.6 evaluates this through PySCF's FCI overlap implementation after embedding the
active CI vectors into the full core+active determinant space.

---

# 8. Why the full embedding matters

Suppose a previous core orbital has nonzero overlap with a current active orbital.

Then the cross-MO overlap contains blocks

$$
S^{AB}_{ca}\neq0,
\qquad
S^{AB}_{ac}\neq0.
$$

A core-times-active factorization discards these terms.

The full determinant overlap does not.

That difference is explicitly unit tested in v0.6 with a two-orbital model where the
exact overlap is changed by core-active orbital rotation.

---

# 9. Global root assignment

For $N$ states, define

$$
O_{ij}
=
\langle
\Psi_i^{\rm previous}
|
\Psi_j^{\rm current}
\rangle.
$$

A greedy assignment can fail if two states both have substantial overlap with the same
new root.

v0.6 instead considers a permutation

$$
\pi:\{0,\ldots,N-1\}\rightarrow\{0,\ldots,N-1\}
$$

and maximizes

$$
\boxed{
\mathcal S(\pi)
=
\sum_i
|O_{i,\pi(i)}|^2.
}
$$

For the small state manifolds used in these educational/direct-dynamics examples, an
exhaustive permutation search is transparent and deterministic.

The result

$$
\pi(i)
$$

means:

> tracked state $i$ is represented by raw current root $\pi(i)$.

---

# 10. Phase/sign correction

After assignment, let

$$
z_i
=
O_{i,\pi(i)}.
$$

For a complex electronic calculation, the phase that makes the overlap positive real
would be

$$
\boxed{
p_i
=
\frac{z_i^*}{|z_i|}.
}
$$

Then

$$
p_iz_i=|z_i|.
$$

The current v0.6 PySCF backend is deliberately restricted to the real RHF/ROHF
SA-CASSCF case.

Therefore it uses

$$
\boxed{
p_i=\operatorname{sign}(z_i)\in\{+1,-1\}.
}
$$

An appreciably complex assigned overlap is treated as incompatible with that
real-gauge contract rather than silently discarded.

---

# 11. Transformation of state-resolved properties

Let raw current roots be transformed as

$$
\boxed{
|\Phi_i^{\rm tracked}\rangle
=
p_i
|\Phi_{\pi(i)}^{\rm raw}\rangle.
}
$$

Then

$$
\boxed{
E_i^{\rm tracked}
=
E_{\pi(i)}^{\rm raw},
}
$$

and

$$
\boxed{
\nabla E_i^{\rm tracked}
=
\nabla E_{\pi(i)}^{\rm raw}.
}
$$

The NAC tensor transforms as

$$
\boxed{
\mathbf d_{ij}^{\rm tracked}
=
p_i^*p_j
\mathbf d_{\pi(i),\pi(j)}^{\rm raw}.
}
$$

For real signs,

$$
p_i^*p_j=p_ip_j.
$$

The CI vectors stored for the next overlap calculation are also reordered and
multiplied by the same phases.

Thus the selected electronic gauge is propagated recursively along the path.

---

# 12. Tracking ambiguity

Maximum overlap does not guarantee that the electronic-state identity is physically
well defined.

Two failure signals are used.

## 12.1 Small assigned overlap

If

$$
\boxed{
\min_i |O_{i,\pi(i)}|
<
O_{\min},
}
$$

the previous tracked state has poor overlap with every plausible current root.

Possible reasons include:

- geometry step too large;
- active-space character change;
- root loss;
- insufficient state manifold;
- severe orbital/state reorganization.

## 12.2 Nearly tied assignments

Let

$$
\mathcal S_1
$$

and

$$
\mathcal S_2
$$

be the best and second-best global assignment scores.

If

$$
\boxed{
\mathcal S_1-\mathcal S_2
<
\Delta_{\rm assign},
}
$$

the state mapping is intrinsically ambiguous at the selected thresholds.

The default v0.6 policy is to raise rather than silently guess.

---

# 13. Exact degeneracy and subspace identity

At an exact degeneracy, individual adiabatic eigenvectors are not unique.

Any unitary rotation inside the degenerate space produces an equally valid electronic
basis:

$$
\boxed{
\Phi
\rightarrow
\Phi U.
}
$$

In this case the scientifically meaningful object is the subspace projector

$$
\boxed{
P=\Phi\Phi^\dagger,
}
$$

not the individual roots.

For overlap matrix block

$$
O_{\mathcal A\mathcal B},
$$

the singular values

$$
\boxed{
\sigma_k
}
$$

measure principal overlap between the previous and current subspaces.

Values near one mean the subspace is continuous even if the individual states rotate
strongly.

v0.6 records singular-value diagnostics and exposes subspace-overlap tools rather than
pretending maximum-overlap root labels remain unique exactly at a degeneracy.

---

# 14. Unitary Procrustes / local diabatic transport

Let

$$
O
=
\Phi_{\rm prev}^\dagger
\Phi_{\rm curr}.
$$

Take the SVD

$$
\boxed{
O=U\Sigma V^\dagger.
}
$$

We seek a rotation $Q$ of the current subspace that maximizes its alignment with the
previous one:

$$
\Phi_{\rm curr}'
=
\Phi_{\rm curr}Q.
$$

The Procrustes solution is

$$
\boxed{
Q=VU^\dagger.
}
$$

Then

$$
OQ
=
U\Sigma U^\dagger,
$$

which is Hermitian positive semidefinite.

This is the discrete local-diabatization/parallel-transport operation implemented in

```text
overlap_transport.py
```

It is useful for degenerate or near-degenerate subspace analysis.

---

# 15. Overlap-derived directional derivative coupling

Suppose the electronic states vary along scalar path coordinate $s$.

For a small displacement $\Delta s$,

$$
|\Phi_j(s+\Delta s)\rangle
=
|\Phi_j(s)\rangle
+
\Delta s
\frac{\partial|\Phi_j\rangle}{\partial s}
+
\mathcal O(\Delta s^2).
$$

Therefore

$$
O_{ij}
=
\langle\Phi_i(s)|\Phi_j(s+\Delta s)\rangle
$$

becomes

$$
\boxed{
O_{ij}
=
\delta_{ij}
+
\Delta s\,d_{ij}^{(s)}
+
\mathcal O(\Delta s^2),
}
$$

where

$$
d_{ij}^{(s)}
=
\left\langle
\Phi_i
\middle|
\frac{\partial\Phi_j}{\partial s}
\right\rangle.
$$

Since $d$ is anti-Hermitian,

$$
\boxed{
d^{(s)}
\approx
\frac{
O-O^\dagger
}{
2\Delta s
}.
}
$$

v0.6 implements this as a local diagnostic.

It reconstructs only the derivative coupling **along the sampled path direction**.
It does not replace the full Cartesian NAC vector returned by PySCF.

---

# 16. Overlap unitarity defect

If the tracked electronic subspace were complete and the geometries infinitesimally
close, its overlap matrix would approach a unitary transformation.

For a finite selected subspace, define

$$
\boxed{
\epsilon_U
=
\|O^\dagger O-I\|_F.
}
$$

A large value indicates loss of overlap outside the selected state manifold and/or a
large geometry step.

v0.6 records this diagnostic for every tracked PySCF step.

---

# 17. Geometric phase is not removed by tracking

At a conical intersection, no globally single-valued real adiabatic gauge exists
around a loop enclosing the degeneracy.

Local maximum-overlap sign tracking can keep adjacent states continuous, but after one
closed circuit it may produce

$$
\boxed{
|\Phi(2\pi)\rangle=-|\Phi(0)\rangle.
}
$$

That is not a failure of the tracker.

It is the molecular geometric phase already demonstrated explicitly in v0.4.

Therefore:

$$
\boxed{
\text{local gauge continuity}
\neq
\text{global removal of Berry phase}.
}
$$

---

# 18. Warm-start orbitals versus state tracking

The v0.5 PySCF backend optionally reuses the preceding CASSCF orbital matrix as an
initial guess.

This improves optimization continuity but does not solve root identity.

Warm starting acts on the orbital optimization.

State tracking acts on the many-electron state labels and phases.

They are complementary:

$$
\boxed{
\text{warm-start MO continuity}
+
\text{many-electron root overlap}
}
$$

is substantially more robust than either alone.

---

# 19. Sequential path semantics

The tracked PySCF backend is stateful.

Each call is compared with the **previous accepted geometry**:

```text
R0 -> R1 -> R2 -> R3
```

This is scientifically meaningful for a single trajectory or geometry scan.

It is not automatically meaningful for calls in arbitrary order such as

```text
pair centroid A
pair centroid C
trajectory point B
pair centroid A again
```

because the reference state would depend on call order.

v0.6 therefore makes the class name and documentation explicitly path based.

---

# 20. Why this matters for spawned Gaussian dynamics

In an FMS/AIMS-style calculation, multiple TBFs and pair centroids may request
electronic data.

There are three increasingly sophisticated strategies.

## Strategy 1 — one tracking history per trajectory

Each TBF center carries its own sequential electronic-state history.

This is natural for trajectory guidance.

## Strategy 2 — precomputed tracked scan

For a one-dimensional reaction coordinate, evaluate a PySCF scan sequentially once,
track the roots, then save

$$
E_I(q),
\quad
E_I'(q),
\quad
d_{IJ}(q).
$$

The saved tracked scan can then be interpolated in arbitrary order.

v0.6 implements this route.

## Strategy 3 — multidimensional gauge graph

For general multidimensional direct dynamics, electronic states should be connected
through a graph of overlap relations and transported consistently over that graph.

That problem includes loop holonomy/geometric phase and is intentionally left for a
future release.

---

# 21. Tracked 1D scan provider

For ordered points

$$
q_0<q_1<\cdots<q_N,
$$

v0.6 can build a tracked table

$$
\{
E_I(q_n),
\nabla_qE_I(q_n),
d_{IJ}^{(q)}(q_n),
M_q(q_n)
\}.
$$

After tracking is complete, a separate interpolation provider is constructed.

That provider is order independent because the state identities have already been
fixed by the sequential overlap calculation.

This is the recommended v0.6 route for connecting PySCF state tracking to the
gridless 1D Gaussian machinery.

---

# 22. v0.6 scientific contract

The tracking layer is accepted only if:

1. global maximum-overlap assignment recovers known root swaps;
2. electronic sign flips are corrected;
3. NAC signs transform consistently;
4. ambiguous assignments are detected;
5. exact core+active determinant embedding reproduces known many-electron overlaps;
6. core-active cross overlap is retained;
7. subspace Procrustes transport recovers a known unitary rotation;
8. overlap-derived directional NAC converges to the analytic generator for a small
   step;
9. tracked scan interpolation is deterministic and order independent;
10. the entire v0.1-v0.5 regression suite remains passing.

The resulting framework is still not a production AIMS package, but it removes one of
the most important electronic-structure consistency failures that would otherwise make
a PySCF-driven nonadiabatic trajectory unreliable.
