# v0.12 PySCF / Ab Initio Bridge

v0.12's main result comes from an analytic two-state LVC model with a known global
diabatic basis.

Real molecular calculations generally do not provide that luxury.

The purpose of this document is to state precisely what can and cannot be transferred
to the existing PySCF infrastructure.

---

## 1. What already exists in the repository

Earlier releases already provide:

```text
v0.5  PySCF SA-CASSCF energies / gradients / NACs
v0.6  cross-geometry many-electron state overlaps and state tracking
v0.7  gauge graph / unitary transport / Wilson loops
v0.8  incremental time-dependent electronic gauge graph
v0.9  basis conditioning and convergence controls
v0.10 exact-reference and reduced-density benchmark machinery
v0.11 width-diverse adaptive Gaussian basis
```

These pieces remain unchanged and regression tested.

---

## 2. What the analytic global diabatic basis did for v0.12

For the LVC model,

$$
\{|d_a\rangle\}
$$

is fixed globally.

Therefore every Gaussian can carry a complete electronic vector

$$
\mathbf C_i
$$

without derivative couplings or phase conventions.

The exact matrix is simply

$$
H_{ia,jb}
=
T_{ij}\delta_{ab}
+
\langle g_i|V_{ab}(R)|g_j\rangle.
$$

There is no exact analogue of this globally fixed basis supplied automatically by
PySCF for an arbitrary molecule.

---

## 3. Ab initio analogue: locally complete electronic subspaces

The molecular analogue should be a selected electronic subspace

$$
\mathcal S(R)
=
\operatorname{span}
\{
\Phi_1(R),\ldots,\Phi_m(R)
\}.
$$

At neighboring geometries, the many-electron overlap matrix is

$$
\boxed{
O_{IJ}
=
\langle
\Phi_I(R_A)
|
\Phi_J(R_B)
\rangle.
}
$$

v0.6 already computes this for restricted SA-CASSCF snapshots.

The overlap matrix is therefore the natural object for transporting the locally
complete electronic vector.

---

## 4. Why individual root tracking is insufficient for v0.12-style coherence

For population-only trajectory guidance, one may identify one persistent adiabatic
root.

For a complete electronic coefficient vector, especially near a degeneracy, the
entire subspace matters.

If

$$
\Phi_B
=
\Phi_AU,
$$

then coefficients must transform as

$$
\boxed{
c_B
=
U^\dagger c_A.
}
$$

Near an exact degeneracy, no unique individual-state labeling exists.

The v0.7 Procrustes/subspace transport machinery is therefore more appropriate than
forcing a scalar root label.

---

## 5. Discrete local-diabatic propagation

A practical molecular v0.12 analogue would use electronic overlaps rather than
explicitly integrating singular derivative couplings.

For neighboring electronic frames,

$$
O_{n,n+1}
=
\Phi_n^\dagger\Phi_{n+1}.
$$

Take its unitary polar factor

$$
L_{n,n+1}.
$$

Then an electronic vector can be transported discretely across the graph.

This is closely aligned with the local-diabatic philosophy already present in v0.7
and v0.8.

---

## 6. Matrix elements between Gaussian regions

The difficult molecular quantity is not merely state transport.

One needs electronic Hamiltonian information between nuclear Gaussian basis regions.

The release does **not** pretend that

$$
\langle g_i\Phi_I|\hat H|g_j\Phi_J\rangle
$$

is available exactly from standard single-geometry PySCF calls.

Possible controlled approximations include:

1. zeroth-order saddle-point electronic Hamiltonians;
2. first-order Taylor matrix elements;
3. local-diabatization using state overlaps at centroids;
4. quadrature over cached electronic points in a low-dimensional reaction coordinate;
5. fitted diabatic Hamiltonian models.

The existing SPA0/SPA1 infrastructure covers the first two educational layers.

---

## 7. Representation-consistent molecular initial conditions

v0.12 shows that this step must precede dynamics validation.

If the intended initial state is

$$
g(R)\Phi_a(R),
$$

a center-frozen approximation

$$
g(R)\Phi_a(R_0)
$$

can already alter the reduced electronic density.

For a molecular calculation, one should therefore evaluate whether the initial
electronic subspace varies appreciably over the nuclear Wigner/packet support.

Possible diagnostics:

```text
many-electron state overlaps across the initial packet width
principal angles between selected electronic subspaces
overlap unitarity defect
change in local electronic density / character
```

If these vary strongly, one center electronic calculation is not a controlled
representation of the initial packet.

---

## 8. Proposed molecular workflow

A rigorous next-stage PySCF calculation would proceed:

```text
1. Define SA-CASSCF electronic contract
        ↓
2. Sample / discretize initial nuclear support
        ↓
3. Evaluate PySCF snapshots at those points
        ↓
4. Build cross-geometry many-electron overlaps
        ↓
5. Construct a locally complete tracked electronic subspace
        ↓
6. Quantify initial representation error
        ↓
7. Build local-diabatic Gaussian matrix approximations
        ↓
8. Propagate the coupled Gaussian coefficients
        ↓
9. Refine electronic points / TBFs independently
```

This prevents an electronic-representation error from being misdiagnosed as a nuclear
propagation failure.

---

## 9. What v0.12 does not claim

v0.12 does not yet implement:

- exact ab initio Gaussian Hamiltonian integrals;
- a globally diabatic molecular PySCF basis;
- production AIMS forward/backward optimal spawning;
- a multidimensional molecular Born-Huang quadrature;
- an automatic active-space selection procedure.

The current PySCF backend remains a rigorous source of local electronic information
and many-electron overlaps.

The v0.12 analytic benchmark establishes the representation/validation requirements
that a future molecular implementation should satisfy.
