# v0.13 PySCF Residual-Driven Basis Bridge

The v0.13 release benchmark uses an analytic two-state LVC model because the full
target wavefunction and exact TDSE residual can be evaluated on a compact nuclear grid.

A molecular PySCF calculation does not provide those objects directly.

This document states the molecular analogue without overclaiming.

## 1. Existing PySCF infrastructure

The repository already contains:

```text
v0.5  explicit SA-CASSCF energies, gradients, and NAC adapter
v0.6  many-electron cross-geometry overlaps and state tracking
v0.7  graph/subspace gauge transport
v0.8  incremental dynamic electronic graph
v0.9  conditioning and convergence machinery
v0.10-v0.12 exact low-dimensional validation layers
```

These remain the foundation for an ab initio residual method.

## 2. What is unavailable in full-dimensional direct dynamics

For a realistic molecule, one generally cannot evaluate

$$
r(R)
=
\Psi_{\rm target}(R)-P\Psi_{\rm target}(R)
$$

over the full

$$
3N-6
$$

dimensional nuclear configuration space.

Likewise, evaluating

$$
\mathcal R(R,t)
=
i\dot\Psi-H\Psi
$$

globally is not tractable.

Therefore the literal 2D grid residual is a benchmark object, not the proposed
full-dimensional implementation.

## 3. Sampled residual analogue

A practical molecular analogue is to define a quadrature or sample set

$$
\{R_k,w_k\}.
$$

At each sampled geometry, construct a tracked local electronic subspace using the
PySCF many-electron overlap machinery.

Then evaluate a discrete residual norm

$$
\boxed{
\|r\|_{\rm sample}^2
=
\sum_k
w_k
\|r(R_k)\|^2.
}
$$

Candidate Gaussian gains can be evaluated on the same sample set.

The result is a controlled approximation to the global Hilbert residual.

## 4. Electronic transport is essential

At different geometries,

$$
\Phi_I(R_k)
$$

do not have independent meaningful phases or state labels.

Before subtracting or comparing electronic amplitudes, the states must be transported
into a common local gauge/subspace.

The v0.6-v0.8 overlap graph provides the required structure.

A molecular residual computed before state/subspace alignment would mix physical error
with arbitrary electronic gauge choices.

## 5. Complete local electronic blocks

v0.12 showed why one-state-per-Gaussian center representations can obscure coherence.

The molecular analogue should therefore use a selected local electronic subspace

$$
\mathcal S(R)
=
\operatorname{span}\{
\Phi_1(R),\ldots,\Phi_m(R)
\}
$$

and transport the complete coefficient vector between nearby electronic frames.

Near degeneracies, subspace Procrustes/polar transport is preferred over independent
root sign matching.

## 6. Candidate score

On a sampled molecular support, a candidate nuclear Gaussian

$$
g_c(R)
$$

can still be orthogonalized against the existing Gaussian span.

The discrete residual-gain score becomes

$$
\boxed{
\Delta_c^{\rm sample}
=
\frac{
\sum_a
\left|
\sum_k
w_k
g_c^\perp(R_k)^*
r_a(R_k)
\right|^2
}{
\sum_k
w_k
|g_c^\perp(R_k)|^2
}.
}
$$

This is the direct sampled analogue of the v0.13 grid formula.

## 7. Electronic-structure cost control

A molecular implementation should not launch a new electronic-structure calculation
for every candidate in a large dictionary.

A realistic workflow would:

1. generate candidate nuclear Gaussians cheaply;
2. rank candidates using cached/interpolated local electronic information;
3. evaluate new PySCF points only for a small residual shortlist;
4. update the overlap/gauge graph;
5. accept a candidate only if it reduces the measured residual enough to justify its
   cost.

This connects residual selection with the basis-growth efficiency problem that
motivates modern AIMS basis-control methods.

## 8. Dynamic defect analogue

The exact full-dimensional defect is unavailable.

Possible tractable surrogates include:

- residuals on a local quadrature around TBF centroids;
- residuals at Gaussian pair saddle points;
- mismatch between SPA0/SPA1 matrix-element predictions;
- discrepancy between local-diabatic overlap propagation and explicit NAC propagation;
- projected coefficient-equation residuals across a sampled electronic graph.

These should be validated first on analytic models where the full defect is known.

## 9. What v0.13 does not claim

v0.13 does not yet implement:

- full-dimensional molecular TDSE residual evaluation;
- automatic PySCF calls during residual candidate search;
- production residual-triggered AIMS spawning;
- residual-driven pruning in a molecular run;
- a universal molecular error estimator.

The contribution of v0.13 is the audited low-dimensional mathematics and the software
interfaces needed to test those ideas before connecting them to costly electronic
structure.
