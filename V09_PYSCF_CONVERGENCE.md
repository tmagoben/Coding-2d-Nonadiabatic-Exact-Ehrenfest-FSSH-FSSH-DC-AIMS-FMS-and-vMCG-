# v0.9 PySCF Convergence Workflow

v0.5--v0.8 established an explicit PySCF SA-CASSCF backend, many-electron state
tracking, and graph gauge transport.  v0.9 adds a **convergence protocol** for using
that backend in a direct-dynamics calculation.

PySCF is optional and is not required for the analytic test suite.

---

## 1. Converge the electronic calculation before the nuclear dynamics

For a small representative geometry set spanning

```text
Franck-Condon region
approach to crossing
minimum-gap region
post-crossing region
```

repeat the electronic calculation while changing one choice at a time.

### SCF tolerance

Example sequence:

```text
1e-8
1e-10
1e-12
```

Compare

- state energies;
- state gradients;
- NAC vectors;
- SCF iteration count.

The default repository backend already fails loudly if SCF does not converge.

---

## 2. CASSCF optimization tolerance

Refine

```text
mc_conv_tol
mc_conv_tol_grad
```

separately.

Do not accept a looser orbital-gradient threshold merely because the state energy has
stopped changing: direct dynamics uses gradients and NACs, not energy alone.

---

## 3. Active-space convergence

The most important qualitative check is orbital/state character.

For every candidate active space, save

```text
active orbital specification
active electron count
state-average weights
tracked energies
tracked state-gradient norms
tracked NAC norms
many-electron overlap matrices between neighboring geometries
```

A state manifold that cannot be tracked continuously through the intended region is
not an adequate direct-dynamics model, even if every individual CASSCF point converges.

---

## 4. Number of state-averaged roots

Repeat with a larger state manifold where feasible.

If the selected two-state subspace is leaking into a third state, the v0.6/v0.7
overlap diagnostics will show reduced singular values or a growing overlap-unitarity
defect.

Converge the **state manifold**, not only the active orbitals.

---

## 5. Geometry-step convergence of state tracking

For a path coordinate $q$, compare tracked scans with

```text
Delta q
Delta q / 2
Delta q / 4
```

Monitor

$$
|\langle\Psi_I(q)|\Psi_I(q+\Delta q)\rangle|,
$$

assignment-score margin, overlap singular values, and the overlap-derived directional
NAC

$$
d^{(q)}\approx\frac{O-O^\dagger}{2\Delta q}.
$$

The latter should approach the projected analytical PySCF NAC as $\Delta q\to0$.

---

## 6. ETF convention

The backend makes `use_etfs` explicit.

Treat calculations with

```text
use_etfs=False
```

and

```text
use_etfs=True
```

as different NAC conventions.  Record the choice in every production manifest.

Do not tune this flag merely to make two curves look similar.

---

## 7. Basis-set convergence

At a representative set of geometries, compare at least two basis levels appropriate
to the target system.

Use tracked states before comparing NACs.  A raw root-by-root comparison can confuse
basis-set effects with root permutation or sign changes.

---

## 8. Centroid electronic-structure convergence

The v0.9 SPA0/SPA1 layer evaluates electronic quantities at Gaussian-pair centroids.

For selected strongly coupled pairs, validate the local approximation by evaluating
additional backend points displaced around the centroid.

Check whether

$$
H_e(q_c+\delta q)
\approx
H_e(q_c)
+
\sum_\alpha
\frac{\partial H_e}{\partial q_\alpha}\delta q_\alpha
$$

is adequate over the region where $|g_i^*g_j|$ is appreciable.

A large v0.9 SPA0/SPA1 difference is a warning that centroid-only electronic
information may be insufficient.

---

## 9. Cache provenance

For a production run, persist

```text
geometry hash
PySCF version
basis
charge/spin
SCF reference
active space
state-average weights
SCF/CASSCF tolerances
ETF choice
tracked root permutation/sign metadata
energies
gradients
NACs
```

alongside every cached electronic point.

Never reuse a cached point under a different electronic-structure contract.

---

## 10. Recommended run order

A reproducible small-system workflow is:

```text
1. validate PySCF backend on isolated geometries
2. converge active space + state manifold
3. converge tracking geometry spacing
4. create tracked electronic-data scan/graph
5. converge exact 2D reference
6. converge Gaussian dt
7. converge spawning action threshold
8. converge max TBF basis size
9. compare SPA0 and first-order Taylor layer
10. report final exact-reference error
```

This order prevents an expensive dynamics calculation from hiding an unconverged
electronic-structure model.
