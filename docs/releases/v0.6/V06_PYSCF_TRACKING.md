# PySCF State Tracking in v0.6

The v0.6 implementation is split into four layers.

```text
PySCF SA-CASSCF calculation
        |
        v
CASSCFWavefunctionSnapshot
        |
        v
cross-geometry many-electron overlap O_ij
        |
        v
maximum-overlap assignment + gauge
        |
        v
tracked E, gradients, NACs, CI roots
```

## Primary classes/modules

```text
gaussian_dynamics/pyscf_tracked_backend_v06.py
gaussian_dynamics/pyscf_wavefunction_overlap.py
gaussian_dynamics/state_tracking.py
gaussian_dynamics/overlap_transport.py
gaussian_dynamics/tracked_scan.py
```

---

## 1. PySCF APIs used

The overlap engine uses PySCF's cross-AO overlap integral:

```python
gto.intor_cross(
    "int1e_ovlp_sph",
    mol_previous,
    mol_current,
)
```

and its FCI overlap helper:

```python
fci.addons.overlap(
    ci_previous,
    ci_current,
    norb,
    nelec,
    s=mo_cross_overlap,
)
```

The SA-CASSCF object supplies:

```python
mc.mo_coeff
mc.ci
mc.ncore
mc.ncas
mc.nelecas
mc.e_states
```

PySCF documents `mc.nelecas` as the active
`(nalpha, nbeta)` electron tuple and `ncore` as the number of core orbitals for
restricted CASSCF.

---

## 2. Why `mc.ci` is retained

At every accepted geometry, v0.6 stores a snapshot containing

```text
molecule
MO coefficients
CI vector for every root
ncore
ncas
nelecas
```

The current raw CI roots are overlapped with the **tracked and phase-corrected**
previous roots.

Thus the selected electronic gauge is propagated from one geometry to the next.

---

## 3. First geometry

The first geometry has no predecessor.

Its raw PySCF root order defines the initial labels:

```text
tracked state 0 = raw root 0
tracked state 1 = raw root 1
...
```

and all initial state phases are set to `+1`.

Every later point is defined relative to that propagated history.

---

## 4. Root crossing behavior

If the raw current overlap matrix is, for example,

```text
[[ 0.05, -0.96],
 [ 0.94,  0.03]]
```

then the best assignment is

```text
tracked 0 <- raw 1
tracked 1 <- raw 0
```

and the first assigned current CI vector is multiplied by `-1`.

As a result, tracked energies may no longer be sorted numerically:

```text
raw energies:     [E_low, E_high]
tracked energies: [E_character_A, E_character_B]
```

That is expected.

The purpose of tracking is persistent state identity, not sorted output.

---

## 5. Ambiguity policy

The backend accepts

```python
ambiguity_policy="raise"
ambiguity_policy="warn"
ambiguity_policy="accept"
```

The default is `raise`.

Two thresholds are exposed:

```python
minimum_overlap
minimum_score_margin
```

If either continuity is poor or the best root assignment is nearly tied with another
assignment, the backend does not silently decide that the states are well defined.

---

## 6. Degenerate subspaces

The state tracker is intentionally conservative at exact degeneracy.

Use

```python
subspace_overlap_singular_values(...)
current_to_previous_procrustes(...)
```

to analyze the entire degenerate block.

The correct object at exact degeneracy is the subspace, not an arbitrary individual
eigenvector.

---

## 7. Sequential use

Typical geometry scan:

```python
backend = PySCFTrackedSACASSCFBackend(config)

for geometry in geometries_in_path_order:
    point = backend.evaluate(geometry)
```

Typical generalized-coordinate path:

```python
cart_backend = PySCFTrackedSACASSCFBackend(config)

provider = GeneralizedCoordinateProvider(
    cart_backend,
    geometry_map,
)

scan = run_tracked_scan(
    q_path,
    provider,
)
```

Afterward:

```python
interp = TrackedScan1DProvider(scan)
```

can be queried in arbitrary order.

---

## 8. Do not share one sequential tracker among unrelated centroid calls

This is intentionally called out because it would create order-dependent science.

Do **not** do:

```python
tracked_backend.evaluate(pair_centroid_A)
tracked_backend.evaluate(unrelated_TBF_B)
tracked_backend.evaluate(pair_centroid_C)
```

and assume all labels now belong to one physically meaningful path.

For the current v0.6 framework:

- sequential TBF center -> one tracking history;
- 1D Gaussian pair data -> use a pretracked scan/interpolator;
- general multidimensional pair graph -> future graph-gauge release.

---

## 9. PySCF runtime status in the build environment

The build environment used to create this package does not contain the PySCF binary.

Therefore the release validates the PySCF state-tracking layer through:

1. official PySCF API/source alignment;
2. a deterministic fake-PySCF call-contract test inherited from v0.5;
3. a deterministic nonorthogonal FCI-overlap test engine;
4. backend-independent root crossing/sign/NAC tests.

For a real installation, run the supplied PySCF examples after installing:

```bash
pip install -e ".[pyscf]"
```
