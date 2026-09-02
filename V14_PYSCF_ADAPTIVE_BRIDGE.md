# v0.14 PySCF Adaptive-Control Bridge

v0.14's complete TDSE defect is available because the analytic LVC benchmark can be
reconstructed on a compact nuclear grid.

A full molecular PySCF direct-dynamics calculation cannot evaluate the global
many-dimensional TDSE residual in this way.

The adaptive-control architecture is nevertheless transferable.

## 1. Reusable control layer

The following concepts do not depend on the analytic LVC Hamiltonian:

```text
separate add/remove thresholds
adaptation cooldown
basis-size budget
minimum independent candidate norm
conditioning budget
exact low-loss pruning of the represented Gaussian state
event provenance
complexity ledger
```

These can be reused directly.

## 2. Residual surrogate required for molecules

The analytic quantity

$$
\mathcal R=i\dot\Psi-H\Psi
$$

must be replaced by a controlled molecular surrogate.

Possible choices include:

- sampled residuals over local Gaussian quadrature points;
- SPA0 versus SPA1 matrix-element discrepancy;
- local-diabatization overlap residuals;
- projected electronic-subspace transport mismatch;
- coefficient-equation residuals on a local electronic gauge graph.

The surrogate must first be calibrated against analytic models where the true TDSE
defect is known.

## 3. Candidate generation

The v0.14 separation

```text
physics-based candidate generation
        +
residual-based candidate ranking
```

is particularly useful for ab initio dynamics.

Expensive PySCF evaluations should be reserved for a small shortlist rather than every
candidate.

A molecular workflow could therefore use:

```text
cheap geometric candidate dictionary
        ↓
cached/interpolated residual pre-screen
        ↓
top-K candidates
        ↓
new PySCF electronic-structure points
        ↓
many-electron overlap/gauge update
        ↓
final acceptance
```

## 4. Complexity difference

For the analytic benchmark, Gaussian algebra dominates the measured runtime.

For PySCF SA-CASSCF, the electronic-structure evaluation may dominate by orders of
magnitude.

Therefore a molecular complexity ledger should add categories such as:

```text
SCF seconds
CASSCF seconds
gradient seconds
NAC seconds
many-electron overlap seconds
state/gauge alignment seconds
cache hit/miss counts
```

The Gaussian dynamics ledger in v0.14 is intentionally separate from those future
provider costs.

## 5. State/gauge consistency

Any molecular residual comparing electronic amplitudes at different geometries must
first put them into a consistent tracked subspace.

The v0.6-v0.8 machinery already provides:

- many-electron overlaps;
- state assignment;
- Procrustes/polar subspace alignment;
- dynamic gauge graph transport.

Residual adaptation must sit **after** that alignment layer.

Otherwise the residual would contain arbitrary electronic phase/root-ordering error.

## 6. What is not claimed

v0.14 does not claim:

- full-dimensional molecular TDSE-defect evaluation;
- online PySCF residual spawning;
- production AIMS basis management;
- automatic active-space refinement;
- linear-scaling molecular dynamics.

It provides the tested numerical controller that those future pieces would plug into.
