# v0.15 PySCF Cost-Aware Adaptation Bridge

v0.15 optimizes the analytic Gaussian layer.

For a molecular PySCF calculation, the electronic-structure provider can be far more
expensive than the Gaussian pair algebra. The cost-aware controller therefore needs a
broader cost vector before it becomes useful for ab initio production work.

## 1. Reusable v0.15 pieces

The following ideas transfer directly:

```text
candidate residual shortlist
condition-aware utility
cost horizon
add/remove hysteresis
adaptation cooldown
incremental Gaussian matrix updates
cache hit/miss accounting
separate physical and computational acceptance criteria
```

## 2. Additional molecular costs

A PySCF-backed utility should include measured or predicted costs for:

```text
RHF/ROHF
CASSCF macro/micro iterations
state-averaged roots
analytic gradients
NAC evaluation
cross-geometry AO overlaps
many-electron wavefunction overlaps
state assignment / gauge alignment
cache misses for new geometries
```

These terms can dominate the Gaussian dynamics cost.

## 3. Geometry-level electronic cache

The existing repository already has geometry fingerprints and electronic-structure
caching.

A molecular cost-aware candidate could therefore distinguish:

```text
candidate near an already cached electronic geometry
candidate requiring a new expensive SA-CASSCF point
```

even if both add one Gaussian to the quantum basis.

This is where cost-aware ranking becomes more physically consequential than in the
analytic LVC benchmark, where all candidates have similar structural compute cost.

## 4. State/gauge requirement

Any cached electronic data must remain attached to the correct tracked electronic
subspace.

The v0.6-v0.8 many-electron overlap and gauge-graph machinery must therefore be part of
the cache key/validation path.

Reusing a cheap but gauge-inconsistent electronic point would be an invalid
optimization.

## 5. Suggested future molecular utility

One possible structure is

$$
U_c
=
\frac{
\text{predicted residual reduction}
}{
\text{Gaussian dynamics cost}
+
\text{electronic-structure cost}
+
\text{conditioning penalty}
}.
$$

The denominator should be calibrated from the actual running backend, not hard-coded
from generic complexity theory.

## 6. What v0.15 does not claim

v0.15 does not yet implement:

- online PySCF wall-time prediction;
- electronic-structure-aware candidate ranking;
- persistent molecular pair caches across arbitrary asynchronous TBF graphs;
- sparse/local coefficient solves;
- production molecular residual-AIMS.

The release establishes the audited analytic architecture those features can build on.
