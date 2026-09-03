# v0.18 Build Validation Report

Validated on 2026-08-13.

## Source validation

```text
292 Python files parsed successfully with Python AST.
```

## Cumulative automated regression suite

```text
212 passed in 9.69 s
```

The suite is cumulative from v0.1 through v0.18.

New v0.18 coverage includes:

```text
global-phase-invariant full-wavefunction metrics
nuclear-density and moment diagnostics
batched candidate ranking equivalence with v0.15 dense ranking
deterministic sampled omitted-edge audits
sampled-audit violation detection
sampled-audit runner behavior
initial/final dense sentinel behavior
physical-time adaptive cadence resolution
Richardson self-convergence helpers
basis/sparsity/growth acceptance logic
```

## Canonical convergence campaign

Machine-readable output:

```text
results/v018_convergence_complete_campaign.json
```

Independent coordinate results:

```text
results/v018_partials/
```

All 12 convergence coordinates were executed as fresh-process worker jobs.

## Release acceptance

```text
passed = True
```

Every configured acceptance check passes.

## Full-wavefunction result

Projected exact reference at `t=0.6`:

```text
fidelity:
0.982566093411826

phase-aligned L2:
0.13232747836123407

nuclear-density L2:
0.052341235444456596

nuclear-density total variation:
0.050223497471400924

centroid error:
0.001269391116081437

covariance error:
0.01269760151753039
```

Reduced electronic density:

```text
projected-reference Frobenius error:
0.00010573932284646514

original-target Frobenius error:
0.03329249794783041

target population error:
0.028073109470748484

target coherence phase error:
0.0028907634670896944
```

## Representation-error separation

```text
initial projection fidelity:
0.8822514544600691

exact projected-target final fidelity:
0.8822514544600707

maximum exact overlap drift:
1.6653345369377348e-15
```

The exact projected-target overlap is conserved to approximately machine precision.

## Timestep self-convergence

```text
||Psi_0.010 - Psi_0.005||:
0.0006101685310847837

||Psi_0.005 - Psi_0.0025||:
0.00015281576549629837

observed order:
1.9974143869640382
```

The observed order is effectively second order.

Adaptive event cadence is normalized in physical time.

## Basis convergence

```text
Nmax=10:
0.19119002020634157

Nmax=11:
0.17750196410133656

Nmax=12:
0.14492430172716658

Nmax=13:
0.13232747836123407
```

Relative 10 -> 13 improvement:

```text
30.79 %
```

The ladder is strictly improving.

## Sparse-edge-budget convergence

```text
B_local=0.030:
0.14573374595371563

B_local=0.010:
0.13232747836123407

B_local=0:
0.13212429667373482
```

The error is nonincreasing as sparse truncation is relaxed.

## Adaptive-growth sensitivity

```text
threshold 0.050:
N=10
L2=0.19119002020634157

threshold 0.035:
N=12
L2=0.18008929313945907

threshold 0.030:
N=13
L2=0.13820391719791275

threshold 0.025:
N=13
L2=0.13232747836123407

threshold 0.015:
N=13
L2=0.13232747836123407
```

The final two thresholds select the same growth history and form a numerical plateau.

## Sparse audit architecture

```text
normal sampled audits:
6

sampled pairs scored:
6

sampled audit failures:
0

full dense sentinels:
2

v0.18 dense sentinel pair factorizations:
146

v0.17 dense audit pair factorizations:
506

dense audit pair-work reduction:
71.15 %
```

The dense sentinel errors remain below the v0.18 matrix-error limits.

## Candidate-grid batching

```text
maximum unbatched candidate-grid elements:
1044800

batched peak elements:
25600

peak candidate-grid reduction:
97.55 %
```

The batched and dense residual-ranking paths are regression-tested to return the same
candidate ordering and capture fractions.

## Conditioning and conservation

```text
maximum condition number:
6509.218903498147

release ceiling:
10000

runtime candidate condition limit:
100000

generalized norm drift:
1.2434515149761793e-06
```

The larger 13-Gaussian basis is retained because it materially improves the complete
wavefunction while remaining numerically stable.

## Representative examples executed

```text
examples/78_v018_full_wavefunction.py
examples/79_v018_trajectory_fidelity.py
examples/80_v018_basis_convergence.py
examples/81_v018_dt_convergence.py
examples/82_v018_edge_budget.py
examples/83_v018_growth_trigger.py
examples/84_v018_audit_and_memory.py
```

The fresh-process worker architecture used to generate the release partials also
validates the same code path exposed by:

```text
examples/85_v018_coordinate_worker.py
```

## PySCF and SOC scope

The inherited PySCF, many-electron overlap, state tracking, and gauge-graph
infrastructure remains in the cumulative repository.

v0.18 does not claim:

```text
production molecular residual-adaptive dynamics
a calibrated molecular sampled-audit model
spin-orbit-coupled dynamics
production AIMS
```

See `V18_PYSCF_CONVERGENCE_BRIDGE.md`.
