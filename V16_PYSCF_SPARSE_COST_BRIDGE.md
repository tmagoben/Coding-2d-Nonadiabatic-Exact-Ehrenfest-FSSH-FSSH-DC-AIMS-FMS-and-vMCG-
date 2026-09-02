# v0.16 PySCF Sparse/Cost Bridge

v0.16 introduces two pieces that are directly relevant to molecular direct dynamics:

1. local Gaussian graph sparsity;
2. explicit electronic-structure cache cost in the basis-growth utility.

The release still validates them on the analytic LVC benchmark.

## 1. Local Gaussian graph for molecular TBFs

The locality screen depends only on nuclear Gaussian centers and widths.

Therefore it can be evaluated before a new electronic-structure calculation.

A molecular workflow can use:

```text
cheap nuclear locality screen
        ↓
candidate local degree
        ↓
predicted Gaussian matrix cost
        ↓
electronic geometry-cache estimate
        ↓
residual shortlist
        ↓
only then request expensive new ab-initio data
```

## 2. Electronic cache cost

The current `GeometryCacheElectronicCostModel` is intentionally simple.

It assigns low cost to a candidate near a registered geometry and high cost to a new
geometry.

A real PySCF implementation should replace generic units with observed provider costs:

```text
RHF/ROHF time
CASSCF time
gradient time
NAC time
cross-geometry overlap time
state/gauge alignment time
```

The cost should be attached to the actual provider cache and geometry fingerprint.

## 3. Gauge consistency

Electronic cache reuse is valid only if the electronic frame is also valid.

The inherited v0.6-v0.8 infrastructure already tracks:

- many-electron overlaps;
- root assignment;
- phase/subspace alignment;
- graph gauge transport.

A molecular cost model must not classify a geometry as a cheap cache hit if the stored
electronic state cannot be consistently transported into the current graph frame.

## 4. Sparse graph and electronic graph are related but not identical

The v0.16 nuclear locality graph answers:

> Which Gaussian basis functions are close enough to retain projected S/H/T blocks?

The v0.7-v0.8 electronic gauge graph answers:

> How are local electronic frames transported consistently between electronic
> structure points?

These graphs should eventually be coupled, but they should not be conflated.

A nuclear pair can be locally negligible even though its electronic points remain
useful for gauge transport, and vice versa.

## 5. Expected molecular cost hierarchy

For the analytic LVC benchmark:

```text
Gaussian pair algebra
sparse solve
grid residual
```

are the relevant costs.

For SA-CASSCF direct dynamics, the likely hierarchy can become:

```text
new electronic-structure point
        >>
Gaussian pair algebra
        >>
cached electronic lookup
```

This is exactly why the v0.16 utility contains an explicit provider-cost term.

## 6. What is not implemented

v0.16 does not yet provide:

- measured PySCF cost calibration;
- persistent electronic-structure task scheduling;
- asynchronous electronic evaluations;
- sparse gauge-graph pruning;
- molecular TDSE-defect evaluation;
- automatic sparse-threshold calibration from ab-initio data.

Those are natural later extensions once the analytic sparse architecture is stable.
