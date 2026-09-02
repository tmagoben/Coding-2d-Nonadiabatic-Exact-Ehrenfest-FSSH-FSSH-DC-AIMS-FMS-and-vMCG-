# v0.13 Release Notes

Version 0.13 is the **residual-driven Gaussian basis refinement release**.

## New modules

```text
residual_basis_v13.py
tdse_defect_v13.py
v13_benchmark.py
```

## New capabilities

- exact one-candidate Hilbert residual gain;
- deterministic Gaussian candidate dictionaries;
- vectorized/prepared dictionary screening;
- monotonic pure residual-greedy basis construction;
- top-K residual plus initial-density screening;
- overlap-conditioning rejection before basis insertion;
- instantaneous TDSE/Galerkin defect reconstruction;
- defect-capture ranking;
- zero-coefficient defect-driven basis enrichment;
- predicted-versus-actual defect reduction validation;
- v0.12 -> v0.13 representation-consistent benchmark comparison.

## Release reference

```text
basis size: 11
projection fidelity: 0.8902521956060818
relative residual: 0.10974780439391818
initial density error: 0.03209140317550961

projected-state dynamics error:
0.00011354880287339317

original-target density error:
0.03178630139393256

population error:
0.025521902605714804

coherence phase error / rad:
0.0023799927838891125

maximum norm drift:
1.0593904686828637e-06

maximum condition number:
3465.8914579773386
```

## Defect enrichment

```text
candidate:
dq=(0.0, -0.4);dp=(0.0, 0.0);width_scale=4

defect before:
0.31502411763651

defect after:
0.28652498129723825

predicted squared reduction:
0.01714362854505765

actual squared reduction:
0.017143629785278974
```

The new candidate enters with exactly zero electronic coefficient.

Thus basis expansion changes the available Galerkin tangent space without changing the
instantaneous represented wavefunction.

## Important terminology

The release should be described as:

> residual-driven spinor-complete Gaussian basis refinement and TDSE-defect enrichment
> on an analytic LVC benchmark.

It should **not** be described as:

- a production residual-AIMS implementation;
- full vMCG;
- a generic molecular TDSE error estimator;
- a completed PySCF residual-spawning engine.

## Automated validation

The final cumulative suite reports:

```text
158 passed
```

Full details are recorded in `V13_BUILD_VALIDATION.md`.
