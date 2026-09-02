# v0.8 Build Validation Report

Validated on 2026-08-12.

## Compilation

All Python files in `gaussian_dynamics/`, `examples/`, and `tests/` compiled successfully.

## Automated suite

```text
88 passed
```

This is the cumulative v0.1-v0.8 suite.

New v0.8 coverage includes:

- independent-endpoint $U(2)$ gauge covariance of overlap/local-diabatic propagation;
- convergence of explicit-NAC and overlap propagation to the same analytic CI result;
- temporal-link unitarity;
- incremental center/centroid graph cycles;
- exact discrete metric compatibility,
  $\dot S=T+T^\dagger$;
- preservation of the seed anti-Hermitian moving-basis connection;
- changing-metric norm conservation;
- dynamic graph growth;
- dynamic zero-amplitude spawning;
- nonzero child amplitude from later coupled propagation;
- incremental many-electron snapshot graph construction;
- raw overlap singular-value and unitarity-defect diagnostics.

## Explicit NAC versus overlap/local-diabatic benchmark

```text
Explicit-NAC coefficients:
[ 0.63054048+9.04367112e-03j -0.77610367+6.90860172e-05j]

Overlap/local-diabatic coefficients:
[ 0.63054096+9.04366945e-03j -0.77610329+6.90864174e-05j]

Fidelity:
0.9999999999996232

Norms:
0.9999999999999993
1.0000000000000064
```

The two propagation formulations therefore agree to essentially machine-visible precision on the refined analytic CI path.

## Dynamic graph-AIMS-style benchmark

```text
Spawn events:
[{'step': 1, 'time': 0.0002, 'parent_uid': 0, 'child_uid': 1, 'target_state': 0}]

Final coefficients:
[1.00000001e+00-1.93233209e-04j
 1.26334734e-04+1.81839129e-08j]

Final time:
0.01 au

Final norm:
1.0

Final basis size:
2

Final graph:
151 nodes
200 edges
50 independent cycles

Final overlap-matrix condition number:
1.0005054666351296
```

The child is inserted with zero coefficient and later acquires nonzero amplitude through the coupled moving-basis equations.

## Metric-compatible moving basis

Representative residual:

```text
||Sdot - T - T^dagger||_F
= 1.962615573354719e-17
```

## PySCF runtime status

PySCF is optional and is not installed in the build environment.

The v0.8 PySCF path adds a public raw point+snapshot entry point and an incremental many-electron snapshot graph builder. Its actual PySCF binary execution must be validated on a PySCF-enabled machine with:

```bash
pip install -e ".[pyscf]"
python examples/25_pyscf_incremental_snapshot_graph.py
```

The inherited v0.5-v0.7 fake-backend and many-electron-overlap tests remain passing.

## Interpretation

The tests establish internal numerical, gauge, and graph consistency for the analytic benchmark. They do not establish chemical convergence for a real molecule or make this repository a production AIMS package.
