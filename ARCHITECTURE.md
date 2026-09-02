# Architecture

```text
Electronic backend
    |
    |  H, operator derivatives, derivative connection,
    |  cross-geometry overlap representation
    v
Electronic operator provider
    |
    |  optional complex gauge / subspace alignment / caching
    v
Gaussian nuclear basis
    |
    +--> geometric locality graph
    |       |
    |       +--> candidate pair centroids
    |       +--> block S/H/T score
    |       +--> omitted-score budget
    |
    v
Sparse block S/H/T matrices
    |
    +--> metric-compatible moving-basis connection
    |
    v
Sparse midpoint/Cayley coefficient propagation
    |
    +--> norm / conditioning / gauge / dense-reference diagnostics
    |
    v
Validation and convergence campaign
```

The electronic backend and Gaussian propagation engine are intentionally separated.

A future SOC backend may change the contents of $H$ without changing the Gaussian
algebra.

## v0.24.0 external-SOC evidence plane

```text
OpenMolcas raw bundle + independent validation
                  |
                  v
strict digest/parser boundary <--- caller-owned trust policy
                  |
          +-------+-------+
          |               |
          v               v
 transported H_sf/H_SOC   accuracy / convergence /
 derivatives + manifolds  frame / tracking audits
          |               |
          +-------+-------+
                  v
       external snapshot admission
                  |
                  v
 admission-bound electronic dynamics
```

Protocol fixtures terminate at diagnostics. They cannot enter the admitted dynamics
path. See `V240_PROGRAM_ARCHITECTURE.md` for the complete diagram.
