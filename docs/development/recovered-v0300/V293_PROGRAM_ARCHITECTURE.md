# v0.29.3 Program Architecture

```text
coordinate-dependent A_mu(R)
        |
        v
explicit path gamma(q -> R)
        |
        v
path-ordered Wilson transport W_gamma
        |
        +------> curvature / plaquette holonomy audit
        |
        +------> gauge transform + patch transition audit
        v
curved packet section g_I(R) W_gamma(R<-q_I) c_I
        |
        v
complete finite-difference real TDVP tangent space
        |
        v
v0.29.2 Schur-complement residual spawning
        |
        v
full-SVD projection + merge/prune/spawn lifecycle
```

## New modules

- `gaussian_dynamics/curved_geometry_v293.py`
  - coordinate-dependent anti-Hermitian connection providers;
  - exact/ODE Wilson transport;
  - straight and polyline/dogleg paths;
  - curvature and plaquette holonomy;
  - coordinate-dependent gauge transforms;
  - patch transitions and same-path cross-patch transport.

- `gaussian_dynamics/curved_basis_adaptation_v293.py`
  - coordinate-dependent curved coefficient basis;
  - full real finite-difference tangent Jacobian;
  - TDVP snapshots;
  - residual-driven candidate evaluation;
  - full-SVD projection and spawning;
  - curved merge/prune/spawn lifecycle under an explicit path family.

- `gaussian_dynamics/v293_benchmark.py`
  - deterministic scientific evidence and fingerprint.

- `tests/test_v0293_coordinate_dependent_geometry.py`
  - geometry, covariance, holonomy, patch, reduction and lifecycle gates.

## Inherited exact boundaries

- constant connection + straight path reduces to v0.29.2 section arrays;
- v0.29.2 residual-gain equation is unchanged;
- v0.27 correlated Gaussian state and lifecycle thresholds remain inherited;
- no second spawn while a dormant newborn exists;
- one topology event per adaptation checkpoint remains the default.
