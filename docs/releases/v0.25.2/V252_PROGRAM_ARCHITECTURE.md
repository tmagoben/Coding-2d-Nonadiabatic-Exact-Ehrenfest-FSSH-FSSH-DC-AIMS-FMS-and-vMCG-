# v0.25.2 program architecture

```mermaid
flowchart TD
    P[Complete fixed-frame analytic SOC provider] --> A[Inherited quadratic and provenance audit]
    A --> H[H0, H1, H2, mass]
    S[Normalized state: C, q, p, alpha, beta] --> L[eta = log alpha]
    L --> T[Degree 0/1/2 adaptive tangents]
    H --> M[Exact complex moments M0 through M4]
    T --> M
    M --> G[McLachlan G and b]
    G --> V[Full SVD pseudoinverse]
    V --> Q{PSD, rank, condition, null-RHS, residual}
    Q -->|fail| X[Fail closed]
    Q -->|pass| F[Adaptive TDVP vector field]
    F --> I[Implicit midpoint nonlinear residual]
    I --> W{HYBR success, residual, width/chirp domain, log-step gate}
    W -->|fail| X
    W -->|pass| E[Adaptive endpoint]
    E --> R[Bound metric/SVD/nonlinear/norm/energy receipt]
    R --> N[Next signed step]
    R --> D[70 validation + 25 core + 535 inherited = 630 gates]
```

Text fallback:

```text
fixed-frame quadratic H + {C,q,p,alpha,beta}
                |          alpha=exp(eta)>0
                +------------------+
                                   v
                 exact chirped moments M0...M4
                                   |
                         G(theta) theta_dot=b(theta)
                                   |
                  full SVD + compatible-null audit
                                   |
              theta1-theta0-h*v((theta0+theta1)/2)=0
                                   |
             HYBR success + residual + width/chirp gates
                                   |
             endpoint and fingerprint-bound step receipt
```

## Module ownership

- `adaptive_multigaussian_tdvp_v252.py`: adaptive state, exact moments/matrices,
  tangents, metric, implicit step, receipts, trajectory, and frozen claims.
- `adaptive_multigaussian_tdvp_validation_v252.py`: even/odd, reversal, adaptive
  motion, covariance, null-space, harmonic, coherent reduction, zero-SOC, and order
  evidence.
- `v252_benchmark.py`: 25 adversarial gates and inheritance of all 535 v0.25.1
  checks.

Spawning/pruning, multidimensional/full-matrix widths, coordinate-dependent
electronic frames, and real molecular providers have no route into this runner.

