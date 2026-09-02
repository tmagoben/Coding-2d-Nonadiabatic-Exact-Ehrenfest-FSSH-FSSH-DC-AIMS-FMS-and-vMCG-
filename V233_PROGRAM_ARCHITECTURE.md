# v0.23.3 program architecture

```mermaid
flowchart LR
    Q[Geometry / configuration] --> P{Electronic provider}
    P --> A[Analytic SOC fixtures]
    P --> Y[PySCF 2.13.1 spin-free]
    P --> R[Replay format 2]
    Y --> N[NAC convention identity]
    R --> N
    N --> O[Raw finite-manifold overlap]
    O --> C{Physical contraction?}
    C -->|yes| U[Unitary polar transport]
    U --> M[Complete manifold audit]
    M --> I[Convention-complete provider identity]
    A --> S[SOC matrix convention audit]
    R --> S
    S --> I
    I --> E[Electronic operator point]
    E --> G[Gaussian S / H / T assembly]
    G --> D[Dense or block-sparse dynamics]
    D --> K[Checkpoint / observables / results]

    L[Legacy replay] --> T{Explicit NAC attestation}
    T -->|known corrected| R
    T -->|unknown or wrong sign| X[Quarantine]

    Y -. no turnkey molecular SOC .-> B[Future method-specific SOC backend]
    B -. v0.23.3 closed .-> Z[External/live admission]
```

## Contract flow

The raw overlap is physical evidence; the unitary polar factor is the coefficient
transport. Complete manifold/time-reversal checks and exact NAC/SOC convention
fingerprints are composed into the provider numerical identity. Replays, caches,
and checkpoints may only cross this boundary when the full identity matches.

The runtime profile is orthogonal to the physics contracts: release-locked means
canonical execution identity, while scientifically compatible means supported
versions without a byte-identity claim.
