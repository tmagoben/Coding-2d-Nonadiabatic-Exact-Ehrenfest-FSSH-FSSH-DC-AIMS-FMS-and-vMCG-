# v0.25.3 program architecture

```mermaid
flowchart TD
    A[Normalized adaptive state<br/>stable IDs and ages] --> B[v0.25.2 implicit TDVP step]
    B --> C{Dormant packet?}
    C -- yes --> D[Reduced coefficient-activation metric<br/>freeze dormant q,p,eta,beta]
    C -- no --> E[Full v0.25.2 metric]
    D --> F[Completed fixed-basis endpoint]
    E --> F
    F --> G{Adaptation checkpoint?}
    G -- no --> N[No-op receipt]
    G -- yes --> H{Admissible merge?}
    H -- yes --> P[Full-SVD projection<br/>merge-to-survivor]
    H -- no --> I{Admissible prune?}
    I -- yes --> Q[Full-SVD projection<br/>remove one packet]
    I -- no --> J[Generate 4 candidates per parent]
    J --> K[Analytic residual coupling]
    K --> L[Novelty, rank, condition,<br/>score and packet-cap gates]
    L --> M{Candidate admitted?}
    M -- no --> N
    M -- yes --> R[Full-SVD enlarged-basis projection]
    P --> S[Lifecycle event receipt]
    Q --> S
    R --> S
    N --> S
    S --> T[Next state, IDs, ages, serial]
```

Text fallback:

```text
state -> fixed-basis TDVP/activation -> checkpoint
      -> merge gate -> prune gate -> residual-spawn gate -> no-op
      -> one SVD-bound event -> updated state + stable metadata
```

Primary modules:

- `adaptive_multigaussian_tdvp_v252.py`: unchanged fixed-basis variational kernel.
- `controlled_basis_adaptation_v253.py`: candidates, projection, activation,
  lifecycle events, and variable-basis trajectory.
- `controlled_basis_validation_v253.py`: 60 deterministic scientific gates.
- `v253_benchmark.py`: 715-gate cumulative campaign.
