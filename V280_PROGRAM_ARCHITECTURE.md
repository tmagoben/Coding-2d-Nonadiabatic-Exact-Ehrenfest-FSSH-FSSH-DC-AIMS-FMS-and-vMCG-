# v0.28.0 Program Architecture

```mermaid
flowchart TD
    F["Flat analytic moving frame G(R)"] --> D["Connection D = G† dG"]
    F --> W["Packet transporter W(R,q) = G†(R)G(q)"]
    S["Center-frame packet coefficients c_I"] --> W
    W --> P["Parallel-transported electronic section"]
    P --> T["Exact trivialization c_ref = G(q)c"]
    T --> V27["Sealed v0.27 correlated TDVP + lifecycle"]
    V27 --> MV["Moving-frame velocity / midpoint / event"]
    F --> L["Independent gauge-link lattice oracle"]
    L --> E["Hermiticity + similarity + action + propagation checks"]
    MV --> EVID["50-gate v0.28 development evidence"]
    E --> EVID
```

New modules are `moving_frame_v280.py`, `moving_frame_validation_v280.py`, and
`moving_frame_evidence_v280.py`. The v0.27 mathematical modules remain unchanged and
serve as the fixed-frame physical reference for the admitted flat gauge class.
