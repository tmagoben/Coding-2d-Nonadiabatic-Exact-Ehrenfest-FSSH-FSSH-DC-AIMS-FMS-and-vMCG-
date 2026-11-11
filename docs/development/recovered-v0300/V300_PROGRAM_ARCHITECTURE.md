# v0.30.0 Program Architecture

```text
molecular provider snapshot
  H, physical dH, NAC/connection, mass, overlap
                 |
                 v
       strict structural/provenance gate
                 |
        +--------+--------+
        |                 |
        v                 v
 NAC Wilson ODE     physical overlap
                         |
                         v
                 certified polar factor
        |                 |
        +--------+--------+
                 v
        local transport agreement audit
                 |
        +--------+------------------+
        |                           |
        v                           v
 v0.29.3 curved sections     v0.29.4 exact oracle
                                    |
                                    v
                       invariant molecular context
```

## New modules

- `gaussian_dynamics/molecular_frame_bridge_v300.py`
  - fixed-dimension provider and provenance admission;
  - NAC/connection `ConnectionFieldV293` adapter;
  - certified reciprocal overlap-polar transport;
  - local Wilson-versus-overlap transport receipts;
  - electronic-only pointwise Hamiltonian action with an explicit non-full-Hamiltonian name;
  - gauge-safe molecular context records for v0.29.4 oracle data.
- `gaussian_dynamics/v300_benchmark.py`
  - deterministic gauge-covariance, transport, provenance, rejection, and data-contract evidence.
- `tests/test_v0300_molecular_frame_bridge.py`
  - positive and adversarial gates for the new interface.

## Next implementation slice

The next v0.30.x slice should connect this bridge to the full nuclear kinetic plus
electronic Hamiltonian action used by the curved TDVP residual, then execute a short,
replayable molecular trajectory with stepwise overlap/NAC/SOC receipts and dataset
emission. ML training remains later than that trajectory-level certification.
