# v0.29.4 to ML Roadmap

1. v0.29.4: freeze oracle/data contracts and leakage policy.
2. v0.30.x: connect live molecular electronic frames/NAC/SOC and emit
   `spawn-dataset-v1` records during trajectories. **Started in v0.30.0:** the
   provider/frame bridge, independent NAC-versus-overlap transport audit, and
   gauge-safe molecular context contract are implemented. Full trajectory binding
   and trajectory-time emission remain open.
3. Expand group splits from trajectory IDs to molecule/model/trajectory families and hold out entire physical regimes.
4. Establish exact-oracle cost baselines and candidate-count scaling.
5. Train the first continuous surrogate for \(\Delta_c\), initially a small interpretable model.
6. Use ML only for shortlist/ranking; run the exact oracle on shortlisted candidates.
7. Exact-audit a deterministic/random subset of ML rejections to estimate false-negative rates.
8. Add uncertainty/OOD gating before any broader control authority is considered.
