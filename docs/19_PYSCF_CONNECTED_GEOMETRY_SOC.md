# PySCF connected-geometry SOC differential preview

v0.24.2 connects the static PySCF BP-SOMF matrices introduced in v0.24.1 across a
controlled OH bond-coordinate stencil. Its purpose is to validate scalable integral
contraction, exact many-electron overlap transport, component separation, provenance,
and finite-difference behavior before those data can enter nuclear dynamics.

## Reproduce the evidence

Use the pinned versions and hashes in
`requirements-pyscf-v242-linux-x86_64-py312.txt`, fix every listed BLAS/OpenMP thread
variable to one, and run:

```bash
python examples/129_recompute_v0242_pyscf_differential_soc.py
python examples/130_recompute_v0242_campaign.py
```

The first command writes the 60-gate runtime evidence. The second writes the full
400-gate campaign and rewrites the same differential evidence from that campaign.

## Reading a record

Each of the three derivative records contains:

- the step and `H_z` coordinate label;
- SHA-256 identities for the common center and exact minus/plus snapshots;
- raw complete-multiplet contractions from center to both endpoints;
- certified unitary right-to-left polar transports;
- transported endpoint `H_spin_free` and `H_soc` matrices;
- independent `K_spin_free`, `K_soc`, and `K_total` centered differences;
- a polar-gauge connection preview and the separate Hermitian contraction slope;
- Hermiticity, component-decomposition, unitarity, and time-reversal residuals.

The scan additionally stores all six compact endpoint receipts, so their geometry,
root data, calculation input, runtime, wavefunction hash, and SOC matrices can be
checked against the derivative record fingerprints.

## Correct interpretation

The raw overlap is not used as an operator transport. Its unitary polar factor is.
This is especially important because the first two equal-weight OH roots are
degenerate and rotate substantially between independent calculations. Root phases
cannot resolve a multidimensional degenerate gauge.

The resulting `K_soc` is a transported finite-difference preview along one OH bond
coordinate. It is not a full molecular SOC gradient. The recorded polar-gauge
connection is not claimed to be a continuous physical derivative coupling. No v0.24.2
object satisfies the trajectory-ready molecular-SOC admission contract.

For the complete derivation, see `V242_PYSCF_DIFFERENTIAL_SOC.md`; for data flow and
claim boundaries, see `V242_PROGRAM_ARCHITECTURE.md`.
