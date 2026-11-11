# v0.30.0 Development Checkpoint Notes

## Outcome

v0.30.0 safely begins the molecular integration layer. It connects
coordinate-dependent provider NAC/derivative connections to the v0.29.3 Wilson
geometry, independently constructs transport from certified physical overlaps, and
requires local agreement before those data enter the curved Gaussian basis.

The v0.29.4 exact residual-spawning oracle remains unchanged. Molecular information is
attached in a separate, fingerprint-bound context containing gauge-invariant spectra
and norms rather than raw frame arrays.

## Validation summary

- portable scientific suite: **611/611 passed**;
- v0.27 + v0.29.2-v0.29.4 + v0.30.0 focused integration: **58/58 passed**;
- v0.30.0 focused tests: **10/10 passed**;
- deterministic v0.30.0 evidence gates: **12/12 passed**;
- maximum NAC-Wilson/overlap-polar disagreement: `3.6911e-16`;
- maximum transport unitarity residual: `6.6603e-16`.

The only excluded cumulative meta-tests demand the exact Python executable used to
seal v0.27.0. The current environment satisfies the scientific compatibility profile
but correctly does not claim byte-identical release-machine identity.

## Next concrete milestone

Bind the bridge to the complete nuclear-kinetic plus electronic Hamiltonian action in
the curved TDVP residual, then run a short replayable molecular trajectory with
stepwise NAC/SOC/overlap receipts and deterministic `spawn-dataset-v1` emission.

ML remains restricted to later candidate ranking/shortlisting and may not override the
exact admission oracle.
