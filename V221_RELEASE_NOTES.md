# v0.22.1 release notes

v0.22.1 is a corrective hardening release for the analytic spin-orbit coupling
framework introduced in v0.22.0. It does not add a molecular SOC backend. It closes
contract, derivative-audit, exact-grid, and convergence gaps that should not be carried
into the planned v0.23 backend integration.

## Corrected

- Spin-free and SOC derivatives are now checked independently as transported full
  matrices for every nuclear coordinate and every finite-difference step. Opposite
  component errors can no longer cancel in the total derivative.
- SOC admission now requires one electron-parity sector, one charge sector, complete
  multiplets, a unitary time-reversal representation with the declared square, complete
  orthogonal physical projectors, zero external magnetic field, and numerical symmetry
  identity in provider provenance.
- Physical-SOC audits derive state and coordinate dimensions from emitted data and no
  longer depend on a provider-specific `.config` object.
- The one-dimensional exact-grid oracle reads mass from the electronic operator
  contract, certifies a constant scalar mass and fixed electronic frame, rejects moving
  frames, precomputes the static split operators once, and always stores the true final
  step even when `steps` is not divisible by `store_every`.
- Gauge and projector utilities now fail closed on nonfinite, nonunitary, non-Hermitian,
  or non-idempotent inputs.
- Campaign check values are canonical native Python booleans.

## Added validation

- An arbitrary three-state, two-coordinate, provider-config-free SOC fixture.
- A negative fixture whose spin-free and SOC derivative errors cancel in the total K
  and in the sampled scalar force.
- Mixed electron-parity, nonunitary time-reversal, and numerical-provenance mismatch
  controls.
- A prescribed 1/3/5-Gaussian SOC population convergence ladder.
- A four-level sparse-threshold SOC convergence ladder against the dense trajectory.
- Exact-grid endpoint, mass-contract, fixed-frame, and precomputation equivalence tests.

## Acceptance result

The canonical campaign passes **67/67 gates**: all 53 v0.22.0 gates plus 14 corrective
v0.22.1 gates. Exact measured results are recorded in
`results/v0221_corrective_hardening_campaign.json` and summarized in
`V221_VALIDATION.md`.

## Compatibility and scope

- The physical analytic singlet–triplet and two-Kramers-doublet models remain.
- Spin-free dynamics remains a permanent first-class mode.
- The v0.22.0 archive is preserved unchanged.
- No ab-initio SOC, molecular SOC accuracy, real PySCF SOC runtime, external magnetic
  field, or production AIMS claim is made.

