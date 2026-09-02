# v0.25.3 validation

The release adds 60 deterministic scientific gates and 25 adversarial/core gates.
Together with 630 inherited v0.25.2 gates, cumulative acceptance is 715/715.

Independent checks include:

- dense-grid reconstruction of analytic candidate residual couplings;
- dense-grid reconstruction of projection error;
- exact enlarged-basis wavefunction, norm, fidelity, and energy preservation;
- newborn coefficient growth with exactly frozen dormant shape variables;
- active-subspace extraction from the full v0.25.2 metric;
- no-event endpoint reduction to v0.25.2;
- age/projection-gated pruning and overlap/projection-gated merge-to-survivor;
- packet permutation and constant complex electronic-gauge covariance;
- even singlet/triplet and odd doublet complete-spin models;
- zero-SOC toggle equivalence;
- rank, conditioning, packet cap, activation, identity, and tamper failures.

Canonical evidence:

- `results/v0253_controlled_basis_evidence.json`
- `results/v0253_controlled_basis_campaign.json`
