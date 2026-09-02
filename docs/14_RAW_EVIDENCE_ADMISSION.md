# Raw-evidence molecular SOC admission

v0.23.1 closes the gap between an evidence summary and the calculations that produced
it.

## Why v0.23.0 was not enough for a real source

v0.23.0 correctly required reference, basis, method, frame, and tracking evidence, but
the contract stored their final errors and thresholds. Those values were fingerprinted,
yet the framework could not recompute them from raw observations. A synthetically
relabelled dataset could therefore satisfy the summary contract without proving that a
method-specific electronic-structure parser had examined genuine outputs.

v0.23.1 retains the v0.23.0 gates and adds three layers:

1. raw calculation inputs/outputs with SHA-256 and size records;
2. derived evidence and exact replay/receipt/provenance binding;
3. executable backend validation for external or live admission.

## What the framework derives

- independent-reference error from computed and reference arrays;
- every adjacent basis and method change;
- translation and proper-rotation residuals;
- minimum singular-value retention within complete physical manifolds;
- spectral leakage and assignment margin against competing manifolds.

The tracking record graph must connect every replay geometry, and the manifold groups
must partition every state. This handles complete Kramers doublets and triplet
subspaces without depending on arbitrary component orientation inside a degenerate
manifold.

## What remains backend specific

Only a method-specific parser can establish that raw output bytes genuinely encode the
declared SOC method, derivatives, overlaps, and convergence state. Real admission calls
that parser through an executable validator whose name and version must match the
runtime attestation. Hashes and a manually authored attestation cannot bypass it.

## Current status

Both analytic parity-sector fixtures pass the raw-evidence protocol and all 123 release
gates pass. Neither fixture is admitted as real. PySCF is absent, so there is no live
runtime attestation or executable PySCF artifact validation.

The next milestone should supply one pinned, method-specific implementation plus a
real reference molecule and raw outputs. Only then can an external snapshot or live
backend claim become true.
