# v0.23.2 algorithmic complexity

Let R be the number of geometry records, s the retained electronic-state count,
d the nuclear-coordinate count, F the number of installed distribution files,
and B their total bytes.

## Runtime provenance

Reading and validating package metadata is O(F); verifying every RECORD-hashed
file is O(B) time with O(F) manifest storage and bounded streaming hash memory.
This occurs during evidence generation, not every dynamics step.

## Overlap validation

There are R self blocks and R(R-1) directed cross blocks. Identity and
reciprocity checks cost O(R^2 s^2). Singular-value contraction checks dominate at

$$
O(R^2s^3)
$$

time and O(R^2s^2) stored overlap memory. Diagnostics reuse the same decompositions.

## NAC finite-difference certification

For h central-difference step sizes and d selected Cartesian directions, the
certification requires 2hd displaced electronic calculations plus the reference
calculation. Framework-side phase alignment and overlap differentiation cost
O(h d s^3), normally negligible beside SA-CASSCF.

## Admission

Structure, capability, and convergence checks are linear in declared fields.
Artifact hashing is linear in raw bytes. Receipt/replay/dossier binding is linear
in their manifests, while the inherited physical matrix audit retains its dense
state-space and finite-difference costs. Method-specific electronic calculations
remain the dominant production cost.
