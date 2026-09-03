# v0.23.3 molecular-SOC convention freeze

Before a real molecular-SOC source can be admitted, its matrix meaning must match
a caller-trusted convention exactly. v0.23.3 fingerprints:

- the SOC operator family and one-/two-electron treatment;
- mean-field approximation and numerical prefactor convention;
- scalar-relativistic treatment;
- source and target basis representations;
- exact ordered state labels and electron parity;
- Cartesian component order and spin quantization axis;
- hartree energy and bohr coordinate units;
- physical fixed-frame SOC derivative semantics;
- complete-multiplet and zero-external-field requirements.

Matrix and derivative Hermiticity, state count/order, symmetry projectors,
electron parity, units, method label, and scalar-relativistic choice are audited
against provider provenance. A hidden prefactor change or reordered state basis
crosses the trust boundary and is rejected even if dimensions still match.

The release exercises this contract on complete analytic singlet/triplet and
two-Kramers-doublet fixtures. That validates the framework convention, not an
ab-initio molecular-SOC implementation.
