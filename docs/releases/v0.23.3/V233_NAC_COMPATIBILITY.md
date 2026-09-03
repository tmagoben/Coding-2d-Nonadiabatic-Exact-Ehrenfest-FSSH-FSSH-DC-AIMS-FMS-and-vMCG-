# v0.23.3 derivative-coupling compatibility

The internal convention is frozen as

$$
d_{ij}=\langle\Phi_i|\nabla_R\Phi_j\rangle
$$

in inverse bohr. Its numerical identity includes the source backend and version,
source mapping, ETF policy, energy-difference scaling policy, and coordinate unit.

For PySCF 2.13.1, production full-overlap dynamics uses `state=(i,j)`,
`use_etfs=False`, and `mult_ediff=False`. This is the corrected mapping certified
in v0.23.2 by phase-aligned many-electron-overlap finite differences. ETF-enabled
couplings have a distinct diagnostic identity and cannot be substituted silently.

Every v0.23.3 replay snapshot exposes the convention fingerprint. Missing or
mismatched identities are rejected. Provider numerical identities bind this NAC
fingerprint together with provenance, overlap/transport contracts, overlap policy,
and replay version, so caches and checkpoints cannot cross convention boundaries.
