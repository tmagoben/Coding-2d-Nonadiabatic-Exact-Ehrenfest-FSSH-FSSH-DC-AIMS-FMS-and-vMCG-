# Molecular backends

The repository separates electronic calculation, cross-geometry state information, and Gaussian propagation.

## Electronic calculation

A backend supplies electronic data at one molecular geometry. The intended real route remains SA-CASSCF/PySCF, while analytic and synthetic providers are used for deterministic validation.

## v0.21 operator conversion

A spin-free adiabatic point with energies \(E_i\), gradients, and derivative couplings \(D_a\) is lifted to the general complex operator contract.

In the adiabatic representation,

$$
H_{ij}=E_i\delta_{ij},
$$

$$
(K_a)_{ii}=\partial_aE_i,
$$

and for \(i\neq j\),

$$
\boxed{(K_a)_{ij}=(E_j-E_i)(D_a)_{ij}.}
$$

The block propagator itself does not require \(H\) to be diagonal.

## PySCF derivative-coupling convention

The internal convention is

$$
d_{ij}=\langle\Phi_i|\nabla_R\Phi_j\rangle.
$$

PySCF documentation describes its tuple as `(ket, bra)`. Real PySCF 2.13.1
overlap finite differences establish that the production mapping required for
internal `d[i,j]` is `state=(i,j)` with `mult_ediff=False` and `use_etfs=False`.
The earlier literal `state=(j,i)` interpretation had the opposite sign. The
`mult_ediff=True` and ETF-corrected quantities remain distinct diagnostics rather
than the full-overlap dynamics NAC. See the v0.23.2 erratum.

## State identity and warm starts

Warm-start orbitals can reduce SCF/CASSCF iteration cost. They do not determine root identity, phase/gauge continuity, or degenerate-subspace orientation. Those require cross-geometry overlap tracking.

## Degenerate manifolds

At exact degeneracy, individual adiabatic eigenvectors are non-unique. v0.21 therefore exposes full-subspace Procrustes alignment instead of forcing root-by-root identity.

## v0.23.0 molecular SOC admission boundary

Molecular SOC providers now declare a `static_soc` or derived `trajectory_ready`
capability tier. Moving-nuclear dynamics requires component derivatives, derivative
connections, and cross-geometry overlaps in addition to the static Hamiltonian.

A real source must also bind molecule/isotope/geometry identity, calculation and
environment hashes, and independent reference, basis, method, frame-invariance, and
tracking evidence. Deterministic replay transports these data without changing their
source classification. See `13_MOLECULAR_SOC_ADMISSION.md` and the root-level
`V230_MOLECULAR_SOC_CONTRACT.md`.

## v0.23.1 raw-evidence boundary

Evidence summaries are now derived from fingerprinted raw observations and linked to
per-calculation receipts. Degenerate multiplets are tracked by subspace singular values
and competing-manifold leakage on a connected record graph. External or live admission
must execute a method-specific raw-artifact validator; hashes and a stored attestation
alone are insufficient. See `14_RAW_EVIDENCE_ADMISSION.md`.

## v0.23.2 runtime boundary

PySCF 2.13.1 is installed and content-verified in the release runtime. A real
spin-free SA-CASSCF fixture validates energies, analytic gradients, NACs, and
many-electron overlaps. This does not establish a method-specific SOC Hamiltonian,
physical SOC derivatives, complete spin manifolds, SOC raw-output semantics, or
molecular-SOC accuracy. No external or live molecular-SOC backend is admitted.
