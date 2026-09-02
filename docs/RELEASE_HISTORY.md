# Release history

- **v0.1–v0.4:** Heller Gaussians, moving nonorthogonal bases, exact-grid references, conical-intersection models, Berry/Wilson topology, and spawning prototypes.
- **v0.5–v0.8:** provider-neutral molecular interfaces, PySCF-oriented backends, many-electron overlap tracking, and electronic gauge graphs.
- **v0.9–v0.10:** convergence controls, conditioning/pruning, repeated spawning, and corrected reduced electronic density matrices.
- **v0.12:** representation consistency and separation of initial representation error from projected-state dynamics error.
- **v0.13:** residual-driven basis refinement and TDSE-defect capture.
- **v0.14:** time-adaptive residual control, growth/pruning, and complexity profiling.
- **v0.15–v0.17:** caching, sparse locality, cost-aware adaptation, sparse propagation, and error-controlled audits.
- **v0.18:** full-wavefunction reconstruction and combined basis/timestep/sparsity convergence.
- **v0.19:** molecular direct-dynamics integration, overlap-capable snapshots, indexed state tracking, and PySCF snapshot bridge.
- **v0.20:** end-to-end sparse molecular machinery with active-edge molecular S/H/T and sampled audits.
- **v0.21:** general complex electronic operators, arbitrary U(s) covariance, full electronic blocks per Gaussian, gauge-invariant block sparsity, Procrustes subspace tracking, Wilson-loop validation, arbitrary-state scaling, and dynamic topology. No spin physics is introduced.

## v0.21.2 — pre-SOC integration hardening

Generalized the v0.21 complex/block framework through the remaining integration layers
needed before physical spin terms: unequal Gaussian widths, coefficient-coupled nuclear
guidance, zero-block adaptive insertion, metric-projected block pruning, generic
Hermitian electronic observables, full-subspace continuity diagnostics, and an explicit
complex-dtype audit. The packaged release also uses explicit Python-package discovery
so the numerical results directory remains release data and clean editable installs are
reproducible. No SOC Hamiltonian is introduced.

## v0.21.3 — SOC-contract freeze

Closed the remaining procedural defects before physical SOC: structural matrices now
use explicit scaled Frobenius residuals rather than `allclose` defaults; exactly
degenerate or empty electronic blocks use transported density-matrix guidance rather
than an arbitrary eigenvector; model-space dimension, multiplet completeness, internal
units, physical derivative semantics, and provider provenance are explicit; grid
projection is general in electronic and nuclear dimension; and fixed-frame complex
cache entries are separated by a complete provenance fingerprint. All 20 release gates
pass and inherit v0.21.2 acceptance. No physical SOC Hamiltonian is introduced.

## v0.21.4 — differential-provider and deterministic-restart certification

Added centered cross-geometry certification that distinguishes pointwise structural
validity from actual H/K/D differential consistency, including a coordinate-dependent
complex frame and wrong-K/wrong-D negative controls. Added versioned checkpoint/restart
for Gaussian state, density-guide memory, sparse-edge hysteresis, global adaptive-step
semantics, provider/settings identity, and full-state SHA-256 integrity. The future
operator-composition seam is rehearsed with explicit zero H_SOC and K_SOC and reproduces
spin-free dynamics. All 21 release gates pass and inherit v0.21.3 acceptance. No physical
SOC Hamiltonian or derivative is introduced.

## v0.22.0 — first physical analytic SOC

Introduced nonzero physical spin-orbit Hamiltonians and analytic physical SOC derivative
operators in two separate, exactly reproducible model families: an even-electron
singlet plus complete triplet and an odd-electron pair of complete Kramers doublets.
Added explicit time-reversal representations, quaternionic doublet coupling,
gauge-covariant physical projectors, independent vector exact-grid propagation, and
SOC-active restart/failure controls. All 53 release gates pass, including all 21
v0.21.4 gates. The release does not claim ab-initio or molecular SOC validation.

## v0.22.1 — corrective SOC hardening

Closed pre-backend weaknesses in the analytic SOC milestone. Spin-free and SOC
derivatives are now audited separately as transported full matrices in arbitrary state
and coordinate dimensions; an adversarial fixture proves cancellation cannot evade the
audit. Added explicit electron-parity/charge admission, independent time-reversal
unitarity, numerical symmetry provenance, exact-grid frame/mass/endpoint corrections,
and Gaussian-basis plus sparse-threshold SOC convergence. All 67 gates pass, including
all 53 v0.22.0 gates. No molecular SOC backend is admitted.

## v0.23.0 — molecular SOC admission protocol

Added explicit static and trajectory-ready capability tiers, traceable molecular and
environment identity, five independent evidence families, and deterministic
component-resolved replay with exact-coordinate and integrity enforcement. Protocol
validity is now separate from real-backend admission, so analytic even-electron and
odd-electron fixtures can certify the framework without claiming molecular accuracy.
Added a fail-closed boundary for a future method-specific PySCF provider. All 93 gates
pass, including all 67 v0.22.1 gates. No real molecular SOC backend, ab-initio SOC
accuracy, or live PySCF SOC runtime is validated.

## v0.23.1 — raw-evidence admission hardening

Replaced self-contained evidence summaries with deterministic dossiers containing raw
artifact hashes, per-calculation receipts, and derived reference, convergence, frame,
and connected-manifold tracking metrics. Bound every trajectory receipt to an exact
replay coordinate and required an executable method-specific artifact validator in
addition to runtime attestation for external or live admission. A synthetic relabel
control demonstrates that v0.23.0 summary evidence alone is now insufficient. All 123
gates pass, including all 93 v0.23.0 gates. No external or live molecular SOC backend,
ab-initio SOC accuracy, or live PySCF runtime is validated.

## v0.23.2 — real PySCF spin-free runtime and contract correction

Installed and content-verified the pinned PySCF 2.13.1 runtime and executed a real
H3+ SA-CASSCF spin-free fixture for energies, analytic gradients, NACs, and
many-electron overlaps. Corrected the production NAC tuple mapping to
`state=(i,j)` using phase-aligned overlap finite differences, separated full
overlap and ETF semantics, and replaced exact cross-geometry isometry with the
physical finite-manifold contraction contract. Strengthened future SOC admission
with callable structure checks, exact trusted method identity, seven convergence
stages, and typed parser/execution proof. All 168 gates pass, including all 123
v0.23.1 gates. No external/live molecular-SOC backend or ab-initio SOC claim is
admitted.

## v0.24.0 — fail-closed OpenMolcas RASSI-SO snapshot intake

Frozen the first method-specific external-SOC protocol and added a strict 55-record
artifact parser, SHA-256 trust chain, transported Cartesian component derivatives,
independent reference/basis/method/frame/tracking evidence, exact caller trust policy,
and admission-bound frozen-snapshot dynamics. All 256 gates pass, including all 208
v0.23.3 gates. The campaign uses a conspicuously synthetic protocol fixture; no
OpenMolcas execution, external/live molecular-SOC source, or ab-initio SOC result is
admitted.

## v0.25.0 — symmetric restricted-TDVP SOC dynamics

Connected the complete representation-neutral SOC operator contract to a symmetric
single-canonical-packet trajectory. Constant-mass nuclear variables use velocity
Verlet, while the complete electronic spinor uses endpoint Strang propagation and
the unitary polar factor computed by SVD from the cross-geometry overlap. Even and
odd analytic SOC manifolds, complex-gauge covariance, signed-step reversal,
zero-SOC equivalence, norm preservation, and second-order convergence are validated.
All 460 gates pass, including all 400 v0.24.2 gates. Full multi-Gaussian/adaptive-
width TDVP, general noncanonical Verlet, real PySCF SOC trajectories, and general
ab-initio dynamics accuracy remain false.

## v0.25.1 — frozen-width multi-Gaussian TDVP metric layer

Replaced the single-packet restriction with a coupled one-dimensional sum of
frozen-width Gaussian packets carrying complete spinors. Exact analytic moments
build the real McLachlan metric and forcing; a full SVD provides rank-aware,
compatible-null minimum-norm velocities; and fully implicit midpoint advances every
coefficient, center, and momentum together. Even/odd analytic SOC, signed reversal,
packet permutation, constant complex gauge, duplicate-packet null spaces, harmonic
reduction, zero SOC, and second-order refinement are validated. All 535 gates pass,
including all 460 v0.25.0 gates. Adaptive widths, spawning/pruning,
coordinate-dependent frames, multidimensional multi-Gaussian motion, real PySCF SOC
trajectories, and general ab-initio accuracy remain false.

## v0.25.2 — adaptive log-width/quadratic-chirp TDVP

Promoted every one-dimensional packet width into `eta=log(alpha)` and paired it
with a real quadratic chirp. Exact complex moments through degree four now build the
adaptive McLachlan metric; full SVD and implicit midpoint propagate coefficients,
centers, momenta, widths, and chirps together. Closed-form harmonic breathing,
coherent reduction to v0.25.1, even/odd SOC, signed reversal, packet permutation,
constant complex gauge covariance, duplicate-packet null spaces, zero SOC, and
second-order refinement are validated. All 630 gates pass, including all 535
v0.25.1 gates. Spawning/pruning, multidimensional/full width matrices, moving
electronic frames, real molecular SOC trajectories, and general accuracy remain
false.

## v0.25.3 — controlled adaptive-basis lifecycle

Wrapped the unchanged v0.25.2 adaptive-width TDVP kernel in a conservative
one-dimensional topology controller. Candidate Gaussians are scored by their exact
analytic coupling to the McLachlan residual after orthogonalization against the
current nuclear basis. Novelty, full-rank overlap, condition number, SVD projection
loss, fidelity, norm, and energy-jump gates bind every accepted spawn. Newborn
electronic coefficients activate immediately while their initially undefined shape
coordinates remain frozen until population is sufficient for a stable full metric.
Age-eligible low-population pruning and high-overlap merge-to-survivor events use the
same projection receipts. All 715 gates pass, including all 630 v0.25.2 gates.
General/multidimensional AIMS, moving electronic frames, real molecular-SOC
trajectories, and general accuracy remain false.

## v0.26.0 — reference-first multidimensional CI+SOC dynamics

Generalized centers, momenta, coordinate-diagonal log widths, chirps, and positive-
definite mass matrices to multiple nuclear dimensions. Added an implementation-
independent 2D FFT/Strang CI+SOC oracle, complete doublet and singlet/triplet analytic
models, multidimensional residual spawning, and metric-gated newborn activation. All
825 gates pass, including all 715 v0.25.3 gates. Full correlated width matrices,
arbitrary rotations, moving electronic frames, and live molecular SOC remain false.

## v0.27.0 — full correlated-width rotation-covariant TDVP

Promoted every packet to full symmetric width and chirp matrices with structural SPD
widths through `Gamma=exp(E)` and Frobenius-orthonormal `svec` coordinates. Exact
multivariate moments, matrix-exponential Frechet tangents, full-SVD McLachlan dynamics,
implicit midpoint, and the controlled lifecycle now operate on the complete matrix
manifold. Dense FFT quadrature, exact matrix Riccati dynamics, proper rotations,
reflections, gauge/permutation covariance, and intrinsic-axis spawning are validated.
All 960 gates pass, including all 825 v0.26.0 gates. Moving electronic frames, live
molecular SOC trajectories, degenerate-direction optimization, and full AIMS remain
false.

## v0.24.2 — connected-geometry direct-JK PySCF SOC differentials

Moved the production two-electron BP-SOMF contraction onto PySCF's direct JK driver,
avoiding stored rank-five AO integrals while retaining the explicit small-system route
as a validation oracle. Added seven independently converged OH SA-CASSCF snapshots,
exact many-electron overlaps, complete-doublet polar transport, separate transported
spin-free/SOC centered differences, a frozen second-order step ladder, and six
serialized endpoint receipts. All 400 gates pass, including all 315 v0.24.1 gates.
The result is a connected-geometry differential preview only: continuous physical
connections, full Cartesian/analytic derivatives, real mixed multiplicity, general
accuracy, and trajectory-ready admission remain false.

## v0.24.1 — direct static PySCF BP-SOMF molecular SOC

Added an independent state-interaction SOC implementation on the pinned PySCF 2.13.1
runtime. Converged common-orbital SA-CASSCF roots now produce direct complex
`H_spin_free`, `H_soc`, and `H_total` matrices in complete `|root,S,M_S>` multiplets
using all-electron Breit-Pauli one-electron integrals, a two-electron spin-orbit
mean-field contraction, spin-resolved transition densities, and project-native spin
algebra. A real OH three-doublet calculation validates Hermiticity, time reversal,
Kramers pairing, a nonzero SOC signal, and agreement of the explicit SOMF contraction
with an independent PySCF JK path. All 315 gates pass, including all 256 v0.24.0 gates.
The provider is intentionally static-only: physical SOC derivatives, cross-geometry
tracking, trajectory-ready admission, general accuracy, and external OpenMolcas
admission remain false.

## v0.23.3 — finite-manifold transport and compatibility hardening

Separated raw finite-manifold contractions from unitary polar coefficient
transport and added independent physical-consistency and trajectory-readiness
criteria. Introduced deterministic replay format 2, fingerprinted NAC conventions,
attested migration, and fail-closed quarantine for ambiguous legacy signs. Added
complete singlet/triplet and Kramers-doublet endpoint-gauge transport audits,
frozen molecular-SOC matrix semantics, convention-complete provider identity, and
separate release-locked/scientifically-compatible runtime profiles. All 208 gates
pass, including all 168 v0.23.2 gates. No external/live molecular-SOC backend or
ab-initio SOC claim is admitted.
