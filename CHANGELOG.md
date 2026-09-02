# Changelog

## 0.27.0 — 2026-08-25

- Replaced coordinate-diagonal packet widths and chirps with full real symmetric
  matrices, using `Gamma=exp(E)` to preserve positive definiteness structurally.
- Added Frobenius-orthonormal `svec` coordinates and the exact matrix-exponential
  Frechet derivative required by full-width McLachlan tangent vectors.
- Added multivariate complex-normal Wick moments through degree four and analytic
  correlated overlap, kinetic, quadratic-potential, metric, and residual elements.
- Extended the full-SVD real McLachlan solve and fully implicit midpoint integrator
  to every matrix degree of freedom with nonlinear, rank, null-space, norm, energy,
  width-domain, and exact dormant-shape receipts.
- Validated arbitrary proper rotations and reflections of the model, wavefunction,
  matrices, metric, forcing, velocity, and implicit midpoint update.
- Added an independent dense FFT quadrature oracle and an exact matrix Riccati oracle;
  the latter demonstrates second-order convergence and correlation generation that
  cannot be represented by the diagonal-width v0.26.0 manifold.
- Generalized spawn, merge, prune, and projection to correlated packets. Spawning
  uses signed intrinsic principal axes and fails closed for degenerate eigenspaces.
- Preserved exact one-dimensional reduction to v0.26.0 and retained all 825 inherited
  gates; added 135 new gates for 960/960 cumulative acceptance.
- Kept coordinate-dependent electronic frames, live molecular-SOC trajectories,
  degenerate-eigenspace direction optimization, and full AIMS branching closed.

## 0.26.0 — 2026-08-24

- Generalized Gaussian centres, momenta, logarithmic widths, and chirps from one
  nuclear coordinate to vector coordinates while retaining exact reduction to the
  v0.25.2 one-dimensional kernel.
- Added exact multidimensional complex Gaussian moments through degree four,
  positive-definite full constant mass matrices, and a fully implicit-midpoint
  McLachlan TDVP solve for coordinate-diagonal adaptive widths.
- Added an implementation-independent two-dimensional FFT/Strang reference solver
  with norm, energy, boundary, reversibility, and second-order refinement receipts.
- Added analytic two-state CI+SOC, complete four-state Kramers-doublet CI+SOC, and
  two-singlet/complete-triplet five-state benchmark Hamiltonians.
- Generalized residual-driven spawning, full-SVD projection, pruning,
  merge-to-survivor, stable packet identity, and one-event checkpoints to multiple
  nuclear coordinates.
- Strengthened newborn activation: coefficient evolution is immediate, while shape
  directions require population, retained-condition, and velocity-amplification
  gates; additional spawning waits while a newborn remains dormant.
- Validated independent dense-grid matrix elements, exact-grid/TDVP comparisons,
  adaptive error improvement, packet/gauge/signed-coordinate covariance, complete
  spin manifolds, Kramers symmetry, zero SOC, and fail-closed receipt tampering.
- Kept full correlated width matrices, arbitrary rotational covariance,
  coordinate-dependent electronic frames, full AIMS branching, real PySCF SOC
  trajectories, and general ab-initio accuracy closed.
- Added 110 new gates for 825/825 cumulative acceptance.

## 0.25.3 — 2026-08-24

- Added analytic residual coupling of orthogonalized Gaussian candidates to
  `dPsi/dt+i*H*Psi`, with explicit novelty, overlap-rank, and conditioning gates.
- Added deterministic position/momentum candidate generation and stable packet IDs,
  ages, monotone serials, and one-event-per-checkpoint lifecycle receipts.
- Added full-SVD fixed-time projection for enlarged and reduced bases, binding
  projection loss, normalized fidelity, linear residual, norm, and energy jump.
- Added coefficient-only newborn activation: electronic amplitudes evolve while an
  initially zero-amplitude packet's center, momentum, log-width, and chirp remain
  frozen until the population gate is crossed.
- Added projection-guarded low-population pruning and high-overlap
  merge-to-survivor events in the frozen order merge, prune, spawn.
- Validated even/odd complete-spin SOC, constant electronic gauge and packet
  permutation covariance, zero-SOC equivalence, exact no-event reduction to
  v0.25.2, and fail-closed receipt tampering.
- Kept general/multidimensional AIMS branching, coordinate-dependent electronic
  frames, real PySCF SOC trajectories, and general accuracy closed.
- Added 85 new gates for 715/715 cumulative acceptance.

## 0.25.2 — 2026-08-24

- Promoted every one-dimensional packet width into a variational logarithmic
  coordinate and added its conjugate real quadratic chirp.
- Added exact complex unequal-width/chirp Gaussian moments through degree four for
  overlap, Hamiltonian, tangent metric, and McLachlan forcing.
- Preserved positive widths structurally with `eta=log(alpha)` and added explicit
  minimum/maximum width, chirp, and per-step logarithmic-width gates.
- Extended fully implicit midpoint to coefficient, center, momentum, log-width, and
  chirp variables with recomputable SVD/nonlinear receipts.
- Validated exact single-packet harmonic thawed-Gaussian equations, a closed-form
  breathing oracle, and coherent-state reduction to the v0.25.1 frozen solver.
- Validated even/odd analytic SOC, adaptive width/chirp evolution, signed reversal,
  packet permutation, constant complex gauge covariance, duplicate-packet null
  spaces, zero SOC, and second-order refinement.
- Kept spawning/pruning, multidimensional/full width matrices, coordinate-dependent
  electronic frames, real PySCF SOC trajectories, and general accuracy closed.
- Added 95 new gates for 630/630 cumulative acceptance.

## 0.25.1 — 2026-08-24

- Added a genuinely coupled, frozen-width multi-Gaussian/complete-spinor ansatz in
  one nuclear coordinate and a fixed electronic frame.
- Added exact analytic overlap, kinetic, quadratic-potential, tangent-metric, and
  McLachlan right-hand-side matrix elements through the required Gaussian moments.
- Added a real-parameter McLachlan metric solve using a full SVD pseudoinverse,
  explicit rank/condition receipts, and fail-closed compatible-null-space auditing.
- Added a fully implicit midpoint residual solve with recomputable nonlinear
  convergence, metric, endpoint, norm, energy, model, and settings receipts.
- Validated complete even singlet/triplet and odd doublet SOC spinors, signed-step
  reversal, Gaussian-permutation covariance, constant complex electronic-gauge
  covariance, zero-SOC equivalence, and second-order timestep convergence.
- Validated the exact continuous reduction to canonical one-packet harmonic
  equations and the compatible rank deficiency of duplicate Gaussian packets.
- Kept adaptive widths, spawning/pruning, coordinate-dependent electronic gauges,
  multidimensional multi-Gaussian motion, real PySCF SOC trajectories, and general
  ab-initio dynamics accuracy closed.
- Added 75 new gates for 535/535 cumulative acceptance.

## 0.25.0 — 2026-08-24

- Added the single-canonical-nuclear-packet / complete-spinor TDVP restriction.
- Added a signed, self-adjoint constant-mass velocity-Verlet nuclear update coupled
  to endpoint-Hamiltonian Strang electronic propagation.
- Computed the physical unitary polar transport from the raw cross-geometry overlap
  by SVD, retaining singular values, conditioning, and principal-angle evidence.
- Added fingerprint-bound step receipts that recompute masses, forces, `H`, physical
  `K`, overlaps, polar factors, energies, endpoints, and time.
- Validated complete even singlet/triplet and odd Kramers-doublet SOC dynamics,
  coordinate-dependent complex-gauge covariance, signed-step reversal, zero-SOC
  equivalence, norm preservation, and second-order convergence.
- Rejected full multi-Gaussian TDVP, adaptive widths, coordinate-dependent-mass
  Verlet, static PySCF SOC snapshots, bad retained manifolds, and receipt tampering.
- Kept real PySCF molecular-SOC trajectories and general ab-initio dynamics accuracy
  false; the recommended future full-TDVP integrator is implicit midpoint/discrete
  variational, not plain Verlet.
- Added 60 new gates for 460/460 cumulative acceptance.

## 0.24.2 — 2026-08-24

- Replaced production rank-five `int2e_p1vxp1` materialization with PySCF direct-JK
  SOMF contractions; retained the explicit route only as a small-system oracle.
- Added overlap-capable, fingerprinted SOC snapshots at one center and six displaced
  OH geometries.
- Added exact restricted-CASSCF root overlaps lifted into complete doublet
  microstates and certified degenerate-safe unitary polar transport.
- Added separate transported `H_spin_free` and `H_soc` centered differences over the
  `0.08/0.04/0.02` bohr second-order ladder.
- Bound every component derivative to four stored transported endpoint matrices and
  every record to serialized minus/plus snapshot receipts.
- Added contraction, retention, time-reversal, Kramers, Hermiticity, Richardson,
  endpoint-provenance, and adversarial tamper controls.
- Explicitly kept the continuous physical connection, full Cartesian and analytic
  derivatives, real mixed multiplicity, trajectory admission, and accuracy false.
- Added 85 new gates for 400/400 cumulative acceptance.

## 0.24.1 — 2026-08-24

- Added direct fixed-geometry PySCF 2.13.1 BP-SOMF state-interaction molecular SOC.
- Added explicit one-/two-electron integral identities, one `1/(2c^2)` prefactor,
  state-averaged SOMF density, and independent PySCF JK contraction cross-check.
- Added project-native Wigner--Eckart/Clebsch--Gordan reconstruction with integer
  twice-quantum numbers and PySCF determinant spin ladders for zero reference CGs.
- Added complete singlet/triplet, doublet/quartet, and Kramers-doublet controls.
- Added a real OH three-doublet SA-CASSCF(5e,4o)/STO-3G evidence calculation that
  returns the direct complex six-state SOC matrix.
- Kept physical SOC derivatives, derivative connections, cross-geometry overlaps,
  trajectory readiness, general ab-initio accuracy, and OpenMolcas admission false.
- Added 59 new gates for 315/315 cumulative acceptance.

## 0.24.0 — 2026-08-21

- Frozen OpenMolcas 26.06 H2O RASSCF/CASPT2/RASSI-SO/AMFI external-SOC
  protocol with exact method, state-order, unit, and convention identity.
- Strict 55-record raw-artifact schema binding native input, output, HDF5,
  exported matrices, cross-geometry overlaps, convergence, environment, and
  independent validation through SHA-256 inventories.
- Separate transported finite-difference certification of spin-free and SOC
  derivatives over all Cartesian coordinates and three step sizes.
- Independent reference, basis, method, rigid-frame, and connected-manifold
  evidence gates.
- Exact parser-type and caller-owned protocol/convention/manifest/environment/
  exporter trust anchors, with relabel, corruption, unknown-file, and subclass
  negative controls.
- Admission-bound frozen-snapshot dynamics, zero-SOC mode, and deterministic restart.
- A 256-gate campaign containing all 208 v0.23.3 gates and 48 new gates.
- The deterministic protocol fixture is not OpenMolcas output; external/live
  molecular SOC and ab-initio SOC claims remain false.
- Native OpenMolcas HDF5/text numerical cross-parsing remains deliberately closed;
  the external admission gate cannot open in this release.

## 0.23.3 — 2026-08-21

- Added certified unitary polar transport while preserving raw finite-manifold
  overlaps as physical contractions.
- Added independent overlap quality gates for retained singular value, condition
  number, and principal angle.
- Added deterministic replay format 2 and explicit NAC convention fingerprints.
- Added evidence-bound legacy replay migration; unknown and wrong-sign NAC data
  are quarantined without automatic sign correction.
- Added convention-complete provider identities for replay/cache/checkpoint use.
- Added complete singlet/triplet and Kramers-doublet transport auditing under
  arbitrary independent endpoint gauges.
- Froze SOC matrix prefactor, operator treatment, state order, units, and physical
  derivative semantics.
- Split runtime claims into canonical `release_locked` and broader
  `scientifically_compatible` profiles.
- Added 40 new gates for 208/208 total acceptance. External/live molecular SOC
  and ab-initio SOC remain unadmitted.

## 0.23.2 — 2026-08-21

- Hash-locked and content-verified PySCF 2.13.1 runtime on CPython 3.12/Linux x86-64.
- Real H3+ SA-CASSCF(2e,3o)/STO-3G energies, analytic gradients, NACs, and
  many-electron overlaps across two geometries.
- Corrected production NAC mapping: internal `d[i,j]` uses PySCF `state=(i,j)`,
  certified by phase-aligned overlap finite differences through 1e-4 bohr.
- Explicit separation of full-overlap (`use_etfs=False`) and ETF-corrected NACs.
- Finite-manifold overlap contract based on identity self overlaps, adjoint
  reciprocity, and contractive cross-overlap singular values rather than isometry.
- Runtime admission hardened with structure-first callable checks, exact trusted
  method identity, seven convergence stages, and typed parser/execution proof.
- A 168-gate campaign containing all 123 v0.23.1 gates and 45 new gates, including
  28 real-runtime checks.
- No external or live molecular-SOC backend admitted; no ab-initio or PySCF SOC claim.

## 0.23.1 — 2026-08-21

- Canonical raw-evidence dossier with SHA-256 and size verification for templates,
  environments, references, calculation inputs/outputs, and runtime probes.
- Seventeen per-sector calculation receipts covering nine trajectory records, three
  basis levels, two method levels, and three rigid-frame calculations.
- Reference errors, adjacent convergence changes, and frame residuals derived from
  stored observations instead of trusted summary values.
- Connected physical-manifold tracking based on minimum singular values, competing
  spectral leakage, and assignment margins.
- Exact dossier-to-replay, receipt-to-coordinate, and raw-summary-to-provenance binding.
- Separate external-snapshot and live-backend outcomes.
- Mandatory executable backend-specific raw-artifact validation in addition to a
  fingerprinted runtime attestation.
- Pinned optional PySCF 2.13.1 target and frozen SA-CASSCF NAC orientation.
- Synthetic relabel, corruption, traversal, duplicate-output, convergence, coordinate,
  environment, evidence, and disconnected-tracking negative controls.
- A 123-gate campaign containing all 93 v0.23.0 gates and 30 new gates.
- No external or live molecular SOC backend admitted; no ab-initio or PySCF SOC claim.

## 0.23.0 — 2026-08-21

- Explicit `static_soc` and derived `trajectory_ready` molecular-SOC capability tiers;
  static-only providers fail closed for moving nuclei.
- Fingerprinted backend, method, unit, electron, nuclear, geometry, calculation-input,
  environment, coordinate, and tracking identity.
- Independent real-admission gates for reference agreement, basis and method
  convergence, translational/rotational invariance, and state-tracking quality.
- Deterministic manifest/NPZ replay with component operators, connections, masses, all
  pair overlaps, convergence flags, time reversal, projectors, and full provenance.
- Exact-coordinate-only replay plus array, manifest, overlap, and cross-dataset
  corruption controls.
- Separate protocol-valid and real-backend-admitted results, preventing analytic
  fixtures from acquiring an ab-initio claim.
- Even singlet–triplet and odd two-Kramers-doublet deterministic reference replays.
- Fail-closed PySCF bridge requiring a live runtime, a method-specific complete
  provider, and affirmative SCF/correlated/SOC convergence.
- A 93-gate release campaign containing all 67 v0.22.1 gates and 26 new gates.
- No real molecular SOC backend admitted; no ab-initio or live PySCF SOC validation.

## 0.22.1 — 2026-08-20

- Independent transported full-matrix finite-difference audits for spin-free and SOC
  derivatives across every coordinate and three displacement sizes.
- Negative cancellation fixture proving that total-K and scalar-force checks cannot
  hide opposite component errors.
- Dimension-neutral SOC audits with no provider `.config` dependency.
- Explicit SOC symmetry admission for electron parity, charge sector, complete
  multiplets, time-reversal unitarity and square, complete projectors, zero external
  field, and numerical symmetry provenance.
- Exact-grid endpoint retention, operator-contract mass, fixed-frame/constant-mass
  certification, moving-frame rejection, and static split-operator precomputation.
- Gaussian-basis and sparse-threshold physical SOC convergence ladders.
- A 67-gate release campaign containing all 53 v0.22.0 gates and 14 corrective gates.
- No molecular or ab-initio SOC backend admitted.

## 0.22.0 — 2026-08-20

- First physical analytic spin-orbit Hamiltonians and physical SOC derivative
  operators in the framework.
- Separate even-electron singlet/complete-triplet and odd-electron two-Kramers-doublet
  model spaces; electron-number sectors are never mixed.
- Explicit even- and odd-electron time-reversal conventions, including transformed
  antiunitary representations under general complex gauges.
- Zero, constant, and coordinate-dependent analytic SOC limits.
- Gauge-covariant physical multiplet/root projectors.
- Independent one-dimensional vector exact-grid Strang propagator with norm, energy,
  timestep, grid-spacing, box-size, and Gaussian-population validation.
- SOC-active dense, sparse, moving-frame, provenance-mismatch, and corruption restart
  checks.
- A 53-gate release campaign containing all 21 v0.21.4 gates and 32 physical-SOC gates.
- No ab-initio SOC or real PySCF runtime claim; spin-free dynamics remains permanent.

## 0.21.4 — 2026-08-20

- Centered cross-geometry provider audit that independently certifies physical
  Hamiltonian derivatives, derivative connections, overlap isometry, structural
  invariants, and provenance identity.
- Negative fixtures proving that pointwise-valid but geometry-inconsistent K or D data
  are rejected.
- Versioned deterministic checkpoint/restart for self-consistent dense and block-sparse
  dynamics, with SHA-256 integrity, provider/settings fingerprints, density-guide
  history, sparse-edge hysteresis state, and stable Gaussian UIDs.
- Global-step continuity for adaptive insertion/pruning after restart.
- Explicit zero-SOC decomposition rehearsal through the future H/K composition path,
  with exact spin-free operator and trajectory equivalence.
- A 21-gate campaign inheriting all v0.21.3 acceptance checks.
- No physical SOC Hamiltonian, no physical SOC derivatives, and no ab-initio SOC claim.

## 0.21.3 — 2026-08-13

- Frobenius-residual structural validation with no implicit relative tolerance.
- Fixed-dimension electronic model-space and complete-multiplet declarations.
- Fingerprinted operator provenance, explicit hartree/bohr internal units, and frozen
  physical Hamiltonian-derivative semantics.
- Gauge-covariant density-matrix nuclear guidance at exact degeneracy and vanishing
  local electronic amplitude.
- Explicit retirement of the gauge-dependent lowest-eigenvector fallback.
- Arbitrary-state, arbitrary-nuclear-dimension fixed-frame initial projection.
- Fixed-frame complex-operator cache keyed by complete provider provenance.
- A 20-gate SOC-contract-freeze campaign inheriting all v0.21.2 acceptance checks.
- No physical SOC Hamiltonian and no ab-initio SOC claim.

## 0.21.2 — 2026-08-13

- Unequal-width complex block Gaussian algebra.
- Self-consistent representation-neutral block nuclear guidance.
- Adaptive zero-block insertion and metric-projected pruning.
- Generic electronic observable interface.
- Full-subspace continuity provider diagnostics.
- Complex-dtype pre-SOC audit.
- Explicit package discovery for reproducible clean editable installs.
- No SOC Hamiltonian added.

## 0.21.0 — 2026-08-13

### Added

- General complex `ElectronicOperatorPointV21`.
- Physical Hamiltonian-derivative operator matrices.
- Complex non-Abelian derivative connection contract.
- Smooth analytic complex $U(s)$ gauge-transform provider.
- Arbitrary-state synthetic operator provider.
- Full-subspace Procrustes alignment.
- Wilson-loop gauge-spectrum utilities.
- Block-sparse molecular Gaussian $S/H/T$ matrices.
- Gauge-invariant block edge score.
- Prescribed moving-basis block propagator.
- Time-dependent complex-gauge equivalence campaign.
- Dynamic sparse topology stress test.
- 2/4/8-state block scaling validation.
- Curated mathematical and algorithmic documentation.

### Changed

- The core representation is no longer defined by real adiabatic roots.
- Sparse block importance no longer uses the norm of the full gauge connection $T$,
  because that norm is representation dependent under coordinate-dependent gauges.
- Documentation now separates implementation claims, validation claims, and future
  physics.

### Not added

- No SOC Hamiltonian.
- No spin-specific propagation assumptions.
