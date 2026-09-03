# Gaussian Nonadiabatic Dynamics

**Current release: v0.27.0 — rotation-covariant full-correlated-width Gaussian TDVP.**

v0.27.0 promotes every multidimensional packet to a full real symmetric positive-
definite width matrix and a full real symmetric chirp matrix. Width positivity is
structural through `Gamma=exp(E)`, and the orthonormal `svec` parameterization makes
the real McLachlan metric and fully implicit midpoint update covariant under arbitrary
orthogonal coordinate transformations, including reflections.

Canonical v0.27.0 acceptance: **960/960 checks passing**—all 825 v0.26.0 gates
plus 100 correlated-width scientific-validation gates and 35 adversarial/core controls.

See [`V270_RELEASE_NOTES.md`](V270_RELEASE_NOTES.md),
[`V270_CORRELATED_WIDTH_TDVP.md`](V270_CORRELATED_WIDTH_TDVP.md),
[`V270_ROTATION_COVARIANCE.md`](V270_ROTATION_COVARIANCE.md), and
[`V270_PROGRAM_ARCHITECTURE.md`](V270_PROGRAM_ARCHITECTURE.md).


## v0.28.0 development branch

The `develop` branch is extending the sealed v0.27.0 correlated-width manifold to
**coordinate-dependent electronic frames**. The admitted first milestone uses exact
parallel-transported electronic sections

```text
Psi(R) = sum_I g_I(R) Phi(R) W(R,q_I) c_I,
W(R,q) = G(R)^dagger G(q),
D_a(R) = G(R)^dagger partial_a G(R),
```

for analytically trivializable flat pure-gauge connections. Packet coefficients live
in the electronic frame at their own centers. The current development evidence passes
50/50 frozen moving-frame gates and an independent gauge-link lattice oracle.

This development milestone does **not** admit nonzero-curvature connections, live
molecular-SOC trajectories, general ab-initio SOC-dynamics accuracy, or full AIMS
branching semantics. See `V280_MOVING_FRAME.md` and `V280_VALIDATION.md`.

---

A research-oriented Python framework for **Gaussian-basis nonadiabatic quantum dynamics**
with explicit treatment of nonorthogonal moving bases, electronic gauge freedom,
adaptive/sparse Gaussian coupling graphs, molecular electronic-structure providers,
and exact-model validation.

The current release is **v0.27.0**.

> **v0.27.0 validates full correlated width/chirp matrices, exact multivariate
> complex-normal moments, arbitrary orthogonal coordinate covariance, and controlled
> intrinsic-principal-axis basis adaptation on fixed-frame quadratic spin models.**
> It is not full AIMS, a coordinate-dependent electronic-frame formulation, or a
> trajectory-ready molecular-SOC backend. Those remain closed gates.

## What v0.21 adds

The central electronic data contract is now representation neutral:

$$
H(q)\in\mathbb{C}^{s\times s},
\qquad
\partial_a H_{\mathrm{op}}(q)\in\mathbb{C}^{s\times s},
\qquad
D_a(q)=\langle\Phi|\partial_a\Phi\rangle\in\mathbb{C}^{s\times s}.
$$

The framework requires

$$
H=H^\dagger,
\qquad
\partial_a H_{\mathrm{op}} =
(\partial_a H_{\mathrm{op}})^\dagger,
\qquad
D_a=-D_a^\dagger.
$$

A local electronic basis may be transformed by any smooth

$$
G(q)\in U(s),
$$

with

$$
H' = G^\dagger H G,
$$

$$
(\partial_a H_{\mathrm{op}})' =
G^\dagger(\partial_a H_{\mathrm{op}})G,
$$

and

$$
D'_a =
G^\dagger D_a G + G^\dagger\partial_a G.
$$

v0.27.0 inherits numerical validation of these identities and propagates full
**electronic blocks** on each Gaussian rather than hard-coding a two-state or
one-active-state structure.

## What v0.27.0 adds

The nuclear packet is now

$$
g_I(\mathbf R)=
\left(\frac{\det\Gamma_I}{\pi^D}\right)^{1/4}
\exp\!\left[-\frac12\mathbf x_I^T(\Gamma_I-iB_I)\mathbf x_I
+i\mathbf p_I^T\mathbf x_I\right],
\qquad \mathbf x_I=\mathbf R-\mathbf q_I,
$$

where $\Gamma_I=\exp(E_I)$ is positive definite by construction and both $E_I$
and $B_I$ are real symmetric. Symmetric matrices use an orthonormal packing:
diagonal entries are unchanged and off-diagonal entries are multiplied by
$\sqrt 2$. Exact correlated Gaussian moments through total degree four provide
all overlap, Hamiltonian, tangent-metric, and residual-spawn matrix elements.

Independent validation includes dense FFT quadrature of correlated matrix elements,
the exact matrix Riccati equations for a rotated quadratic Hamiltonian, second-order
implicit-midpoint convergence, proper rotations, reflections, constant electronic
gauges, packet permutations, and exact reduction to v0.26.0 in one dimension.

The controlled lifecycle now proposes signed position and momentum displacements
along the intrinsic eigenvectors of each nondegenerate width matrix. Degenerate or
near-degenerate eigenspaces have no unique intrinsic directions, so candidate
generation fails closed instead of silently selecting laboratory axes. Merge,
prune, novelty, conditioning, projection, residual, stable-identity, and dormant
full-shape activation gates remain mandatory.

## What v0.26.0 adds

The released multidimensional packet is

$$
g_I(\mathbf R)=
\prod_{\mu=1}^{D}\left(\frac{\alpha_{I\mu}}{\pi}\right)^{1/4}
\exp\!\left[-\frac{\alpha_{I\mu}-i\beta_{I\mu}}{2}
(R_\mu-q_{I\mu})^2+i p_{I\mu}(R_\mu-q_{I\mu})\right],
$$

with every positive width represented by `eta[I,mu]=log(alpha[I,mu])`.  Exact
complex-normal moments through total operator/tangent degree four build the
McLachlan metric, Hamiltonian, and residual candidate scores analytically.  The
full real parameter vector is advanced by an implicit-midpoint nonlinear solve;
independent velocity-Verlet updates are not used on this noncanonical manifold.

The reference solver independently propagates

$$
i\partial_t\boldsymbol\Psi(\mathbf R,t)=
\left[-\frac12\nabla^T M^{-1}\nabla\,\mathbf 1
+\mathbf V(\mathbf R)\right]\boldsymbol\Psi(\mathbf R,t)
$$

on a periodic two-dimensional FFT grid.  It uses pointwise Hermitian matrix
exponentials and a unitary Strang kinetic step, and shares no Gaussian moment,
metric, projection, or spawning implementation with the variational solver.

The multidimensional lifecycle retains the frozen **merge, prune, spawn** order.
New candidates are displaced along every signed coordinate and momentum axis and
must pass residual, novelty, full-rank overlap, and conditioning gates.  Newborn
coefficients evolve immediately.  Their shape coordinates activate only after both
the population threshold and two additional metric-safety gates pass: retained
condition number below $10^8$ and velocity amplification below $100$.

The diagonal-width manifold is covariant under signed coordinate permutations but
not under arbitrary rotations of anisotropic widths.  General rotational covariance
requires full correlated complex width matrices and is deliberately deferred rather
than approximated silently.

## What v0.25.3 adds

At each configured checkpoint, the controller permits at most one topology event in
the frozen order **merge, prune, spawn**. A trial spawn packet is admitted only when
its component orthogonal to the existing nuclear basis has sufficient novelty, the
enlarged overlap remains full rank and conditioned, and its analytic coupling to

$$
R=\dot\Psi+i\hat H\Psi
$$

exceeds the residual threshold. Every accepted topology change solves the full-SVD
least-squares projection and records projection loss, fidelity, norm, energy jump,
rank, condition number, packet IDs, and packet ages.

An exactly projected newborn has zero coefficient, making its shape coordinates
temporarily undefined. v0.25.3 therefore evolves all newborn electronic coefficients
immediately while freezing only that packet's `q`, `p`, `log(width)`, and chirp until
its coefficient-row population crosses the audited activation gate. This avoids
silently regularizing the null shape metric.

## What v0.25.2 adds

Each packet is now

$$
g_I(x)=\left(\frac{\alpha_I}{\pi}\right)^{1/4}
\exp\left[-\frac{\alpha_I}{2}(x-q_I)^2
+\frac{i\beta_I}{2}(x-q_I)^2+ip_I(x-q_I)\right],
$$

with `eta_I=log(alpha_I)` and real chirp `beta_I` propagated alongside `q_I`,
`p_I`, and the complete electronic coefficients. Log-width coordinates make
`alpha_I>0` structural; the chirp supplies the quadratic phase required for genuine
Gaussian breathing. Tangents now reach polynomial degree two, so the analytic metric
uses exact complex cross moments through degree four.

For a scalar harmonic potential, the continuous equations reduce exactly to

$$
\dot q=\frac{p}{m},\qquad \dot p=-m\omega^2q,
\qquad \dot\eta=-\frac{2\beta}{m},\qquad
\dot\beta=\frac{\alpha^2-\beta^2}{m}-m\omega^2.
$$

The release compares these equations with a closed-form harmonic breathing solution
and proves that `alpha=m omega`, `beta=0` reduces to the v0.25.1 frozen coherent-state
trajectory. Even/odd SOC, width/chirp reversal, packet permutation, constant complex
electronic gauges, compatible duplicate-packet null spaces, zero SOC, and second-order
timestep refinement are also validated.

## What v0.25.1 adds

For frozen normalized nuclear Gaussians `g_I(x)` and complete electronic spinors
`C_I`, the released wavefunction is

$$
\Psi(x,t)=\sum_{I=1}^{N_g}g_I(x;q_I,p_I,\alpha_I)C_I(t),
\qquad \dot\alpha_I=0.
$$

All real degrees of freedom are propagated together. With
`theta=(Re C, Im C, q, p)`, McLachlan variation gives

$$
G_{\mu\nu}\dot\theta_\nu=b_\mu,
\quad G_{\mu\nu}=\operatorname{Re}\langle\partial_\mu\Psi|\partial_\nu\Psi\rangle,
\quad b_\mu=\operatorname{Im}\langle\partial_\mu\Psi|\hat H\Psi\rangle.
$$

The metric is solved by full SVD. Directions below the absolute/relative cutoff are
retained as explicit null-space evidence; the projected right-hand side must be
compatible or the solve fails closed. The timestep solves

$$
R(\theta_{n+1})=\theta_{n+1}-\theta_n-h\,v\!\left(
\frac{\theta_n+\theta_{n+1}}2\right)=0
$$

with a nonlinear implicit-midpoint solve and stores a receipt that independently
rebuilds the midpoint metric, SVD spectrum, velocity, residual, norm, and energy.

The exact released Hamiltonian contract is one-dimensional and fixed-frame,
`H(x)=H0+x H1+x^2 H2`, with positive constant nuclear mass and a complete spin
manifold. Complete even singlet/triplet and odd Kramers-doublet analytic SOC models
are covered. Gaussian permutation, constant complex electronic gauges, exact signed
reversal, zero SOC, compatible duplicate-packet null spaces, canonical harmonic
reduction, and second-order timestep refinement are validated.

## What v0.25.0 adds

The new trajectory is the single-canonical-nuclear-packet / complete-electronic-
spinor restriction of a time-dependent variational principle. Its equations are

$$
\dot q=M^{-1}p,\qquad
\dot p_a=-c^\dagger K_a c,\qquad
i\dot c=(H-i\dot q^aD_a)c.
$$

For constant `M`, one step applies a nuclear half kick, a drift, an endpoint
electronic Strang step, and the final force half kick. Electronic frame motion is
handled without an uncontrolled derivative-coupling estimate: if the raw retained-
manifold overlap is `O=U Sigma V^dagger`, the unitary polar transport is
`W=UV^dagger`. SVD therefore computes the physical polar factor and simultaneously
provides the singular-value, condition-number, and principal-angle quality gates.

Signed-step reversibility, even singlet/triplet and odd Kramers-doublet SOC models,
coordinate-dependent complex gauges, norm preservation, raw-contraction/polar
separation, zero-SOC equivalence, and second-order convergence are validated.

This does not promote velocity Verlet to the full TDVP solver. The future coupled
multi-Gaussian/adaptive-width manifold requires an implicit midpoint or discrete
variational solve with explicit metric, constraint, and gauge-null-space controls.
The real PySCF SOC objects remain static/differential evidence and are not trajectory
admitted.

## What v0.24.2 adds

The production BP-SOMF construction now evaluates the `int2e_p1vxp1` Coulomb- and
exchange-like contractions with PySCF's direct JK driver. Only three AO output
matrices are retained instead of a `(3,n,n,n,n)` tensor. The explicit tensor remains
only as a small central-geometry oracle and agrees with the direct path to machine
precision.

The canonical OH case solves independent ROHF/equal-weight three-root
SA-CASSCF(5e,4o)/STO-3G wavefunctions at the center and at
`q0 +/- {0.08,0.04,0.02}` bohr. Exact restricted-CASSCF root overlaps are lifted into
the six complete doublet microstates. Their certified unitary polar factors transport
`H_spin_free` and `H_soc` separately to the center before differencing. This handles
the large arbitrary rotation observed inside the degenerate first-two-root subspace;
root-by-root phases are insufficient.

All six endpoint receipts are serialized and fingerprint-bound to the three
derivative records. Each record stores its transported endpoint component matrices,
so equal-and-opposite `K_spin_free`/`K_soc` tampering is rejected even if `K_total`
is unchanged. Both component differences exhibit the expected second-order plateau.

The anti-Hermitian polar-aligned overlap slope is explicitly only a local discrete
gauge preview. v0.24.2 does not claim a continuous physical derivative connection,
a full Cartesian/analytic derivative, mixed-multiplicity runtime coverage, general
accuracy, or trajectory admission.

## What v0.24.1 adds

`PySCFStateInteractionSOCProviderV241` consumes converged common-orbital SA-CASSCF
roots and returns direct `H_spin_free`, `H_soc`, and `H_total` matrices, their SOC
eigensystem, exact root/multiplicity/`M_S` order, time-reversal operator, root
projectors, integral identities, units, and content-addressed provenance.

The pinned OH radical smoke case uses three SA-CASSCF(5e,4o)/STO-3G doublet roots.
Its six-state SOC matrix is Hermitian, nonzero, time-reversal invariant, and Kramers
paired. The explicit rank-five two-electron contraction agrees with a separate
PySCF JK-driver contraction to about machine precision. A mixed singlet/triplet CI
test exercises the zero-Clebsch--Gordan fallback by constructing a normalized
`M_S=+1` triplet component with PySCF determinant ladder operations.

The provider is deliberately static-only. Calls to trajectory `components(q)`,
snapshot evaluation, or cross-geometry overlap methods raise. The v0.23.2 NAC
mapping (`state=(i,j)` for internal `d[i,j]`) is untouched. Prism was useful as an
external development reference but is neither imported nor distributed; the release
uses only PySCF plus the project's own spin algebra.

## What v0.24.0 adds

The first external intake target is OpenMolcas 26.06 RASSCF/CASPT2/RASSI-SO with
SEWARD AMFI for neutral water in one singlet plus one complete triplet. The exact
basis, active space, scalar-relativistic method, state order, units, and zero-field
convention are fingerprinted.

One reference and 54 displaced records cover all 9 Cartesian coordinates at three
step sizes. Native inputs, outputs, HDF5 files, exported component matrices, overlap
evidence, convergence flags, independent validation, environment, and manifest are
SHA-256 bound. Raw validation observations also require in-bundle content-addressed
blobs. Cross-geometry contractions are retained while operator derivatives use their
separately certified unitary polar transports.

The supplied generator creates a conspicuously synthetic protocol fixture. It passes
protocol tests but has a mandatory non-execution marker and a deliberately non-HDF5
placeholder, so it cannot be relabeled as external evidence. Production frozen-
snapshot dynamics requires a typed admission proof and rejects this fixture. Native
OpenMolcas HDF5/text numerical cross-parsing is not yet implemented, so that admission
prerequisite is hard-coded false even for a structurally complete submitted bundle.

## What v0.23.3 hardens

For a retained finite electronic manifold, the raw cross-geometry overlap is a
physical contraction and is kept as evidence. Coefficient transport uses its
separately certified unitary polar factor. Singular-value retention, condition
number, and principal angle decide trajectory readiness without redefining
physical consistency.

Replay format 2 stores both objects, the exact quality policy, singular values,
and a fingerprinted NAC convention. Legacy data need an exact migration
attestation; unknown or wrong-sign NAC records are quarantined without automatic
repair. Provider identities bind provenance, NAC mapping, overlap/transport
contracts, policy, and replay version for cache and checkpoint safety.

Complete singlet/triplet and Kramers-doublet manifolds are tracked through
projectors under independent endpoint gauges, with leakage and time-reversal
controls. The SOC matrix convention now freezes operator treatment, prefactor,
scalar relativity, state order, units, fixed-frame derivative semantics, complete
multiplets, and zero-field assumptions. Runtime reports separately distinguish
the canonical byte-locked build from supported scientific environments.

## What v0.23.2 corrects and validates

The optional PySCF target is now exercised in a hash-locked 2.13.1 runtime. A real
H3+ SA-CASSCF(2e,3o)/STO-3G calculation validates energies, analytic gradients,
nonadiabatic couplings, and many-electron overlaps for three singlet roots.

The internal derivative coupling `d[i,j]=<Phi_i|d Phi_j>` now uses PySCF
`state=(i,j)`, as certified by phase-aligned overlap central differences. The earlier
`state=(j,i)` adapter mapping had the opposite sign. Full-overlap derivatives use
`use_etfs=False`; ETF-corrected values remain a separate diagnostic.

Cross-geometry overlaps in a retained finite state manifold are now validated as
physical contractions: self identity and adjoint reciprocity are mandatory, while
cross singular values may be below one because omitted roots carry amplitude. Values
above one remain rejected.

Runtime admission is also hardened with structure-before-capability validation,
callable-interface checks, exact method identity, seven convergence stages, and typed,
trust-anchored parser/execution proof. These controls remain fail closed: they do not
admit a molecular-SOC backend.

## What v0.23.1 hardens

Every trajectory, basis, method, and frame calculation receives a fingerprinted input,
output, coordinate, method identity, and independent SCF/correlated/SOC/derivative/
overlap convergence receipt. Reference errors, convergence changes, and frame residuals
are recomputed from stored observations rather than accepted as summary assertions.

Tracking is derived on a connected record graph using minimum singular values within
complete physical manifolds and spectral-norm leakage into competing manifolds. This
keeps the test meaningful for triplet spaces and Kramers doublets where individual
components can rotate inside a degenerate subspace.

Hashes prove byte identity but not scientific meaning. External or live admission also
requires an executable method-specific raw-artifact parser whose identity matches the
runtime attestation. v0.23.1 proves that a synthetic dataset relabelled as external no
longer passes on summarized evidence alone.

## What v0.23.0 adds

Molecular SOC capabilities are split into `static_soc` and a derived
`trajectory_ready` tier. Moving-nuclear use requires spin-free and SOC physical
derivatives, derivative connections, and cross-geometry overlaps; a static matrix is
rejected.

Every real source must freeze the method, basis, electron count, SOC and relativistic
operators, derivative method, atomic/isotope/geometry identity, calculation-input hash,
environment hash, units, coordinate definition, and tracking policy. Real admission
also requires five independent evidence families: reference agreement, basis
convergence, method convergence, translational/rotational invariance, and quantitative
state-tracking quality.

The deterministic replay format stores exact component operators, connections, mass
matrices, every pair overlap, per-record convergence, time reversal, projectors, and
complete fingerprinted provenance. Protocol success and real-backend admission are
reported separately, so a fixture can validate the machinery without becoming a
physical-accuracy claim.

## What v0.22.1 corrects

The provider audit now verifies the transported full-matrix derivatives of
`H_spin_free` and `H_soc` separately for every nuclear coordinate. A regression fixture
proves that equal and opposite component errors are rejected even when the total K and
a sampled scalar force are unchanged.

SOC admission now binds the numerical time-reversal representation and physical
projectors into provenance and requires a single electron-parity/charge sector,
complete multiplets, time-reversal unitarity and square, complete projectors, and zero
external magnetic field. The audit is independent of provider-specific `.config` data
and is exercised on a three-state, two-coordinate provider.

The one-dimensional exact-grid oracle now reads mass from the operator contract,
requires a fixed frame and constant scalar mass, precomputes static split operators, and
always records the final step. Physical SOC convergence is additionally demonstrated
over 1/3/5-Gaussian bases and a four-level sparse-threshold ladder.

## What v0.22.0 adds

The even-electron model uses

$$
\{|S\rangle,|T_{-1}\rangle,|T_0\rangle,|T_{+1}\rangle\},
$$

with a time-reversal-constrained complex singlet–triplet SOC vector. The odd-electron
model uses two complete doublets,

$$
\{|D_1,+\tfrac12\rangle,|D_1,-\tfrac12\rangle,
|D_2,+\tfrac12\rangle,|D_2,-\tfrac12\rangle\},
$$

with a quaternionic SOC block that preserves twofold Kramers degeneracy at zero magnetic
field. Both models provide explicit

$$
H=H_{\mathrm{sf}}+H_{\mathrm{SOC}},\qquad
K_a=K_{a,{\mathrm{sf}}}+K_{a,{\mathrm{SOC}}}.
$$

Physical projector populations, time-reversal representations, SOC forces, differential
providers, exact-grid dynamics, moving complex gauges, and checkpoint identity are all
tested. Spin-free operation remains a permanent mode.

## What v0.21.4 adds

Pointwise structural validity is necessary but not sufficient. The new provider audit
uses centered geometry displacements and provider overlaps to transport neighboring
Hamiltonians back to the center frame. It independently compares the resulting
derivatives with K and the overlap derivative with D. Fixtures with Hermitian but wrong
K, or anti-Hermitian but wrong D, are rejected.

Checkpoint/restart is now a versioned numerical contract. Checkpoints carry the complete
propagated state needed for deterministic continuation, a canonical SHA-256 integrity
digest, and strict provider/settings fingerprints. Dense, sparse, moving-complex-frame,
zero-block density-guidance, and adaptive insertion/pruning restarts are exercised.

The zero-SOC rehearsal explicitly composes complex zero H_SOC and K_SOC with the
spin-free operators. This validates the future optional decomposition without adding
or claiming physical spin-orbit physics.

## Release validation highlights

Complex molecular operator covariance:

```text
H covariance error:
0.0

maximum dH covariance error:
0.0

force covariance error:
8.673617379884035e-19
```

Full block matrix covariance under a coordinate-dependent complex $U(2)$ gauge:

```text
S error:
1.7399964775418827e-16

H error:
1.7933278321452532e-16

T error:
1.7899468791987713e-16

maximum edge-score change:
5.551115123125783e-17
```

Time-dependent gauge-equivalent propagation converges at essentially second order:

| dt | Steps | Gauge-mapped coefficient error | Norm drift |
|---:|---:|---:|---:|
| 0.02000 | 5 | 2.047450e-08 | 6.570370e-09 |
| 0.01000 | 10 | 5.118673e-09 | 1.642602e-09 |
| 0.00500 | 20 | 1.279672e-09 | 4.106520e-10 |

Observed orders:

```text
[1.999986653128339, 1.9999961654805536]
```

The framework is also tested for **2, 4, and 8 electronic states** with no two-state
hard-coding:

| Electronic states | Total dimension | Active Gaussian edges | H nnz | H density |
|---:|---:|---:|---:|---:|
| 2 | 48 | 45 | 456 | 0.197917 |
| 4 | 96 | 45 | 1824 | 0.197917 |
| 8 | 192 | 45 | 7296 | 0.197917 |

## Dynamic sparse topology

The block-sparse graph is not assumed to be static. A six-Gaussian crossing stress test
produces both edge creation and edge deletion:

```text
total entered edges:
15

total exited edges:
9

maximum active edges:
15
```

## Install

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e ".[dev]"
```

For the validated PySCF paths, install the pinned optional dependency with
`python -m pip install -e ".[pyscf,dev]"`; the exact CPython 3.12/Linux x86-64 wheel
hashes for this release are in
`requirements-pyscf-v253-linux-x86_64-py312.txt`.

On Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -U pip
python -m pip install -e ".[dev]"
```

Run the cumulative regression suite:

```bash
python -m pytest -q
```

Recompute the v0.25.3 evidence and cumulative campaign in the pinned PySCF
environment:

```bash
python examples/137_recompute_v0253_controlled_basis.py
python examples/138_recompute_v0253_campaign.py
```

## Documentation

Start here:

1. [`docs/01_MATHEMATICAL_FOUNDATIONS.md`](docs/01_MATHEMATICAL_FOUNDATIONS.md)
2. [`docs/02_COMPLEX_GAUGE_FRAMEWORK.md`](docs/02_COMPLEX_GAUGE_FRAMEWORK.md)
3. [`docs/03_BLOCK_SPARSE_ALGORITHM.md`](docs/03_BLOCK_SPARSE_ALGORITHM.md)
4. [`docs/04_ALGORITHMIC_COMPLEXITY.md`](docs/04_ALGORITHMIC_COMPLEXITY.md)
5. [`docs/05_VALIDATION_STRATEGY.md`](docs/05_VALIDATION_STRATEGY.md)
6. [`docs/06_MOLECULAR_BACKENDS.md`](docs/06_MOLECULAR_BACKENDS.md)
7. [`docs/07_SCOPE_AND_ROADMAP.md`](docs/07_SCOPE_AND_ROADMAP.md)
8. [`docs/08_PRE_SOC_INTEGRATION_HARDENING.md`](docs/08_PRE_SOC_INTEGRATION_HARDENING.md)
9. [`docs/09_SOC_CONTRACT_FREEZE.md`](docs/09_SOC_CONTRACT_FREEZE.md)
10. [`docs/10_DIFFERENTIAL_AND_RESTART_CERTIFICATION.md`](docs/10_DIFFERENTIAL_AND_RESTART_CERTIFICATION.md)
11. [`docs/11_PHYSICAL_ANALYTIC_SOC.md`](docs/11_PHYSICAL_ANALYTIC_SOC.md)
12. [`docs/12_CORRECTIVE_SOC_HARDENING.md`](docs/12_CORRECTIVE_SOC_HARDENING.md)
13. [`docs/13_MOLECULAR_SOC_ADMISSION.md`](docs/13_MOLECULAR_SOC_ADMISSION.md)
14. [`docs/14_RAW_EVIDENCE_ADMISSION.md`](docs/14_RAW_EVIDENCE_ADMISSION.md)
15. [`docs/15_PYSCF_RUNTIME_AND_OVERLAP_CORRECTIONS.md`](docs/15_PYSCF_RUNTIME_AND_OVERLAP_CORRECTIONS.md)
16. [`docs/16_TRANSPORT_AND_COMPATIBILITY_HARDENING.md`](docs/16_TRANSPORT_AND_COMPATIBILITY_HARDENING.md)
17. [`docs/17_OPENMOLCAS_EXTERNAL_SOC_INTAKE.md`](docs/17_OPENMOLCAS_EXTERNAL_SOC_INTAKE.md)
18. [`docs/18_PYSCF_STATIC_MOLECULAR_SOC.md`](docs/18_PYSCF_STATIC_MOLECULAR_SOC.md)
19. [`docs/19_PYSCF_CONNECTED_GEOMETRY_SOC.md`](docs/19_PYSCF_CONNECTED_GEOMETRY_SOC.md)
20. [`docs/20_VARIATIONAL_SOC_DYNAMICS.md`](docs/20_VARIATIONAL_SOC_DYNAMICS.md)
21. [`docs/21_MULTIGAUSSIAN_TDVP.md`](docs/21_MULTIGAUSSIAN_TDVP.md)
22. [`docs/22_ADAPTIVE_MULTIGAUSSIAN_TDVP.md`](docs/22_ADAPTIVE_MULTIGAUSSIAN_TDVP.md)
23. [`docs/23_CONTROLLED_BASIS_ADAPTATION.md`](docs/23_CONTROLLED_BASIS_ADAPTATION.md)
24. [`docs/RELEASE_HISTORY.md`](docs/RELEASE_HISTORY.md)

Release-specific files are also retained at the repository root.

## Scientific scope

The appropriate description of v0.25.3 is:

> **A complex representation-neutral, block-sparse Gaussian nonadiabatic dynamics
> research framework with physical analytic SOC, validated real PySCF spin-free
> SA-CASSCF data, certified finite-manifold transport, complete even/odd manifolds,
> a strict fail-closed OpenMolcas RASSI-SO external-snapshot intake protocol, and
> direct-JK PySCF BP-SOMF state-interaction matrices plus endpoint-bound,
> polar-transported OH bond-coordinate SOC differential evidence, restricted-TDVP
> analytic-SOC trajectories, one-dimensional coupled adaptive
> log-width/quadratic-chirp multi-Gaussian McLachlan dynamics, and now a controlled
> residual-driven topology lifecycle with SVD-projected spawning, coefficient-only
> newborn activation, and conservative prune/merge receipts.**

It is **not** a production implementation of general AIMS and does **not** admit
multidimensional spawning, unrestricted branching, full-matrix adaptive widths,
coordinate-dependent electronic frames in the v0.25.3 solver, a trajectory-ready
external/live molecular SOC backend, a full Cartesian or analytic molecular-SOC
derivative, a continuous
physical molecular SOC derivative connection, an OpenMolcas execution, external
magnetic fields, or general molecular SOC dynamics accuracy.

The framework is designed so that ordinary spin-free nonadiabatic dynamics remains a
first-class use case even after an optional SOC-capable backend is added later.
