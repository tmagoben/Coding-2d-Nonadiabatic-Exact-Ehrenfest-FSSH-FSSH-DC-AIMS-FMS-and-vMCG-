# Gaussian Nonadiabatic Dynamics — Complete Mathematical Walkthrough through v0.29.2

This document explains the project in the order the ideas were introduced. Atomic units are used unless noted, so \(\hbar=1\).

## 0. The problem we are solving

The molecular TDSE is

\[
i\partial_t\Psi(R,t)=\hat H\Psi(R,t).
\]

For several electronic states,

\[
\hat H = \hat T_N\mathbf 1 + H_{\rm el}(R),
\qquad
\hat T_N=-\frac12\nabla_R^TM^{-1}\nabla_R.
\]

The exact nuclear-electronic wavefunction is expensive because it is a multicomponent function over the full nuclear configuration space. The project replaces the nuclear field by a controlled moving Gaussian basis while preserving electronic coherence, nonorthogonality, gauge covariance, spin-orbit structure, and eventually curved electronic bundles.

---

## v0.1 — Gaussian mechanics and moving nonorthogonal bases

### Single Gaussian

Start with

\[
g(x;q,p,\alpha)=\left(\frac{\alpha}{\pi}\right)^{1/4}
\exp\left[-\frac\alpha2(x-q)^2+ip(x-q)\right].
\]

It has

\[
\langle x\rangle=q,\quad \sigma_x^2=\frac1{2\alpha},\qquad
\langle p\rangle=p,\quad \sigma_p^2=\frac\alpha2.
\]

A thawed Heller Gaussian generalizes the exponent to a complex quadratic width. Expanding the potential locally to second order closes the Gaussian equations and is exact for a quadratic potential.

### Multiple moving Gaussians

The wavefunction becomes

\[
\Psi(x,t)=\sum_I C_I(t)g_I(x,t).
\]

The basis is nonorthogonal:

\[
S_{IJ}=\langle g_I|g_J\rangle.
\]

Because the basis itself moves,

\[
T_{IJ}=\langle g_I|\dot g_J\rangle
\]

must enter. Projection of the TDSE gives

\[
iS\dot C=(H-iT)C.
\]

This is the core moving-basis identity that later sparse, gauge, and molecular versions inherit.

---

## v0.2 — Two-state nonadiabatic dynamics, adiabatic couplings, and spawning

Introduce a two-state diabatic Hamiltonian

\[
H_d=-\frac1{2M}\frac{d^2}{dx^2}I_2+
\begin{pmatrix}V_{11}&V_{12}\\V_{12}&V_{22}\end{pmatrix}.
\]

Diagonalizing \(V_d(x)\) gives adiabatic states \(|\phi_a(x)\rangle\). Their position dependence generates derivative couplings

\[
d_{ab}(x)=\langle\phi_a|\partial_x\phi_b\rangle.
\]

A phase change \(|\phi_a\rangle\to e^{i\chi_a(x)}|\phi_a\rangle\) changes the raw derivative coupling, so phase/gauge continuity must be tracked.

A Gaussian on one surface cannot generally represent branching after a crossing. v0.2 therefore introduces dynamic spawning: a child Gaussian is created in a new electronic channel while coherence with the parent is retained.

---

## v0.3 — Electronic-structure provider contract

The dynamics should ask only for physical electronic information at geometry \(R\):

\[
E_a(R),\qquad \nabla_RE_a(R),\qquad d_{ab}(R).
\]

If dynamics uses generalized coordinate \(q\) with Cartesian geometry \(R(q)\),

\[
\frac{dE_a}{dq}=\nabla_RE_a\cdot\frac{dR}{dq}.
\]

The provider abstraction separates dynamics from the electronic-structure program and establishes the bridge toward direct dynamics.

---

## v0.4 — Two-dimensional conical intersections

The nuclear coordinate becomes \(R=(x,y)^T\):

\[
H_d=-\frac1{2M}(\partial_x^2+\partial_y^2)I_2+V_d(x,y).
\]

A conical intersection requires at least a two-dimensional branching plane. Derivative couplings become vectors, Berry phase becomes observable around loops, and spawned Gaussians must carry vector positions and momenta.

---

## v0.5 — Backend-driven molecular direct dynamics

The generalized-coordinate provider is connected to a Cartesian backend. Electronic quantities are computed in \(3N\)-dimensional Cartesian space and projected into the active nuclear coordinates. The architecture becomes

\[
\text{Gaussian dynamics}\to\text{provider}\to\text{electronic backend}.
\]

This is where analytic models and PySCF-style backends can share one dynamics interface.

---

## v0.6 — State tracking and gauge continuity

Energy ordering is not state identity near crossings. Cross-geometry many-electron overlaps

\[
O_{ab}(R,R')=\langle\Phi_a(R)|\Phi_b(R')\rangle
\]

are used to match roots and align phases/subspaces. This prevents a smooth physical state from becoming numerically discontinuous because independent diagonalizations arbitrarily permute or rephase eigenvectors.

---

## v0.7 — Gauge graph and Wilson-loop viewpoint

Electronic points become nodes of a graph; overlaps define edges. Polar decomposition produces unitary link variables. Along a path \(\gamma\), transport is the ordered product

\[
U_\gamma=U_{01}U_{12}\cdots U_{n-1,n}.
\]

Around a closed loop,

\[
W_\gamma=\prod_{(ij)\in\gamma}U_{ij}
\]

is a Wilson loop whose spectrum captures Berry/non-Abelian holonomy. TBF electronic vectors are transported to common pair frames before pair matrix elements are built.

---

## v0.8 — Time-dependent gauge graph and AIMS-like propagation

The graph becomes time dependent. Temporal overlap transports electronic amplitudes between successive frames. The moving-basis equation remains

\[
iS\dot C=(H-iT)C,
\]

but \(T\) now includes both Gaussian motion and frame motion. Dynamic zero-amplitude spawning is introduced so child amplitudes grow through coupling rather than by arbitrary amplitude assignment.

---

## v0.9 — Convergence, SPA layers, basis conditioning, and projection pruning

Electronic matrix elements are approximated by saddle-point/Taylor layers. SPA0 evaluates at the pair centroid. SPA1 adds first-order corrections using Gaussian moments.

Basis stability is measured through singular/eigenvalues of \(S\). Redundant functions are removed by Hilbert-space projection rather than deletion. Spawning becomes timestep-aware through integrated coupling exposure, and exact-grid comparisons are formalized.

---

## v0.10 — Observable-correct benchmark campaigns

The project stops judging methods only by populations. The reduced electronic density matrix is

\[
\rho_e=\operatorname{Tr}_N|\Psi\rangle\langle\Psi|.
\]

From it we compute populations, coherences, purity

\[
\mathcal P=\operatorname{Tr}(\rho_e^2),
\]

and entropy. Exact and Gaussian calculations are compared in a common fixed frame. Wigner initial sampling and explicit multi-axis convergence campaigns become part of the method.

---

## v0.11 — Basis-completeness campaign

The basis size itself becomes a convergence variable. Instead of assuming the dynamics error is solely time-integration error, the project separates representation error from propagation error by growing the initial Gaussian basis and measuring convergence of projected observables.

---

## v0.12 — Representation-consistent spinor Gaussians

A center-frozen electronic approximation \(g_I(R)\Phi_a(q_I)\) can badly misrepresent entanglement. v0.12 reconstructs coordinate-dependent electronic spinors consistently. This exposes why purity/coherence can be wrong even before propagation. The benchmark compares like representations with like representations.

---

## v0.13 — Residual-driven basis refinement

Define the instantaneous TDSE/Galerkin defect

\[
r=\dot\Psi+iH\Psi.
\]

A candidate basis function \(\phi_c\) is useful if it captures a component of this residual that the current basis cannot represent. Candidate dictionaries are ranked by predicted residual reduction, then actual reprojection checks verify the prediction. This is the conceptual ancestor of v0.29.2.

---

## v0.14 — Fully adaptive defect control

Residual checks occur repeatedly in time. Hysteresis prevents add/remove chatter; cooldowns prevent rapid topology oscillation. Pruning uses exact leave-one-out nonorthogonal projection loss. Basis-budget replacement and emergency conditioning policies make adaptation a controlled dynamical process instead of one-time basis construction.

---

## v0.15 — Shared pair cache and incremental algebra

Gaussian pair operations are cached. Since

\[
S_{JI}=S_{IJ}^*,\qquad H_{JI}=H_{IJ}^*,
\]

only canonical pairs need expensive factorization. Candidate acceptance appends blocks incrementally rather than rebuilding dense matrices. Physics is intended to remain equivalent while computational work is reduced.

---

## v0.16 — Persistent locality graph and sparse propagation

Only sufficiently important Gaussian pairs become graph edges. Sparse \(S,H,T\) are assembled, and propagation uses sparse Cayley/midpoint-type solves. The graph has entry/exit hysteresis. The key methodological shift is from dense \(O(N_g^2)\) pair work toward locality-controlled sparse work.

---

## v0.17 — Sparse error control

Edge selection is no longer overlap-only. Importance combines projected contributions to \(S,H,T\). A global omitted-edge budget controls sparsification, periodic dense audits estimate error, and thresholds automatically relax if the sparse approximation violates the audit.

---

## v0.18 — Dynamics/convergence completeness

The full spinor wavefunction is reconstructed and compared using global-phase-aligned fidelity and \(L^2\) error. Nuclear densities, centroids, covariances, timestep convergence, basis convergence, edge-budget convergence, and adaptation-threshold sensitivity are all separated. This turns “looks stable” into an explicit convergence ladder.

---

## v0.19 — Molecular direct-dynamics integration

The provider becomes overlap-capable and cache-aware. Root matching uses polynomial-time assignment rather than factorial permutation enumeration. Center and centroid requests share a molecular gauge graph. Spawned dynamics becomes provider neutral, so an analytic molecular backend and a PySCF bridge can exercise the same interfaces.

---

## v0.20 — End-to-end sparse molecular machinery

A buffered KD-tree finds nearby electronic geometries; a persistent locality graph limits pair work. Sparse molecular \(S/H/T\), metric-compatible moving-basis connections, sampled omitted-edge audits, and search-floor relaxation are combined. This is the first full sparse molecular architecture.

---

## v0.21.0 — Complex, block, gauge-ready generalization

Electronic amplitudes become arbitrary complex blocks. A local gauge is now a general unitary matrix

\[
U(R)\in U(N_s).
\]

Operators transform as

\[
H' = U^\dagger H U.
\]

Subspaces are aligned by unitary Procrustes/polar maps and Wilson-loop spectra validate non-Abelian gauge structure. This release intentionally adds no physical spin yet.

---

## v0.21.2 — Pre-SOC integration hardening

Unequal-width complex pair algebra is generalized. Nuclear guidance is representation neutral. Coefficients and nuclear variables are propagated self-consistently with predictor/corrector structure. Adaptive block insertion/pruning and generic complex observables are hardened before SOC is allowed.

---

## v0.21.3 — SOC contract freeze

The electronic Hamiltonian is decomposed conceptually as

\[
H=H_{\rm sf}+H_{\rm SOC},
\qquad
K_a=\partial_aH=K_{a,\rm sf}+K_{a,\rm SOC}.
\]

Strict Hermiticity, model-space, provenance, unit, and complete-multiplet contracts are introduced. Degeneracy-safe density-matrix guidance replaces basis-dependent “lowest eigenvector” fallback logic.

---

## v0.21.4 — Differential provider and deterministic restart certification

Pointwise-valid providers can still have wrong derivatives. Centered cross-geometry audits therefore check

\[
K_a(R)\approx\frac{H(R+h e_a)-H(R-h e_a)}{2h}
\]

and corresponding connection/overlap structure. Deterministic restart records provider/settings fingerprints, UIDs, graph edges, coefficients, and guidance state.

---

## v0.22.0 — First physical SOC models

Physical nonzero spin-orbit Hamiltonians are introduced in complete spin manifolds. Even-electron singlet/triplet and odd-electron Kramers-doublet models obey time-reversal structure. The total operator is

\[
H(R)=H_{\rm sf}(R)+H_{\rm SOC}(R),
\]

with physical derivatives

\[
\partial_aH=\partial_aH_{\rm sf}+\partial_aH_{\rm SOC}.
\]

An independent exact-grid matrix-wavefunction calculation provides the first SOC reference.

---

## v0.22.1 — Corrective SOC hardening

SOC conventions, complete multiplets, time-reversal/Kramers structure, exact-grid implementation, and Gaussian convergence are tightened. Sparse-threshold convergence to dense SOC dynamics is explicitly checked rather than inferred.

---

## v0.23.0 — Molecular SOC provider contract

The analytic SOC contract is generalized toward molecular providers: complete spin model spaces, SOC matrix provenance, physical derivative meaning, and arbitrary-geometry snapshot expectations are specified before trusting a live backend.

---

## v0.23.1 — External evidence hardening

Provider claims are separated from synthetic fixtures. External/molecular evidence must carry provenance that cannot be satisfied merely by relabeling an analytic test dataset.

---

## v0.23.2 — NAC convention correction

Derivative-coupling conventions and cross-geometry overlap/connection semantics are corrected and revalidated so signs, anti-Hermiticity, and coordinate derivatives agree consistently.

---

## v0.23.3 — Runtime/profile hardening

Runtime profiles and provider failure boundaries are made explicit. Unsupported live SOC or overlap operations fail closed instead of silently substituting synthetic behavior.

---

## v0.24.0 — External SOC intake

A formal intake layer admits externally supplied SOC operators only after structure, units, provenance, complete model-space, and derivative audits pass.

---

## v0.24.1 — Runtime admission hardening

Static or differential provider capabilities are distinguished from true arbitrary-geometry trajectory capability. Missing snapshot/overlap functionality is rejected rather than extrapolated.

---

## v0.24.2 — Gauge preview / moving-frame bridge

The code begins representing coordinate-dependent electronic frames explicitly and previews their connection structure, without yet claiming a fully variational moving-frame Gaussian trajectory.

---

## v0.25.0 — First canonical Gaussian TDVP trajectory

The project moves from partly prescribed Gaussian nuclear motion to a variationally determined Gaussian trajectory. For real parameters \(\theta\), McLachlan minimizes

\[
\mathcal F(\dot\theta)=\|\dot\Psi+iH\Psi\|^2.
\]

Stationarity gives

\[
G\dot\theta=b,
\]

with

\[
G_{\mu\nu}=\Re\langle\partial_\mu\Psi|\partial_\nu\Psi\rangle,
\qquad
b_\mu=\Im\langle\partial_\mu\Psi|H\Psi\rangle.
\]

This is the central equation of the modern project.

---

## v0.25.1 — Coupled multi-Gaussian TDVP and compatible null spaces

The ansatz becomes

\[
\Psi(x)=\sum_{I=1}^{N_g}g_I(x)\sum_{a=1}^{N_s}C_{Ia}|a\rangle.
\]

All coefficients, centers, and momenta are varied together. Because duplicate or nearly redundant packets can make \(G\) singular, full SVD is used:

\[
G=U\Sigma V^T,
\qquad
G^+=V\Sigma^+U^T.
\]

Discarded singular directions are accepted only if the RHS is compatible with the null space. Time propagation is fully implicit midpoint:

\[
\theta_{n+1}-\theta_n
=\Delta t\,v\!\left(\frac{\theta_n+\theta_{n+1}}2\right).
\]

---

## v0.25.2 — Thawed widths and chirps

Each 1D packet is

\[
g_I(x)=\left(\frac{\alpha_I}{\pi}\right)^{1/4}
\exp\left[-\frac{\alpha_I-i\beta_I}{2}(x-q_I)^2+ip_I(x-q_I)\right].
\]

Use \(\eta_I=\log\alpha_I\) so \(\alpha_I>0\) structurally. Width and chirp are conjugate shape variables. Tangents include

\[
\partial_{\eta_I}g_I=\left(\frac14-\frac{\alpha_I}{2}y_I^2\right)g_I,
\qquad
\partial_{\beta_I}g_I=\frac i2y_I^2g_I.
\]

Gaussian cross moments through total degree four close the metric and RHS analytically.

---

## v0.25.3 — Controlled basis lifecycle

At adaptation checkpoints the order is frozen:

\[
\text{merge}\rightarrow\text{prune}\rightarrow\text{spawn}\rightarrow\text{no-op}.
\]

A candidate is orthogonalized against the existing basis and scored by residual capture. Every topology change uses SVD projection. A newborn inserted by exact projection has zero coefficient, so its shape tangents vanish; therefore only its electronic coefficients evolve until its population reaches an activation gate.

---

## v0.26.0 — Multidimensional reference-first TDVP

For diagonal multidimensional shape matrices,

\[
g_I(R)=\prod_\mu\left(\frac{\alpha_{I\mu}}\pi\right)^{1/4}
\exp\left[-\frac{\alpha_{I\mu}-i\beta_{I\mu}}2(R_\mu-q_{I\mu})^2+i p_{I\mu}(R_\mu-q_{I\mu})\right].
\]

A 2D exact FFT/Strang CI+SOC solver becomes an independent reference. Newborn activation now requires all three:

\[
P_I\ge10^{-6},\qquad \kappa_{\rm metric}\le10^8,
\qquad \|v_{\rm active}\|\le100\|v_{\rm dormant}\|.
\]

---

## v0.27.0 — Full correlated Gaussian widths

The packet becomes

\[
g_I(R)=\left(\frac{\det\Gamma_I}{\pi^D}\right)^{1/4}
\exp\left[-\frac12x_I^T(\Gamma_I-iB_I)x_I+i p_I^Tx_I\right],
\]

with

\[
\Gamma_I=\exp(E_I)\succ0,
\]

where \(E_I\) and \(B_I\) are full real symmetric matrices. Symmetric matrices are stored in a Frobenius-orthonormal `svec` basis. The width derivative is the Fréchet derivative

\[
\delta\Gamma=D\exp(E)[\delta E].
\]

This allows rotated/correlated packets and honest arbitrary orthogonal coordinate covariance.

---

## v0.28.0 — Coordinate-dependent electronic moving frames, flat connection

A coordinate-dependent electronic frame \(U(R)\) introduces a matrix connection

\[
A_\mu(R)=U^\dagger(R)\partial_\mu U(R),
\]

which is anti-Hermitian. The covariant derivative is

\[
D_\mu=\partial_\mu+A_\mu.
\]

v0.28 validates moving-frame covariance, discrete transport, lattice Hermiticity, and zero/flat-curvature consistency. The key conceptual shift is that electronic basis change is treated as geometry, not as a post-processing phase convention.

---

## v0.29.0 — Constant non-Abelian connection and curved Gaussian sections

A genuinely noncommuting constant connection is admitted. For packet center \(q_I\), define a straight Wilson section

\[
L_I(R)=\exp[-(R-q_I)^\mu D_\mu].
\]

The packet is no longer merely \(g_I(R)c_I\); it is the electronic-fiber section

\[
g_I(R)L_I(R)c_I.
\]

Noncommuting generators produce nonzero curvature/holonomy even though the generators are spatially constant.

---

## v0.29.1 — Native curved full-width TDVP

The complete state is

\[
\Psi(R;\theta)=\sum_I g_I(R)L_I(R)c_I.
\]

The center tangent must differentiate both the nuclear Gaussian and the Wilson transporter:

\[
\frac{\partial\Psi}{\partial q_{I\mu}}
=g_I\left[((P_Ix_I)_\mu-ip_{I\mu})u_I
+\mathcal L_{\exp}(A_I,D_\mu)c_I\right],
\]

where \(P_I=\Gamma_I-iB_I\), \(u_I=L_Ic_I\), and \(\mathcal L_{\exp}\) is the Fréchet derivative of the matrix exponential.

The Hamiltonian kinetic term is handled covariantly,

\[
H=-\frac12(M^{-1})_{\mu\nu}D_\mu D_\nu+V(R),
\]

and the weak projection is

\[
\langle J_j|H\Psi\rangle
=\frac12(M^{-1})_{\mu\nu}\langle D_\mu J_j|D_\nu\Psi\rangle
+\langle J_j|V|\Psi\rangle.
\]

The same SVD McLachlan solve and implicit midpoint then operate on the curved manifold.

---

## v0.29.2 — Residual-driven curved basis adaptation

This checkpoint asks the correct spawning question:

> Does the candidate add a physical direction that the *entire existing curved TDVP tangent manifold* cannot already generate?

Let \(J\) be the current complete real tangent matrix and \(K_c\) the candidate's coefficient-tangent block,

\[
K_c=\{g_c(R)L_c(R)e_a\}_{a=1}^{N_s}.
\]

With quadrature weight operator \(W\), define

\[
M=\Re(J^\dagger WJ),\qquad
X=\Re(J^\dagger WK_c),\qquad
C=\Re(K_c^\dagger WK_c).
\]

Project candidate directions out of the current tangent space. The metric of the genuinely new directions is the Schur complement

\[
\boxed{S_c=C-X^TM^+X}.
\]

For the current McLachlan defect

\[
r=\dot\Psi+iH\Psi,
\]

define

\[
g=\Re(J^\dagger Wr),\qquad
q_c=\Re(K_c^\dagger Wr),
\]

and the residual component in the new candidate subspace

\[
\boxed{y_c=q_c-X^TM^+g}.
\]

The predicted squared-defect reduction is

\[
\boxed{\Delta_c=y_c^TS_c^+y_c}.
\]

A large \(\Delta_c\) means the candidate supplies a missing dynamical direction. A small \(\Delta_c\) means existing coefficient/center/momentum/width/chirp/Wilson-center tangents can already absorb that motion.

### Candidate geometry

Candidates inherit the correlated principal-axis rule. If

\[
\Gamma_Iu_a=\lambda_au_a,
\]

then trial centers/momenta are

\[
q'=q\pm\frac{c_q}{\sqrt{\lambda_a}}u_a,
\qquad
p'=p\pm c_p\sqrt{\lambda_a}u_a.
\]

Degenerate width eigenspaces are rejected because no unique intrinsic axis exists.

### Full curved projection

Every accepted topology change solves a full least-squares projection in the curved coefficient basis. If \(A\) contains all target coefficient columns and \(\psi\) is the source field,

\[
\min_c\|Ac-\psi\|_W^2
\]

is solved by SVD. Spawn therefore preserves the wavefunction and gives a zero-amplitude newborn to numerical precision.

### Merge and prune

Two packets are merge candidates only when their *curved electronic coefficient subspaces* are nearly identical. Their overlap block is

\[
S_{IJ}^{ab}=\int dR\,[g_IL_Ie_a]^\dagger[g_JL_Je_b].
\]

The minimum singular value of this block is used as a strict redundancy score. Merge/prune are accepted only if full curved projection loss and energy jump stay below their gates.

### Newborn activation

The newborn coefficients evolve immediately, but its shape coordinates remain dormant. Activation still requires

\[
P_I\ge10^{-6},\qquad
\kappa\le10^8,
\qquad
\frac{\|v_{\rm trial}\|}{\|v_{\rm dormant}\|}\le100.
\]

A second spawn is forbidden while a dormant newborn remains.

### Exact no-event boundary

If no lifecycle event occurs, v0.29.2 returns the existing state without numerical rewriting. This isolates basis-adaptation effects from propagation effects.

### Validation receipt

For the canonical development fixture:

\[
\Delta_{\rm predicted}=1.1783841637889198\times10^{-2},
\]

\[
\Delta_{\rm independent}=1.1783841637889186\times10^{-2},
\]

with relative disagreement

\[
1.03\times10^{-15}.
\]

Spawn projection loss is \(4.36\times10^{-31}\), fidelity is 1 to displayed precision, and the newborn coefficient norm is \(4.02\times10^{-16}\).

---

# What the latest project is doing in one continuous chain

1. Choose a complete electronic model space.
2. Obtain or define the electronic Hamiltonian and its physical coordinate derivatives.
3. Maintain state/subspace identity across geometry with overlaps and gauge transport.
4. Represent nuclear motion by nonorthogonal correlated Gaussian packets.
5. Attach a local electronic fiber section to each packet through Wilson transport.
6. Build the physical total wavefunction by coherently summing all curved Gaussian spinors.
7. Differentiate that full state with respect to every real variational parameter.
8. Build the McLachlan metric \(G=\Re(J^\dagger J)\).
9. Build the force/RHS \(b=\Im(J^\dagger H\Psi)\).
10. Solve \(G\dot\theta=b\) with a compatible-null full SVD.
11. Advance \(\theta\) by a nonlinear implicit-midpoint step.
12. Reconstruct norm, energy, width positivity, covariance, and solver receipts.
13. At an adaptation checkpoint, try merge first, then prune, then spawn.
14. For spawning, generate intrinsic correlated principal-axis candidates.
15. Turn each nuclear candidate into a curved Wilson section.
16. Project each candidate against the *full current tangent space*.
17. Compute its Schur-complement residual reduction \(\Delta_c\).
18. Reject candidates that add no independent tangent direction, are ill conditioned, or remove too little defect.
19. Select the largest admitted \(\Delta_c\), with deterministic tie-breaking.
20. Project the complete old curved wavefunction into the enlarged basis using full SVD.
21. Keep the newborn shape dormant while its coefficient grows.
22. Activate shape only after population, conditioning, and velocity-amplification gates pass.
23. Continue propagation and repeat the cycle.
24. Cross-check the continuum implementation against finite-difference and independent lattice/grid oracles.

The next planned layer after this checkpoint is coordinate-dependent connections/path families/patch transitions, followed by live molecular derivative-coupling and SOC trajectories. ML spawning should initially rank candidates but remain downstream of the exact physics oracle defined here.
