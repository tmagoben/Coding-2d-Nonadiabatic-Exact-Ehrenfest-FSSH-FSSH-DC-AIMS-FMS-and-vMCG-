# Physical analytic spin-orbit coupling

## Stable v0.22.0 contract

v0.22.0 introduces the first nonzero physical SOC operators in the project. A provider
must supply a complete, fixed electronic model space; explicit provenance; the total
Hermitian Hamiltonian; physical Hermitian derivative operators; derivative connections;
and cross-geometry overlaps. The analytic providers additionally expose the decomposition

$$
H=H_{\mathrm{sf}}+H_{\mathrm{SOC}},\qquad
K_a=K_{a,{\mathrm{sf}}}+K_{a,{\mathrm{SOC}}}.
$$

The Gaussian engine consumes the total operators and remains unaware of which physical
terms produced them.

## Model families

The even-electron reference contains a singlet and all three triplet components. Its
time-reversal representation squares to $+I$. The odd-electron reference contains
two complete doublets, uses a quaternionic inter-doublet SOC block, and has a
time-reversal representation that squares to $-I$. At zero field this enforces
Kramers pairing.

The two references are separate models. Electron-number parity is provenance, not a
dynamical label, and the release never combines even- and odd-electron sectors.

## Required invariants

Every physical-SOC path is checked for:

1. Hermiticity of spin-free, SOC, and total $H$ and $K_a$;
2. exact operator composition;
3. correct time-reversal action and square;
4. complete, orthogonal physical projectors;
5. centered finite-difference agreement for SOC forces;
6. cross-geometry consistency of $H/K/D$;
7. general complex-gauge covariance;
8. complete provenance identity across cache and restart boundaries.

For $c'=G^\dagger c$, the non-obvious transformation is

$$
J'=G^\dagger JG^*,
$$

while a physical projector transforms as $P'=G^\dagger PG$. The derivative
connection receives the usual inhomogeneous gauge term; the physical derivative
operator $K_a$ does not.

## Validation hierarchy

The release first tests pointwise structure, then cross-geometry differentials, then
symmetry and gauge covariance, then independent exact-grid propagation, and finally
SOC-active Gaussian restart. Negative fixtures prove that plausible but wrong data are
rejected. All 21 v0.21.4 gates remain part of the 53-gate campaign.

## Supported claims

The supported description is an analytic physical-SOC Gaussian nonadiabatic dynamics
framework with complete singlet/triplet and Kramers-doublet references, physical SOC
derivatives, exact-grid validation, gauge-covariant observables, and deterministic
restart.

The release does not claim an ab-initio SOC backend, PySCF SOC runtime, molecular SOC
accuracy, external magnetic fields, or production AIMS. Those require a separate
backend-validation milestone. Spin-free propagation remains supported indefinitely.

See `docs/releases/v0.22.0/V220_DERIVATIONS.md`,
`docs/releases/v0.22.0/V220_VALIDATION.md`, and
`docs/releases/v0.22.0/V220_ALGORITHM_COMPLEXITY.md` in the version archive for the
detailed equations, measured values, and scaling analysis.
