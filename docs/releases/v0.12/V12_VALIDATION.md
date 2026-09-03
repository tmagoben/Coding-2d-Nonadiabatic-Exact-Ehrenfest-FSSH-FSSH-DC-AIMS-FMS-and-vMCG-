# v0.12 Validation Contract

Version 0.12 changes the validation question.

Earlier releases compared a center-based Gaussian electronic representation directly
with a coordinate-dependent adiabatic exact wavepacket. v0.12 first verifies that the
two calculations represent the same initial state before attributing a discrepancy to
the dynamics.

The validation hierarchy is therefore

```text
representation
    ↓
projection
    ↓
Hamiltonian matrix elements
    ↓
moving-basis metric consistency
    ↓
projected-state dynamics
    ↓
original-target dynamics
```

---

## 1. Cumulative regression

Every retained v0.1-v0.11 test must remain passing.

This includes:

- Heller and exact-grid Gaussian propagation;
- moving nonorthogonal bases;
- variational Gaussian dynamics;
- spawning;
- 2D conical-intersection topology;
- PySCF SA-CASSCF API contracts;
- many-electron state tracking;
- gauge graphs and Wilson loops;
- moving graph propagation;
- SPA0/SPA1;
- basis conditioning/pruning;
- convergence campaigns;
- unequal-width Gaussian algebra;
- optimized-spawning-inspired child placement.

---

## 2. Exact analytic LVC Gaussian integrals

For unequal-width Gaussians, direct numerical 2D quadrature independently checks:

- overlap;
- gradient matrix element;
- kinetic matrix element;
- full LVC potential matrix element;
- total Hamiltonian matrix element.

The exact-LVC Hamiltonian must be Hermitian to numerical precision for mixed
electronic states and mixed Gaussian widths.

---

## 3. Moving electronic frame

For a center-following adiabatic state,

$$
\dot\Phi_j
=
\sum_i
\Phi_i
(\dot R\cdot d)_{ij}.
$$

The implementation is checked against a finite-difference derivative of the analytic
2D electronic frame away from the CI.

The full moving-basis identity

$$
\dot S=T+T^\dagger
$$

is checked independently by finite differences.

---

## 4. Generalized Cayley integrator

The v0.12 midpoint step

$$
\left[
S_m+\frac{\Delta t}{2}(iH_m+T_m)
\right]C_{n+1}
=
\left[
S_m-\frac{\Delta t}{2}(iH_m+T_m)
\right]C_n
$$

must reduce exactly to the ordinary fixed-basis Cayley propagator when

$$
S=I,
\qquad
T=0.
$$

The fixed-basis Cayley operator must be unitary for Hermitian $H$.

---

## 5. Local-diabatic transport

For the complete analytic two-state electronic frame, overlap transport around a small
step must preserve the same physical vector in the global diabatic basis.

The instantaneous-adiabatic reset is kept as a separate ablation and must not be
confused with parallel transport.

---

## 6. Spinor-complete electronic basis

For

$$
\Psi
=
\sum_{ka}
C_{ka}g_k|d_a\rangle,
$$

tests verify:

- $S=S^\dagger$;
- $H=H^\dagger$;
- correct block dimensions;
- one-Gaussian center-adiabatic initialization has unit generalized norm;
- reduced-density trace equals one after normalization;
- short fixed-bank propagation conserves the generalized norm.

Pair-preserving pruning removes an entire nuclear Gaussian together with all of its
electronic components.

---

## 7. Coordinate-dependent Born-Huang representation

For

$$
\Xi_i(R)=g_i(R)\Phi_{a_i}(R),
$$

the tests verify:

1. one TBF reconstructs the exact localized adiabatic initial packet on the same grid;
2. its initial reduced electronic density matches direct grid integration;
3. projected $S$ and $H$ are Hermitian;
4. the moving-basis matrix reproduces the finite-difference derivative of $S$;
5. short projected dynamics conserves the metric norm.

The kinetic operator used in this path is the same periodic FFT operator as the exact
2D benchmark.

---

## 8. Initial projection

The spinor-complete Gaussian projection solves

$$
SC=b.
$$

The validation suite requires a shifted multi-Gaussian bank to improve the projection
fidelity relative to one center Gaussian.

The projection output reports:

- wavefunction fidelity;
- relative residual;
- initial reduced-density error;
- overlap condition number.

---

## 9. Representation-consistent projected-state benchmark

This is the primary v0.12 release criterion.

For each initial Gaussian bank:

1. construct the exact target initial packet;
2. project it into the spinor-complete Gaussian bank;
3. propagate the projected wavefunction exactly on the 2D grid;
4. propagate the same projected state with the Gaussian method;
5. compare the two final reduced electronic density matrices.

Define

$$
\epsilon_{\mathrm{dyn}}
=
\|
\rho_G(t_f)
-
\rho_{\mathrm{projected}}^{\mathrm{exact}}(t_f)
\|_F.
$$

The release reference requires

$$
\boxed{
\epsilon_{\mathrm{dyn}}<10^{-3}.
}
$$

This criterion isolates propagation from initial-state representation.

---

## 10. Original-target benchmark

The Gaussian result is also compared with exact propagation of the intended
coordinate-dependent adiabatic packet:

$$
\epsilon_{\mathrm{target}}
=
\|
\rho_G(t_f)
-
\rho_{\mathrm{target}}^{\mathrm{exact}}(t_f)
\|_F.
$$

The v0.12 reference thresholds are:

```text
initial reduced-density representation error <= 0.05
projected-state dynamics density error        <= 0.001
original-target full-density error            <= 0.05
original-target population L2 error           <= 0.05
coherence phase error                         <= 0.01 rad
generalized norm drift                        <= 1e-4
overlap condition number                      <= 1e5
```

These are release-regression criteria for the analytic benchmark, not universal
chemical-accuracy standards.

---

## 11. Coherence diagnostics

The release reports separately:

- complex $\rho_{01}$;
- $|\rho_{01}|$;
- wrapped phase error;
- trace distance;
- Bloch vector;
- Bloch-vector error;
- purity;
- linear entropy;
- von Neumann entropy.

A population match alone is not sufficient for acceptance.

---

## 12. Required interpretation of v0.10-v0.11

The earlier center-frozen full-density error remains a valid benchmark of that
particular representation.

v0.12 does **not** erase or rewrite the old result.

Instead it adds the missing distinction:

```text
old comparison:
center-frozen / finite basis / dynamics
              versus
coordinate-dependent exact target

v0.12:
projection error        measured separately
projected dynamics      measured separately
original target error   measured separately
```

This correction must remain explicit in future releases.

---

## 13. PySCF status

The v0.12 representation-consistent release benchmark uses the analytic 2D LVC model.

The inherited PySCF backend remains regression-tested at the API and many-electron
overlap level, but v0.12 does not claim that the new spinor-complete or Born-Huang
benchmark path has been turned into a production molecular PySCF/AIMS implementation.

See `V12_PYSCF_BRIDGE.md`.
