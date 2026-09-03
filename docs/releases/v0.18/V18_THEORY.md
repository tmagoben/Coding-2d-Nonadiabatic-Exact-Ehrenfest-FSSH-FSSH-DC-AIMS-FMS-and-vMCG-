# v0.18 Theory: Convergence-Complete Sparse Gaussian Dynamics

v0.18 is the **dynamics and convergence-completeness release**.

The sparse S/H/T controller introduced in v0.17 is retained, but v0.18 changes the
question from

> Is the sparse reduced electronic result acceptable?

to

> Does the complete Gaussian wavefunction converge under independently controlled
> timestep, basis, sparsity, and adaptive-growth coordinates?

The release therefore adds four major layers:

1. full spinor Gaussian wavefunction reconstruction on the exact 2D grid;
2. phase-aligned wavefunction, nuclear-density, centroid, and covariance errors;
3. physically normalized adaptive control times for timestep studies;
4. sampled omitted-edge audits plus only initial/final dense matrix sentinels.

---

## 1. Spinor-complete Gaussian reconstruction

For $N$ nuclear Gaussians and $s$ electronic components,

$$
\boxed{
\Psi_G(R,t)
=
\sum_{i=1}^N
\sum_{a=1}^s
C_{ia}(t)
g_i(R,t)
|a\rangle .
}
$$

For the 2-state analytic LVC benchmark,

$$
s=2.
$$

The reconstruction uses the actual propagated coefficients and current moving Gaussian
parameters.

This allows v0.18 to compare the **complete approximate state**, not only its reduced
electronic density matrix.

---

## 2. Global-phase alignment

A global quantum phase is physically irrelevant.

For normalized reference and candidate wavefunctions, define

$$
z=
\langle\Psi_{ref}|\Psi_G\rangle .
$$

The phase-aligned Gaussian state is

$$
\boxed{
\widetilde\Psi_G
=
e^{-i\arg z}
\Psi_G.
}
$$

The fidelity is

$$
\boxed{
F
=
|\langle\Psi_{ref}|\Psi_G\rangle|^2.
}
$$

The phase-aligned grid error is

$$
\boxed{
\epsilon_\Psi
=
\|\Psi_{ref}-\widetilde\Psi_G\|_2.
}
$$

For normalized states,

$$
\epsilon_\Psi^2
=
2-2\sqrt F.
$$

So fidelity and phase-aligned L2 error are mathematically consistent but are both
reported because they emphasize different regimes.

---

## 3. Nuclear-density observables

The spin-summed nuclear density is

$$
\boxed{
n(R)
=
\sum_a |\Psi_a(R)|^2.
}
$$

v0.18 reports

$$
\epsilon_n
=
\left[
\int
|n_G(R)-n_{ref}(R)|^2dR
\right]^{1/2},
$$

and total-variation distance

$$
\boxed{
D_{TV}
=
\frac12
\int
|n_G(R)-n_{ref}(R)|dR.
}
$$

The first two spatial moments are

$$
\mu
=
\int R n(R)dR,
$$

and

$$
\Sigma
=
\int
(R-\mu)(R-\mu)^T
n(R)dR.
$$

The benchmark therefore records centroid and covariance errors in addition to
electronic populations.

---

## 4. Representation error versus propagation error

The initial exact wavepacket is not represented perfectly by the finite 10-Gaussian
initial basis.

Its initial projection fidelity is

$$
\boxed{
F_{proj}(0)
=
0.882251454460.
}
$$

Let

$$
|\Psi_{target}(0)\rangle
$$

be the original exact initial packet and

$$
|\Psi_{proj}(0)\rangle
$$

its finite-Gaussian projection.

Both exact reference states are propagated under the same Hamiltonian.

Therefore exact unitarity requires

$$
\boxed{
|\langle
\Psi_{target}(t)
|
\Psi_{proj}(t)
\rangle|^2
=
|\langle
\Psi_{target}(0)
|
\Psi_{proj}(0)
\rangle|^2 .
}
$$

Measured exact-grid fidelity drift is only

$$
\boxed{
1.665e-15.
}
$$

This is important: comparing the Gaussian trajectory to the **projected exact**
trajectory isolates propagation/basis-dynamics error.

Comparing it directly to the original exact target mixes propagation error with initial
representation error.

---

## 5. Canonical v0.18 result

Against the projected exact reference at $t=0.6$:

$$
\boxed{
F_G
=
0.982566093412
}
$$

and

$$
\boxed{
\epsilon_\Psi
=
0.132327478361.
}
$$

The nuclear-density L2 error is

$$
\boxed{
\epsilon_n
=
0.052341235444.
}
$$

The centroid error is

$$
0.0012693911,
$$

and the covariance Frobenius error is

$$
0.012697602.
$$

The reduced electronic density error remains much smaller:

$$
\boxed{
\|\rho_G-\rho_{exact,proj}\|_F
=
0.00010573932.
}
$$

This difference between reduced-state and full-wavefunction error is scientifically
useful: a reduced electronic observable can converge earlier than the full nuclear-
electronic wavefunction.

---

## 6. Basis-completeness ladder

v0.18 explicitly varies the maximum adaptive basis size:

| $N_{max}$ | Final $N$ | Phase-aligned L2 error | Fidelity |
|---:|---:|---:|---:|
| 10 | 10 | 0.191190020 | 0.963780418 |
| 11 | 11 | 0.177501964 | 0.968741225 |
| 12 | 12 | 0.144924302 | 0.979107229 |
| 13 | 13 | 0.132327478 | 0.982566093 |

The wavefunction error decreases strictly from

$$
0.191190020
$$

to

$$
0.132327478,
$$

a relative improvement of

$$
\boxed{
30.79\%.
}
$$

This is the clearest v0.18 evidence that the extra Gaussian basis functions are
improving the actual state rather than merely adding parameters.

---

## 7. Timestep self-convergence

All adaptation intervals are expressed in physical time before being converted to step
counts.

Thus the enrichment events occur at the same physical times for

$$
\Delta t
=
0.010,\ 0.005,\ 0.0025.
$$

Define successive solution differences

$$
D_h
=
\|
\Psi_h-\Psi_{h/2}
\|_2,
$$

after global-phase alignment.

Measured:

$$
D_{0.01}
=
0.000610168531085,
$$

$$
D_{0.005}
=
0.000152815765496.
$$

The Richardson-style observed order is

$$
\boxed{
p_{obs}
=
\frac{
\ln(D_h/D_{h/2})
}{
\ln2
}
=
1.997414387.
}
$$

This is essentially second order.

The important point is that this order is inferred from **successive Gaussian
wavefunctions**, not from the nonzero Gaussian-model error relative to the exact grid
solution.

---

## 8. Sparse-edge-budget convergence

The local omitted-edge score budget is independently tightened:

| $B_{local}$ | Wavefunction L2 error | Average graph sparsity |
|---:|---:|---:|
| 0.030 | 0.145733746 | 0.040389 |
| 0.010 | 0.132327478 | 0.013860 |
| 0.000 | 0.132124297 | 0.000000 |

The result approaches the fully retained local graph as

$$
B_{{local}}\rightarrow0.
$$

This gives sparsification a genuine convergence coordinate rather than treating a
single threshold as universal.

---

## 9. Adaptive-growth threshold sensitivity

The TDSE-defect growth trigger is also swept:

| Enrichment threshold | Final basis | Enrichment steps | Wavefunction L2 error |
|---:|---:|---|---:|
| 0.050 | 10 | [] | 0.191190020 |
| 0.035 | 12 | [70, 120] | 0.180089293 |
| 0.030 | 13 | [10, 20, 70] | 0.138203917 |
| 0.025 | 13 | [10, 20, 30] | 0.132327478 |
| 0.015 | 13 | [10, 20, 30] | 0.132327478 |

The error improves as the controller is permitted to enrich earlier and then reaches a
plateau once the same 13-Gaussian basis-growth history is selected.

This is preferable to interpreting one spawn threshold as a physically privileged
constant.

---

## 10. Sampled sparse audits

v0.17 performed complete dense S/H/Snuc audits throughout the trajectory.

v0.18 keeps complete dense audits only at:

```text
t = 0
t = final
```

Normal checkpoints use exact S/H/T scoring on omitted-edge samples.

The sample combines:

1. priority omitted edges near the geometric search boundary;
2. deterministic pseudo-random omitted edges.

A sampled edge is considered a violation if its exact S/H/T score exceeds the current
graph enter threshold.

This is a diagnostic error-control layer, **not a probabilistic proof that no omitted
edge is important**.

The initial/final dense sentinels remain the authoritative release check.

---

## 11. Candidate-grid batching

For $K$ adaptive candidates and $G$ grid points, earlier ranking materialized

$$
O(KG)
$$

candidate grid values simultaneously.

v0.18 processes batches of size $B$:

$$
\boxed{{
\text{{peak candidate-grid storage}}
=
O(BG),
\qquad
B\ll K.
}}
$$

The total contraction work remains approximately $O(KG)$.

The improvement is therefore a **peak-memory reduction**, not a claim of lower
arithmetic complexity.

Measured maximum dense candidate-grid size:

```text
{comp["candidate_max_dense_grid_elements"]} complex elements
```

Measured batched peak:

```text
{comp["candidate_peak_grid_elements"]} complex elements
```

giving

$$
\boxed{{
{100*comp["candidate_peak_memory_reduction_fraction"]:.2f}\%
}}
$$

peak candidate-grid storage reduction.

---

## 12. Conditioning tradeoff

The improved 13-Gaussian wavefunction reaches a maximum overlap-metric condition number

$$
\boxed{{
\kappa(S_{{nuc}})
=
{c["maximum_condition_number"]:.4f}.
}}
$$

This is higher than the smaller v0.17 basis, but still below the v0.18 release ceiling

$$
10^4
$$

and far below the runtime candidate rejection limit

$$
10^5.
$$

The generalized norm drift remains only

$$
{c["maximum_norm_drift"]:.3e}.
$$

v0.18 therefore keeps the more complete basis and documents the conditioning cost
rather than pruning it away purely to improve the condition-number headline.

---

## 13. Scientific label

The appropriate description is:

> **convergence-controlled, sampled-audited sparse spinor-complete Gaussian
> nonadiabatic dynamics with full-wavefunction validation on an analytic 2D
> conical-intersection benchmark.**

It is still not:

- a production AIMS implementation;
- a full molecular direct-dynamics engine;
- a spin-orbit-coupled release;
- a proof that one sparse score is universal across Hamiltonians.
