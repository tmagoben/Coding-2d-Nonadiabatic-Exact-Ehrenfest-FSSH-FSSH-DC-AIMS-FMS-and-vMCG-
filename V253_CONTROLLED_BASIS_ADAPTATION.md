# v0.25.3 controlled basis adaptation: mathematical specification

## 1. Inherited variational state

The wavefunction remains

$$
|\Psi\rangle=\sum_{I=1}^{N_g}g_I(x;q_I,p_I,\alpha_I,\beta_I)
\sum_{a=1}^{N_s}C_{Ia}|a\rangle,
$$

with `eta_I=log(alpha_I)` and the exact v0.25.2 McLachlan/implicit-midpoint
propagator. Topology changes occur only at a completed fixed-basis step endpoint.

## 2. Residual-driven candidate score

Let

$$R=\dot\Psi+i\hat H\Psi$$

use the accepted v0.25.2 minimum-norm TDVP velocity. For a normalized trial packet
`g_c`, define the current nuclear Gram matrix and overlap vector

$$S_{IJ}=\langle g_I|g_J\rangle,\qquad b_I=\langle g_I|g_c\rangle.$$

The projection coefficients, novelty, and normalized orthogonal direction are

$$a=S^{-1}b,\qquad \nu=1-b^\dagger a,$$

$$|\widetilde g_c\rangle=
\frac{|g_c\rangle-\sum_Ia_I|g_I\rangle}{\sqrt\nu}.$$

For electronic component `a`, production evaluates the exact complex quantity

$$r_{c,a}=\langle g_c,a|R\rangle$$

using Gaussian moments through degree four. With
`r_{I,a}=<g_I,a|R>`, the released capture score is

$$
\sigma_c=\frac{\|r_c-a^\dagger r_{\rm current}\|_2}{\sqrt\nu}.
$$

A candidate is admitted only if `nu` exceeds the novelty gate, the enlarged nuclear
overlap is full rank and conditioned, the width/chirp domain is valid, the packet
cap is not reached, and `sigma_c` exceeds its threshold.

## 3. Candidate pool and deterministic selection

Each parent produces four candidates. With
`sigma_q=1/sqrt(2 alpha_I)` and `sigma_p=sqrt(alpha_I/2)`, production displaces
`q_I` by `+-2 sigma_q` or `p_I` by `+-2 sigma_p`, preserving the parent's width and
chirp. The highest admitted score is selected; exact ties use a canonical physical
geometry key rather than input order.

## 4. Full-SVD variational projection

For target packets `h_A`, define

$$S^T_{AB}=\langle h_A|h_B\rangle,\qquad
B_{AI}=\langle h_A|g_I\rangle.$$

For every electronic component, the least-squares coefficients solve

$$S^T D=BC.$$

Production uses a full SVD, requires full retained target rank and bounded condition
number, and records the linear residual. The exact squared projection error is

$$
\epsilon_{\rm proj}^2=
\langle\Psi|\Psi\rangle-2\operatorname{Re}\operatorname{Tr}(D^\dagger BC)
+\operatorname{Tr}(D^\dagger S^T D).
$$

The projected state is normalized. Relative loss, normalized fidelity, and the
energy discontinuity are all independently recomputed by the receipt.

## 5. Newborn activation

Exact enlargement contains the old basis, so a unique linearly independent
projection gives the newborn zero amplitude. Coefficient real/imaginary tangents
remain meaningful, but its shape tangents are proportional to the zero electronic
vector. For a dormant packet, v0.25.3 solves the reduced system

$$G_{AA}\dot\theta_A=b_A,$$

where `A` contains every coefficient coordinate and only the shape coordinates of
packets whose row population

$$P_I=\sum_a|C_{Ia}|^2$$

is at least `1e-6`. Frozen coordinates have exactly zero velocity. The activation
mask is fixed over each implicit midpoint solve and reconsidered at the next step.

## 6. Pruning and merging

An age-eligible low-population packet may be pruned only when projection into the
remaining basis satisfies the prune loss and energy gates. A pair may merge only
when `|S_IJ|>=0.999`; both possible survivors are projected and the lower-loss
choice is retained. “Merge” therefore means conservative merge-to-survivor, not an
unvalidated averaged packet geometry.

## 7. Lifecycle order

At most one event is permitted per checkpoint:

1. merge an admitted redundant pair;
2. otherwise prune an admitted low-population packet;
3. otherwise spawn the highest admitted residual candidate;
4. otherwise record a no-op.

This isolates every discontinuous basis decision in one projection receipt.
