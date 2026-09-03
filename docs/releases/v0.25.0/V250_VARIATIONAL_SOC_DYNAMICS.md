# v0.25.0 restricted variational SOC dynamics

## 1. Released ansatz

The released variational manifold contains one canonical nuclear packet and one
complete electronic spinor,

$$
|\Psi(t)\rangle
= |g(q(t),p(t))\rangle
  \sum_{I=1}^{s} c_I(t)|\Phi_I(q(t))\rangle,
\qquad c^\dagger c=1.
$$

This is the Ehrenfest-like single-packet restriction of the time-dependent
variational principle (TDVP), not the full vMCG/AIMS multi-Gaussian manifold.

For a moving orthonormal electronic frame, define

$$
H_{IJ}(q)=\langle\Phi_I|\hat H_\mathrm{el}|\Phi_J\rangle,
\quad
K_{a,IJ}(q)=\langle\Phi_I|\partial_a\hat H_\mathrm{el}|\Phi_J\rangle,
\quad
D_{a,IJ}(q)=\langle\Phi_I|\partial_a\Phi_J\rangle.
$$

Here `K` is the physical Hamiltonian-operator derivative; it is not the derivative
of the matrix entries in an arbitrary moving frame. A convenient restricted action
is

$$
L = p^T\dot q + i c^\dagger\dot c
    + i\dot q^a c^\dagger D_a c
    - \left[\frac12p^T M^{-1}p+c^\dagger Hc\right].
$$

The corresponding restricted equations are

$$
\dot q=M^{-1}p,
\qquad
\dot p_a=-c^\dagger K_a c,
\qquad
i\dot c=(H-i\dot q^aD_a)c.
$$

The code never estimates `D` by taking an uncontrolled logarithm of an overlap.
Instead, finite-step frame motion is represented directly by certified
cross-geometry transport.

## 2. SVD-polar frame transport

For start and endpoint electronic frames,

$$
O_{01,IJ}=\langle\Phi_I(q_0)|\Phi_J(q_1)\rangle.
$$

A finite retained manifold makes `O_01` a contraction; it need not be unitary. Its
singular-value decomposition is

$$
O_{01}=U\Sigma V^\dagger.
$$

The right-to-left unitary polar factor is

$$
W_{01}=UV^\dagger,
\qquad W_{01}^\dagger W_{01}=I.
$$

Thus `W_01^dagger` transports a coefficient vector from the start frame to the
endpoint frame. `Sigma` is retained to test contraction consistency, minimum
retention, conditioning, and maximum principal angle. The raw contraction is never
silently substituted for unitary coefficient transport.

## 3. One signed symmetric step

Let `h` be positive or negative. From `(q_n,p_n,c_n)`:

1. Evaluate `M`, `H_n`, and `K_n` and form

   $$F_{n,a}=-c_n^\dagger K_{n,a}c_n.$$

2. Apply the first nuclear half kick,

   $$p_{n+1/2}=p_n+\frac h2F_n.$$

3. Drift the canonical coordinate,

   $$q_{n+1}=q_n+hM^{-1}p_{n+1/2}.$$

4. Evaluate the endpoint operator snapshot. The endpoint mass must equal the start
   mass within the frozen tolerance.
5. Compute `O_01 = U Sigma V^dagger`, certify it, and set `W_01 = U V^dagger`.
6. Apply the endpoint-Hamiltonian Strang/polar electronic step,

   $$
   c_{n+1}=e^{-iH_{n+1}h/2}
             W_{01}^\dagger
             e^{-iH_nh/2}c_n.
   $$

7. Form `F_{n+1,a}=-c_{n+1}^dagger K_{n+1,a}c_{n+1}` and finish the kick,

   $$p_{n+1}=p_{n+1/2}+\frac h2F_{n+1}.$$

Hermitian exponentials and unitary polar transport preserve the electronic norm up
to floating-point roundoff. Accepting a signed `h` makes the adjoint calculation an
ordinary application of the same step with reversed order and negative time.

## 4. Why this is not full TDVP

The full Gaussian ansatz would contain several moving packets, complex amplitudes,
possibly moving width matrices, and a nonorthogonal metric. Its Euler--Lagrange
equations have the generic coupled form

$$
i\,\mathcal M(z,z^*)\dot z=\frac{\partial E}{\partial z^*},
$$

with constraints and gauge/null directions. These coordinates are generally not a
separable canonical `(q,p)` system. Plain velocity Verlet would therefore lose its
usual symplectic justification. A future implementation should solve an implicit
midpoint or discrete variational residual, with metric regularization and constraint
handling, before opening the full-TDVP claim.

## 5. Provider boundary

The trajectory provider must return an arbitrary-geometry
`ElectronicOperatorSnapshotV21` containing full `H`, physical `K`, anti-Hermitian
`D`, a positive mass matrix, and state vectors, plus a cross-snapshot overlap method.
The analytic v0.22 models satisfy this contract. The real v0.24.1/v0.24.2 PySCF SOC
objects remain static or differential evidence and are intentionally rejected by
this trajectory entry point.
