# v0.22.0 physical analytic-SOC derivations

## 1. Operator decomposition and force convention

At nuclear geometry \(q\), v0.22.0 composes the electronic operator as

\[
H(q)=H_{\rm sf}(q)+H_{\rm SOC}(q),
\]

and exposes the physical Hamiltonian-derivative operators

\[
K_a(q)=\left(\partial_a H\right)_{\rm op}
=K_{a,{\rm sf}}(q)+K_{a,{\rm SOC}}(q).
\]

The derivative is a physical operator, not the derivative of a particular matrix
representation. For an electronic density matrix \(\rho\), the SOC contribution to
the nuclear force is

\[
F_{a,{\rm SOC}}=-\operatorname{Tr}\!\left[\rho K_{a,{\rm SOC}}\right].
\]

The validation fixes \(\rho\) and independently checks this expression against the
centered finite difference of \(-\operatorname{Tr}[\rho H_{\rm SOC}(q)]\).

## 2. Even-electron singlet–triplet model

The fixed spin-diabatic basis is

\[
\mathcal B_{ST}=
\left(|S\rangle,|T_{-1}\rangle,|T_0\rangle,|T_{+1}\rangle\right).
\]

The spin-free part contains displaced harmonic surfaces,

\[
H_{\rm sf}^{ST}(q)=
\operatorname{diag}\!\left(V_S(q),V_T(q),V_T(q),V_T(q)\right).
\]

With a complex linear \(\lambda(q)\) and a real linear \(\mu(q)\), the SOC block is

\[
H_{\rm SOC}^{ST}(q)=
\begin{pmatrix}
0 & \lambda & i\mu & \lambda^*\\
\lambda^* & 0 & 0 & 0\\
-i\mu & 0 & 0 & 0\\
\lambda & 0 & 0 & 0
\end{pmatrix}.
\]

The analytic \(K_{q,{\rm SOC}}^{ST}\) is obtained by replacing \(\lambda,\mu\) with
their exact derivatives. Setting all coupling gradients to zero therefore leaves a
nonzero constant \(H_{\rm SOC}\) and exactly zero \(K_{q,{\rm SOC}}\).

Writing time reversal as \(\Theta=J_{ST}\mathcal K\), where \(\mathcal K\) is complex
conjugation,

\[
J_{ST}=
\begin{pmatrix}
1&0&0&0\\
0&0&0&1\\
0&0&-1&0\\
0&1&0&0
\end{pmatrix},
\qquad J_{ST}J_{ST}^*=I.
\]

The coupling pattern \((\lambda,i\mu,\lambda^*)\) makes both \(H\) and \(K_q\)
time-reversal invariant:

\[
J_{ST}H^*J_{ST}^\dagger=H,
\qquad
J_{ST}K_q^*J_{ST}^\dagger=K_q.
\]

The physical observables are the complete singlet and triplet projectors,

\[
P_S=\operatorname{diag}(1,0,0,0),\qquad P_T=I-P_S.
\]

## 3. Odd-electron two-doublet model

The fixed spin-diabatic basis is

\[
\mathcal B_D=
\left(|D_1,+\tfrac12\rangle,|D_1,-\tfrac12\rangle,
|D_2,+\tfrac12\rangle,|D_2,-\tfrac12\rangle\right).
\]

The spin-free Hamiltonian is

\[
H_{\rm sf}^{D}(q)=
\begin{pmatrix}E_1(q)I_2&0\\0&E_2(q)I_2\end{pmatrix}.
\]

The inter-doublet SOC block has quaternionic form

\[
B(q)=
\begin{pmatrix}a(q)&b(q)\\-b(q)^*&a(q)^*\end{pmatrix},
\qquad
H_{\rm SOC}^{D}(q)=
\begin{pmatrix}0&B(q)\\B(q)^\dagger&0\end{pmatrix},
\]

where \(a\) and \(b\) are complex linear functions. The physical derivative has the
same structure with \(a,b\) replaced by \(\partial_q a,\partial_q b\).

For

\[
j_2=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\qquad
J_D=I_2\otimes j_2,
\]

one has

\[
J_DJ_D^*=-I,
\qquad J_DH^*J_D^\dagger=H.
\]

At zero magnetic field, \(\Theta^2=-I\) and time-reversal invariance imply Kramers
degeneracy. v0.22.0 verifies the twofold energy pairing across multiple geometries.
The physical root projectors are

\[
P_{D_1}=\operatorname{diag}(1,1,0,0),\qquad P_{D_2}=I-P_{D_1}.
\]

## 4. General complex gauge

For coefficient convention \(c'=G^\dagger c\), with \(G(q)\in U(4)\), physical
operators and projectors transform as

\[
H'=G^\dagger HG,\qquad K_a'=G^\dagger K_aG,
\qquad P'=G^\dagger PG.
\]

Because time reversal is antiunitary, its unitary part does not transform by ordinary
similarity. Instead,

\[
J'=G^\dagger JG^*.
\]

This relation is essential: using the untransformed \(J\) in a general complex frame
is deliberately tested as an invalid control. For a moving frame, the derivative
connection obeys

\[
D_a'=G^\dagger D_aG+G^\dagger\partial_aG,
\]

while \(K_a\) remains the representation of the physical derivative operator and
therefore has no connection term.

## 5. Independent exact-grid reference

The reference propagates a four-component spinor on a periodic uniform grid under

\[
i\partial_t\Psi(x,t)=
\left[-\frac{1}{2M}\partial_x^2 I_4+H(x)\right]\Psi(x,t).
\]

It does not call the Gaussian assembly or Gaussian propagation routines. A second-order
Strang step is

\[
\Psi(t+\Delta t)=
e^{-iH(x)\Delta t/2}
\mathcal F^{-1}e^{-ik^2\Delta t/(2M)}\mathcal F
e^{-iH(x)\Delta t/2}\Psi(t)+O(\Delta t^3).
\]

The matrix exponential at each grid point is evaluated from a Hermitian eigendecomposition.
Norm, energy, physical projector populations, timestep order, grid spacing, and box
size are checked independently of the moving-Gaussian trajectory.

## 6. Restart identity

SOC parameters, model-space order, electron parity, spin convention, derivative method,
and representation are part of the provider provenance fingerprint. A checkpoint may
resume only when that fingerprint and the numerical settings match. The SOC-active
tests exercise dense and sparse segment equivalence, sparse-edge hysteresis, moving
complex frames, parameter changes, and full-state SHA-256 corruption detection.
