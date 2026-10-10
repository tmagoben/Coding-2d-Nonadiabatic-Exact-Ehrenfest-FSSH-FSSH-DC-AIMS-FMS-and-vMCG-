# v0.28.0 Moving Electronic Frame

v0.28.0 development extends the sealed v0.27.0 correlated-width Gaussian TDVP to a
coordinate-dependent electronic frame without changing the finite-dimensional packet
coefficient stored at each Gaussian centre.

## Variational section

The admitted wavefunction is

\[
\Psi(R)=\sum_I g_I(R)\,\Phi(R)\,W(R,q_I)\,c_I,
\]

where `c_I` is stored in the electronic frame at `q_I`. For the first v0.28.0
milestone the moving frame is an analytically trivializable pure gauge

\[
\Phi(R)=\Phi_{\rm ref}G(R),\qquad G(R)=e^{i\theta(R)K}U_0,
\]

with Hermitian `K`, real scalar quadratic `theta(R)`, and constant unitary `U0`.

## Connection and transporter

The frozen convention is

\[
D_a(R)=G^\dagger(R)\partial_aG(R),\qquad \nabla_a^{\rm cov}=\partial_a+D_a,
\]

and

\[
W(R,q)=G^\dagger(R)G(q).
\]

Under a constant right gauge `G -> G U`, `D_a -> U^dagger D_a U`,
`W -> U^dagger W U`, and `c_I -> U^dagger c_I`; the physical packet is unchanged.

## Exact trivialization

For every admitted packet,

\[
\Phi(R)W(R,q_I)c_I=\Phi_{\rm ref}G(q_I)c_I.
\]

Thus the moving-frame state maps exactly to the sealed v0.27.0 fixed-frame state with
`c_I^ref = G(q_I)c_I`. This identity is used as a high-precision oracle for
wavefunctions, TDVP velocities, implicit-midpoint endpoints, and basis lifecycle events.

## Moving coefficient velocity

Differentiating `c_ref = G(q)c` gives

\[
\dot c=G^\dagger(q)\dot c_{\rm ref}-\sum_a \dot q_a D_a(q)c.
\]

## Independent lattice oracle

A separate periodic finite-difference Hamiltonian is built from exact gauge links
`U_ij=G^dagger(R_i)G(R_j)`. Its kinetic off-diagonal blocks are `T_ij U_ij`.
The construction is checked against global unitary similarity of the fixed-frame
lattice Hamiltonian, as well as action and finite-time propagation covariance.

## Claim boundary

The current milestone validates only flat connections with an exact global
trivialization. It fails closed on nonzero curvature or missing trivialization.
It does not claim live molecular SOC, general curved electronic bundles, or full AIMS
branching semantics.
