# v0.29.3 Coordinate-Dependent Non-Abelian Geometry

## Status

v0.29.3 is a development checkpoint. It does not replace the certified v0.27.0 release. It preserves the v0.29.2 residual-driven curved basis controller while replacing the constant-connection / fixed-straight-section assumption with a coordinate-dependent non-Abelian connection field and explicit path/patch transport.

## Connection convention

The covariant derivative is

\[
D_\mu=\partial_\mu+A_\mu(R),
\qquad A_\mu^\dagger=-A_\mu.
\]

For a path \(\gamma:s\mapsto R(s)\), the Wilson transporter obeys

\[
\frac{dW}{ds}
=-\dot R^\mu(s)A_\mu(R(s))W(s),
\qquad W(0)=I.
\]

Therefore

\[
W(b\leftarrow a;\gamma)
=\mathcal P\exp\left[-\int_\gamma A_\mu(R)dR^\mu\right].
\]

The implementation evaluates generic paths by deterministic DOP853 ODE integration. Connection providers may expose an exact segment transporter for analytic fixtures and exact reductions.

## Explicit path families

A section is no longer specified only by its packet centre and endpoint. The path is part of the basis convention:

\[
\Phi_{Ia}(R;\gamma_I)=g_I(R)W(R\leftarrow q_I;\gamma_I)e_a.
\]

The checkpoint supplies straight and polyline/dogleg paths. For nonzero curvature, different paths are intentionally allowed to produce different transporters.

## Curvature and holonomy

With \(D_\mu=\partial_\mu+A_\mu\),

\[
F_{\mu\nu}
=\partial_\mu A_\nu-\partial_\nu A_\mu+[A_\mu,A_\nu].
\]

For a small positively oriented \(\mu\nu\) plaquette of area \(h^2\), the present Wilson convention gives

\[
W_\square=I-F_{\mu\nu}h^2+O(h^3).
\]

This is checked numerically against a finite-difference curvature oracle.

## Coordinate-dependent gauge covariance

For a unitary frame change \(G(R)\),

\[
A'_\mu
=G^\dagger A_\mu G+G^\dagger\partial_\mu G.
\]

The Wilson line must transform as

\[
W'(b\leftarrow a)
=G^\dagger(b)W(b\leftarrow a)G(a).
\]

A packet coefficient stored at its centre transforms as

\[
c'_I=G^\dagger(q_I)c_I,
\]

so its full section transforms pointwise as

\[
W'(R\leftarrow q_I)c'_I
=G^\dagger(R)W(R\leftarrow q_I)c_I.
\]

This is explicitly validated.

## Patch transitions

If patch \(\beta\) uses the frame obtained from patch \(\alpha\) by \(G_{\alpha\beta}(R)\), then

\[
A^{(\beta)}_\mu
=G_{\alpha\beta}^\dagger A^{(\alpha)}_\mu G_{\alpha\beta}
+G_{\alpha\beta}^\dagger\partial_\mu G_{\alpha\beta}.
\]

At a switching point \(R_*\), coefficient coordinates change by

\[
c^{(\beta)}(R_*)
=G_{\alpha\beta}^\dagger(R_*)c^{(\alpha)}(R_*).
\]

For the same geometric path split at \(R_*\),

\[
W^{(\beta)}(b\leftarrow R_*)
G_{\alpha\beta}^\dagger(R_*)
W^{(\alpha)}(R_*\leftarrow a)
=G_{\alpha\beta}^\dagger(b)
W^{(\alpha)}(b\leftarrow a).
\]

The checkpoint tests both the local connection relation and this same-path transport identity.

## Integration with v0.29.2 spawning

The v0.29.2 spawning score is unchanged:

\[
S_c=C-X^TM^+X,
\]

\[
y_c=q_c-X^TM^+g,
\]

\[
\Delta_c=y_c^TS_c^+y_c.
\]

The difference is that every current tangent and candidate coefficient block is now constructed from coordinate-dependent path-ordered sections. Full-SVD projection, merge/prune gates, one-event ordering, and dormant-newborn behavior are preserved.

## Explicit non-claims

This checkpoint does not yet claim:

- automatic atlas construction for arbitrary conical-intersection bundles;
- automatic selection/optimization of path families;
- patch-aware nonlinear midpoint propagation across changing patch membership;
- live molecular derivative-coupling or SOC providers inside this new geometry layer;
- molecular predictive accuracy;
- unrestricted multi-event AIMS branching;
- ML-controlled spawning.

The next method layer should bind this geometry to a live molecular electronic provider and validate derivative couplings/SOC against chemically meaningful reference trajectories.
