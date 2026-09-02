# v0.23.2 finite-manifold overlap contract

Let P(q) and P(q') project onto the selected electronic roots at two geometries.
The stored cross-geometry matrix is the finite-manifold restriction

$$
U(q,q')_{ij}=\langle\Phi_i(q)|\Phi_j(q')\rangle.
$$

It is generally a **contraction**, not a unitary matrix. Selected states at q'
can contain amplitude in roots outside P(q), so requiring `U^dagger U=I` is
physically incorrect for a truncated state manifold.

The v0.23.2 contract requires:

1. exact-shape, finite complex matrices;
2. self overlaps equal identity within tolerance;
3. reciprocity `U(q',q)=U(q,q')^dagger` within tolerance;
4. every cross-overlap singular value no greater than one within tolerance.

Small singular values and a nonzero isometry defect are reported as diagnostics,
not automatically rejected. Singular values above one indicate an impossible
norm expansion and are rejected. Applications can impose a separate minimum
retention threshold appropriate to their selected manifold.

In the canonical real H3+ calculation, the two-geometry overlap has minimum
singular value about 0.999924 and isometry defect about 1.98e-4. This is a valid,
high-retention contraction.
