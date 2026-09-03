# v0.23.3 finite-manifold overlap and transport

For orthonormal retained states at two geometries, v0.23.3 defines

$$
O_{LR,ij}=\langle\Phi_i(q_L)|\Phi_j(q_R)\rangle.
$$

`O_LR` is a physical contraction, not generally a unitary matrix. Singular
values below one measure amplitude lost to omitted roots. Singular values above
one (outside tolerance) violate the finite-orthonormal-manifold contract.

## Transport

If `O_LR = U Sigma V†`, the right-to-left coefficient transport is the unitary
polar factor

$$
W_{L\leftarrow R}=UV^\dagger.
$$

The raw overlap remains available for diagnostics and amplitude information;
`W` is used for basis transport. Under independent endpoint gauges `G_L,G_R`,

$$
O' = G_L^\dagger O G_R,\qquad
W' = G_L^\dagger W G_R.
$$

## Two distinct decisions

Physical consistency requires finite square data, contraction singular values,
unitary polar transport, and a positive-semidefinite polar factor. Trajectory
readiness additionally requires the configured minimum singular value, maximum
condition number, and maximum principal angle. A physically meaningful overlap
can therefore be rejected as too poorly retained for propagation.

The default release policy requires a minimum singular value of 0.5, condition
number at most `1e6`, and principal angle at most `pi/3`. Low-level Procrustes
consumers use a permissive nonzero-retention policy so old well-conditioned use
cases remain compatible while rank loss fails closed.

Directed overlap pairs must satisfy raw adjoint reciprocity and unitary-transport
adjoint reciprocity independently.
