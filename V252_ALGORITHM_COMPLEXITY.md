# v0.25.2 algorithmic complexity

Let `N_g` be the packet count, `N_s` the electronic-state count, and

$$P=2N_gN_s+4N_g$$

the real adaptive parameter count. Let `K` be nonlinear function evaluations per
implicit step.

## One vector-field evaluation

- Chirped pair moments through degree four: `O(N_g^2)` cached complex scalars.
- Dense tangent metric: `O(P^2)` scalar contractions and `O(P^2)` storage.
- Dense RHS: worst-case `O(P N_g N_s^2)`.
- Combined overlap/Hamiltonian: `O(N_g^2 N_s^2)`.
- Full `P x P` SVD: `O(P^3)` time and `O(P^2)` workspace.

The dense reference step is therefore `O(K P^3)` in the asymptotic worst case. The
extra width/chirp pair increases `P` by `2N_g` relative to v0.25.1 and raises the
maximum moment degree from three to four, but does not change the cubic dense-SVD
asymptote.

## Receipt storage

Complete validation receipts store the midpoint metric, RHS, velocity, spectrum,
parameters, endpoint states, and solver diagnostics: `O(P^2)` per step and
`O(N P^2)` for `N` retained steps.

## Before multidimensional widths

Diagonal multidimensional widths require one width/chirp pair per coordinate and
full correlated Gaussians require symmetric matrix pairs. Before opening either, the
framework should introduce structured tangent indexing, pair-moment caching,
condition-aware block solvers, and explicit rotational covariance tests.

