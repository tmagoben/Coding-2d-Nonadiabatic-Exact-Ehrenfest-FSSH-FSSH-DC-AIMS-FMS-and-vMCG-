# v0.25.1 algorithmic complexity

Let `N_g` be the Gaussian count, `N_s` the electronic-state count, and

$$P=2N_gN_s+2N_g$$

the number of real variational parameters. Let `K` be the number of nonlinear
function evaluations in one implicit step.

## One TDVP vector-field evaluation

- Pairwise Gaussian moments through degree three: `O(N_g^2)` time and storage if
  cached for the evaluation.
- Dense tangent metric construction: `O(P^2)` contracted scalar entries and
  `O(P^2)` storage.
- RHS construction: worst-case `O(P N_g N_s^2)` for dense electronic matrices.
- Full dense SVD of the `P x P` metric: `O(P^3)` time and `O(P^2)` workspace.
- Analytic generalized overlap/Hamiltonian construction: `O(N_g^2 N_s^2)` time and
  storage for the explicit combined matrices.

For the present dense reference implementation, SVD dominates as the variational
space grows. One implicit midpoint step costs approximately `K` vector-field
evaluations plus the start and final receipt rebuilds, hence `O(K P^3)` in the dense
worst case.

## Storage

The validation trajectory retains complete receipts, including midpoint metric,
RHS, velocity, and singular spectrum. This is `O(P^2)` per step and `O(N P^2)` for
`N` retained steps. A later production runner may stream or checkpoint compressed
receipts, but v0.25.1 deliberately favors auditability.

## Scaling direction

Adaptive widths add at least two real parameters per one-dimensional packet and
more tangent polynomial degree. Multidimensional correlated widths add matrix-valued
parameters. These should not be enabled until rank-aware block structure, cached
pair moments, iterative nonlinear linear algebra, and controlled basis management
are separately validated.

