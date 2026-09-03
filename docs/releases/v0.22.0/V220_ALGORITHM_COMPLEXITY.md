# v0.22.0 algorithmic complexity

Let

- $s$ be the electronic model-space dimension ($s=4$ in both v0.22.0 models),
- $d$ the nuclear dimension ($d=1$ in the analytic references),
- $N_x$ the exact-grid point count,
- $N_t$ the number of time steps,
- $N$ the number of Gaussian basis functions, and
- $E$ the number of active Gaussian-pair edges.

## Analytic providers

Constructing $H$, $K$, and the physical projectors requires
$O(ds^2)$ time and $O(ds^2)$ storage. The concrete four-state formulas have fixed
cost, but the dimension-dependent notation makes their relation to future backends
explicit.

## Contract and symmetry audits

Hermiticity, composition, time-reversal, and projector residuals require matrix
products and therefore at most $O(ds^3)$ time with $O(ds^2)$ working memory.
Centered differential and force audits multiply that cost by a fixed number of
neighboring geometries.

Kramers certification diagonalizes the Hamiltonian at each sampled geometry, costing
$O(n_qs^3)$ for $n_q$ geometries.

## Exact-grid reference

Potential construction stores $N_x$ dense $s\times s$ matrices:

$$
\text{memory}=O(N_xs^2).
$$

The pointwise Hermitian eigendecompositions used to prepare the potential half-step
cost $O(N_xs^3)$. Each propagation step then costs

$$
O(N_xs^2+sN_x\log N_x),
$$

from the pointwise dense matrix actions and $s$ FFTs. A complete trajectory costs

$$
O\!\left(N_xs^3+N_t(N_xs^2+sN_x\log N_x)\right).
$$

Stored wavefunctions cost $O(N_{\mathrm{save}}N_xs)$; the implementation stores only
requested sampling points.

## Gaussian propagation

SOC does not change the asymptotic Gaussian engine established in earlier releases.
It changes the electronic values inside the same blocks:

- dense block assembly/propagation retains its existing $N$ and $s$ scaling;
- sparse assembly retains $O(Es^2)$ block storage;
- provider evaluation remains attached to Gaussian centers and active pair operations;
- checkpoint data grow only with the already-present electronic block dimension and
  sparse graph state.

The v0.22.0 exact-grid implementation is a validation oracle, not a production
high-dimensional solver. Its exponential dependence on nuclear dimension through a
direct-product grid is intentionally avoided in the Gaussian engine.
