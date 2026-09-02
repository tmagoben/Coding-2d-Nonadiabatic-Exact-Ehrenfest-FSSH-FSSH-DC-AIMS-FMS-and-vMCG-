# v0.25.3 algorithmic complexity

Let `N` be the Gaussian count, `s` the electronic-state count,
`P=2Ns+4N`, and `K=4N` the default candidate count.

- Nuclear overlap/projection matrices require `O(N^2)` analytic packet pairs.
- Full target SVD requires `O(N^3)` work and `O(N^2)` storage; electronic projection
  back-solves add `O(N^2 s)`.
- The inherited full TDVP metric stores `O(P^2)` values and uses a dense SVD with
  `O(P^3)` worst-case work.
- A candidate residual coupling contracts `P` tangent overlaps plus `N` Hamiltonian
  packet actions. The fixed-basis TDVP system is shared across a production candidate
  batch; candidate-specific analytic contractions scale as `O(K(Ps+Ns^2))`.
- Merge/prune admission can test `O(N^2)` possible removals in the worst case, each
  with a target SVD. The packet cap bounds this release path at `N<=8` by default.

This is a correctness-first dense implementation. Batched moment caches,
rank-one factor updates, and screened candidate pools are future optimizations and
are not claimed here.
