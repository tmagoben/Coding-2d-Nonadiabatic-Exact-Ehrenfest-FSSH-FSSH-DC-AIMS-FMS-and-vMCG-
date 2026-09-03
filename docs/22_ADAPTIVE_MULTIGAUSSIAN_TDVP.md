# Adaptive multi-Gaussian TDVP (v0.25.2)

v0.25.2 gives every one-dimensional Gaussian a positive logarithmic width and a real
quadratic chirp. These variables join coefficient real/imaginary parts, centers, and
momenta in the same McLachlan metric and fully implicit midpoint solve.

Exact complex Gaussian moments through degree four support the complete analytic
metric. Validation includes independent grid reconstruction, exact harmonic thawed
equations, closed-form breathing, coherent reduction to v0.25.1, even/odd SOC,
signed reversal, packet permutation, constant electronic gauge covariance,
compatible duplicate-packet null directions, zero SOC, and second-order refinement.

Detailed release documents:

- `docs/releases/v0.25.2/V252_ADAPTIVE_MULTIGAUSSIAN_TDVP.md`
- `docs/releases/v0.25.2/V252_WIDTH_AND_SOLVER_POLICY.md`
- `docs/releases/v0.25.2/V252_PROGRAM_ARCHITECTURE.md`
- `docs/releases/v0.25.2/V252_ALGORITHM_COMPLEXITY.md`
- `docs/releases/v0.25.2/V252_VALIDATION.md`

The scope remains one coordinate, scalar packet widths/chirps, a fixed electronic
frame, constant mass, and a quadratic Hermitian matrix potential. Spawning/pruning,
multidimensional/full width matrices, and real molecular-SOC trajectories remain
closed.

