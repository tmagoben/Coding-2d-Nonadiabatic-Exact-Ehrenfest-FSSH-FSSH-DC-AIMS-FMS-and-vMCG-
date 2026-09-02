# Controlled basis adaptation

v0.25.3 adds a conservative variable-basis controller to the v0.25.2 adaptive
Gaussian TDVP. The full derivation, gates, and architecture are in:

- [`../V253_CONTROLLED_BASIS_ADAPTATION.md`](../V253_CONTROLLED_BASIS_ADAPTATION.md)
- [`../V253_LIFECYCLE_POLICY.md`](../V253_LIFECYCLE_POLICY.md)
- [`../V253_PROGRAM_ARCHITECTURE.md`](../V253_PROGRAM_ARCHITECTURE.md)
- [`../V253_VALIDATION.md`](../V253_VALIDATION.md)

The most important restriction is semantic: “spawning” means one residual-selected,
SVD-projected packet at a checkpoint in the certified one-dimensional quadratic
model. It does not imply general multidimensional AIMS branching.
