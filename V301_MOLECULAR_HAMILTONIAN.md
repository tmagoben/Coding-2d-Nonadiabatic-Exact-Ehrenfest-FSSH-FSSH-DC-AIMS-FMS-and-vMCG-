# v0.30.1 development — molecular Hamiltonian action

The molecular-frame bridge now supplies both covariant nuclear kinetic energy and
the total electronic Hamiltonian to the recovered curved Gaussian TDVP residual.
This is the first implementation slice after the saved v0.30.0 checkpoint. The
validated package release remains v0.27.0.

## Operator and boundary conditions

For rectilinear coordinates with constant diagonal nuclear mass, the declared
projected model is

$$
\hat H=-\sum_a\frac{1}{2m_a}(\partial_a+D_a)^2+H_{\mathrm{el}}(R),
\qquad D_a^\dagger=-D_a.
$$

The electronic operator includes SOC when the provider enables it. The implementation
does not infer the additional Born–Huang scalar operator associated with electronic
states outside a retained manifold. Consequently, “full Hamiltonian” here means
both terms of the model above, not a certification of a complete molecular projection.

Each axis is a finite, strictly increasing uniform grid. Its entries are interior
nodes. The wavefunction is zero at ghost nodes one step outside either end; there
is no periodic wrap. With $W_{j\leftarrow i}$ transporting coefficients from node
$i$ to its neighbour $j$, the kinetic blocks are

$$
T_{ii}=\sum_a\frac{1}{m_a h_a^2}I,\qquad
T_{ij}=-\frac{W_{j\leftarrow i}^\dagger}{2m_a h_a^2},\qquad
T_{ji}=T_{ij}^\dagger.
$$

The grid inner product uses uniform cell-volume weights. Reverse blocks are
assembled as adjoints, giving a Hermitian discrete Hamiltonian. Wilson links use the
NAC ODE, with every edge independently checked against certified reciprocal
overlap-polar transport. A failed local-length or transport-consistency check rejects
construction. Point samples, audit receipts, grid, boundary, mass, and sparse matrix
entries contribute to the operator fingerprint.

Under a coordinate-dependent electronic gauge, the link transforms as
$W'_{j\leftarrow i}=G_j^\dagger W_{j\leftarrow i}G_i$. The whole discrete operator
therefore transforms by block-diagonal unitary similarity, without assuming a
single global electronic frame. Numerical Wilson integration sets the accuracy of
this covariance.

For $N$ tensor-grid points, $d$ coordinates, and $s$ electronic states, the sparse
operator stores $O(N d s^2)$ entries. Assembly evaluates $O(Nd)$ transport segments
in addition to $N$ electronic snapshots; provider and Wilson-integration costs are
additional. With $P$ real packet parameters, the inherited finite-difference TDVP
stores an $Ns$ by $P$ complex Jacobian and uses a dense $P$ by $P$ metric solve.
The grid size is the product of the axis sizes, so this is a small-system validation
operator, not a scalable substitute for high-dimensional molecular matrix elements.

## TDVP integration

`build_molecular_tdvp_snapshot_v301` evaluates the path-dependent packet sections,
applies the assembled kinetic-plus-electronic operator, then calls the recovered
v0.29.3 full-parameter tangent construction and SVD solve. Its result contains the
wavefunction, $H\Psi$, metric/velocity/residual receipt, and Hamiltonian fingerprint.
The provider's grid-point fingerprints are rechecked before constructing the tangent
basis, so changed sampled operators cannot silently reuse a stale grid Hamiltonian.

`MolecularGridHamiltonianV301.hamiltonian_action` also implements the existing
curved-energy/basis-controller callback. It rejects different quadrature coordinates
or weights. The current tests validate the TDVP snapshot integration; they do not
certify an end-to-end adaptive topology trajectory with this new operator.

The underlying provider must be stationary over the evaluation. Mass constancy is
checked at its reference point and all grid points. No finite sample can prove
constancy between samples. Variable masses and off-diagonal mass tensors are outside
this discretization and are rejected when detected.

## Reproduction

From the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -m pytest -q tests/test_v0301_molecular_hamiltonian.py

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python examples/158_recompute_v0301_molecular_hamiltonian.py
```

The example writes `results/v0301_molecular_hamiltonian_evidence.json`. It uses the
analytic singlet/triplet SOC model. Its exponential propagation is an independent
finite-grid reference, not a new Gaussian time integrator. The recorded TDSE defect
is an approximation residual, not an assertion that the Gaussian ansatz is exact.

Tests compare the zero-connection operator with independent 1D and 2D Dirichlet
stencils, check Hermiticity and kinetic positivity, exercise moving-gauge action and
TDVP-residual covariance, and verify second-order kinetic refinement against an
analytic Gaussian derivative. They also compare the free-packet TDVP velocity with
the analytic centre, momentum, log-width, and chirp equations.

Negative controls cover bad grids, nonfinite wavefunctions, variable and off-diagonal
masses, inconsistent NAC/overlap transport, excessive link length, bridge mismatch,
changed provider samples, and mismatched quadrature. Coincident packet centres and
grid nodes have identity transport without displacing the point.

## Recovery and remaining work

The [recovery record](docs/development/recovered-v0300/README.md) documents the
verified source archive and preserves its historical notes. This branch restores
the curved engine but does not promote its previous claims to a current release.
The restored package-level v0.29.4 record schema and the existing top-level
v0.29.4/v0.29.5 dataset schema remain separate.

The [v0.30.2 development slice](V302_MOLECULAR_MIDPOINT.md) adds converged
implicit-midpoint Gaussian trajectories using this action, with replayable stepwise
transport and Hamiltonian receipts. The two dataset contracts still need an explicit audited adapter, real physical
symmetry re-evaluations, and trajectory-level split/leakage checks. Live ab-initio SOC
trajectories, changing atlas membership, molecular accuracy, and learned ranking
remain later milestones.
