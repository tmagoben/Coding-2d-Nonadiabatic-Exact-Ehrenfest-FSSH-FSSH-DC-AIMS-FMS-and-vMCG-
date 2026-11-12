# v0.30.3 development — analytic SOC trajectory campaign

The v0.30.2 solver had a scalar time-refinement campaign and a single short SOC
step. This iteration adds a complete singlet–triplet trajectory in both a fixed
spin-diabatic frame and a coordinate-dependent mixed electronic frame. It uses the
existing implicit-midpoint solver and audited covariant grid Hamiltonian.

## Model and independent comparisons

One correlated-width Gaussian starts as a pure singlet, with centre 0.1, momentum
0.3, width 1.1 and chirp 0.08. The analytic v0.22 singlet–triplet provider uses mass
2 and SOC scale 0.2; its remaining parameters retain their documented defaults.
These deliberately accelerated parameters test the algorithm and are not a fit to
a molecule. The tensor grid has 13 interior points from -3 to 3 and zero Dirichlet
ghost values one spacing beyond the endpoints.

The final time is 0.4 atomic units, ten times the earlier SOC smoke-test step.
Both electronic frames run 2, 4 and 8 steps: dt = 0.2, 0.1 and 0.05. The moving
frame uses a reproducible random unitary (seed 303) and coordinate-dependent phases
with gradients (0.3, -0.2, 0.1, 0.25). Initial coefficients are transformed at the
packet centre. Every final wavefunction is transformed back at each grid point
before comparing physical singlet/triplet populations.

Three comparisons answer different questions:

- Successive Gaussian trajectories measure temporal self-convergence on the same
  grid and manifold. Both frames must show order between 1.8 and 2.2.
- A sparse matrix exponential evolves the initial sampled wavefunction with the
  same complete grid Hamiltonian, without the Gaussian ansatz or midpoint solver.
  It measures the combined Gaussian projection and temporal error against that
  finite-grid model. It is not a continuum reference.
- Physical wavefunctions in the two electronic frames measure numerical gauge
  agreement. Midpoint in nonlinear parameter charts is not exactly covariant at
  finite step size. Gauge error must fall by more than a factor of eight across
  the fourfold timestep refinement.

The whole coarse moving-frame trajectory is serialized and exactly replayed.
An independent run with `soc_enabled=False` and zero SOC scale uses eight steps
and must retain a pure singlet in both Gaussian and grid propagation.

## Acceptance and evidence

`evaluate_soc_gates_v303` recomputes fourteen gates from the measured values.
Stored `passed` or `gates` flags cannot override numerical failures. Fine-step
limits include norm drift below 2e-7, physical gauge wavefunction error below
2e-5, triplet gauge error below 2e-6, triplet error versus the grid reference below
5e-4, and weighted wavefunction error versus that reference below 1e-2. Accepted
nonlinear and tangent-equation residuals are also required. The SOC-enabled case
must transfer more than 1e-3 population; the disabled case must remain below 1e-12.

Norm and energy drift are maxima over the accepted step endpoints, not just the
final endpoint. Energy drift is reported without claiming a universal conservation
threshold for this projected finite-grid model. Both frames' errors, convergence
orders, physical populations, operator fingerprints, runtime versions, and coarse
trajectory receipts are saved in `results/v0303_soc_trajectory_evidence.json`.

Recompute from the repository root:

```bash
python -m pip install -e ".[dev]"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python examples/160_recompute_v0303_soc_trajectory.py
python -m pytest -q tests/test_v0303_soc_trajectory.py
```

The example exits unsuccessfully if any gate fails. CI recomputes the complete
campaign in a dedicated job and uploads a fresh `ci-v0303-soc-trajectory.json`.
The CLI defaults to two worker processes (`--workers 1` selects sequential
execution). Each process owns its provider and random seed; the moving-frame
trajectory and its replay run in the same process. Results are collected in a
fixed order. The fast tests check the stored evidence,
perturb individual measurements to reject false success, and demonstrate that
mixed-gauge component magnitudes invent triplet population unless the physical
spin projectors are transformed consistently. The existing cumulative suites still
run on every PR; historical runtime-identity distinctions remain as documented in
v0.30.2.

## Recorded local result

All 14 campaign gates pass in the committed evidence. The 19 focused campaign
checks pass, as do the 55 targeted Hamiltonian/midpoint/physical-population checks.
The wheel builds and imports in isolation.

| Frame | Time order | Fine norm drift | Fine triplet population |
| --- | ---: | ---: | ---: |
| Fixed | 1.992519 | 1.668e-09 | 0.009467483 |
| Moving | 1.992676 | 1.634e-09 | 0.009467339 |

The physical wavefunction gauge error falls from
1.819e-05 to 1.141e-06
under fourfold timestep refinement. The fine physical triplet-population difference
between gauges is 1.439e-07. The SOC-disabled trajectory
has triplet population 2.879e-32, versus exactly
zero for the grid exponential. Exact same-runtime replay passes.

The fine wavefunction error against the grid exponential remains about 9.96e-3,
while fine triplet-population error is about 8.6e-5. These are measured residual
errors of this particular Gaussian/grid campaign, not estimates of molecular
accuracy. The committed record is local evidence; CI writes its own fresh result.

## Limits and next milestone

This validates one packet in an analytic one-dimensional model over a short
physical time, despite being longer than the preceding smoke test. It does not
establish continuum accuracy, basis completeness, multidimensional molecular SOC,
nonzero-curvature dynamics, adaptive topology, stiff integration, or live ab-initio
provider accuracy. The stable package remains v0.27.0. A smaller timestep alone
cannot remove the measured Gaussian projection defect.

The next data milestone remains an audited adapter from package-level v0.29.4
oracle records into the top-level v0.29.5 certification harness, with physical
symmetry re-evaluations and trajectory-level split/leakage checks. Larger basis
and multidimensional trajectory campaigns remain separate scientific work.
