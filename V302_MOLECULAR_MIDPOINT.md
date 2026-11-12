# v0.30.2 development — molecular-frame Gaussian midpoint trajectories

The curved Gaussian TDVP now propagates coefficients, centres, momenta, full
positive-definite widths and symmetric chirps using the v0.30.1 covariant grid
Hamiltonian. This is a fixed-basis development build. The validated package release
remains v0.27.0; the build does not certify a live ab-initio SOC trajectory.

## Equation and acceptance

Let `theta` be the inherited real parameter vector (real/imaginary coefficients,
centres, momenta, orthonormal symmetric log-width and chirp coordinates). Each step
solves

$$
\theta_{n+1}-\theta_n=\Delta t\,v((\theta_n+\theta_{n+1})/2),
$$

where `v` comes from the full real TDVP metric and SVD solve. Every evaluation uses
both nuclear kinetic energy and the provider's electronic Hamiltonian. The provider
must be stationary and pass the v0.30.1 grid fingerprint checks.

The solver starts from an Euler predictor and applies fixed-point iteration. It
accepts only when

$$
\max_i \frac{|\theta_{n+1,i}-\theta_{n,i}-\Delta t\,v_i|}
 {a+r\max(|\theta_{n,i}|,|\theta_{n+1,i}|)}\le 1.
$$

Defaults are `a=2e-10`, `r=2e-9`, tangent finite-difference step `2e-6`, and at most
24 iterations. Failure raises an exception; the caller may retry a smaller step.
There is no automatic time-step substitution or renormalization. Negative time
steps use the same equation. Widths remain positive through the existing matrix
exponential parameterization. Rank, conditioning, SVD residual and incompatible
null-space gates remain those of the inherited tangent solver.

An accepted step contains before/after states, nonlinear and tangent-solve residuals,
metric rank/conditioning, TDSE projection defect, finite-grid norm and energy,
operator and bridge fingerprints, and an overlap-versus-NAC audit for each centre
displacement. Norm and energy are diagnostics, not separate universal acceptance
thresholds. Packet membership and atlas patches stay fixed.

## Replay and cost

`run_molecular_midpoint_v302` returns a tuple of accepted steps.
`serialize_molecular_trajectory_v302` checks receipts, state consistency and step
continuity before emitting canonical JSON. Actual width matrices are serialized;
reconstructing them through a second log/exp would introduce rounding changes.
`replay_molecular_trajectory_v302` reruns each nonlinear solve and transport audit,
then requires exact receipt equality. It detects altered results even when the JSON
hash has been recomputed. Hashes are integrity checks, not signatures.

Exact replay is a same-numerical-runtime guarantee. A different BLAS, dependency or
interpreter may change low-order bits and correctly fail exact replay. Agreement
across runtimes must be assessed numerically, separately from byte identity.

Wilson segments are cached only within one tangent evaluation. Coefficient,
momentum, width and chirp perturbations reuse the same centre-to-grid transport;
centre perturbations retain distinct segments. Cache entries never survive the
evaluation. A direct uncached moving-SOC tangent solve is a regression control.
The finite-grid tensor-product and dense parameter-metric costs remain; the
fixed-point iteration is not a robust stiff-regime integrator.

## Reproduction and current evidence

```bash
python -m pip install -e ".[dev]"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -m pytest -q tests/test_v0301_molecular_hamiltonian.py tests/test_v0302_molecular_midpoint.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python examples/159_recompute_v0302_molecular_midpoint.py
```

The example writes `results/v0302_molecular_midpoint_evidence.json`. Its analytic
zero-potential fixture uses one packet, mass 2, 13 interior grid points from -3 to 3,
initial centre 0.1, momentum 0.3, width 1.1, chirp 0.08, and final time 0.16.
All seven campaign gates pass in the recorded local run:

| Time step | Steps | Norm drift | Energy drift (hartree) |
| --- | ---: | ---: | ---: |
| 0.08 | 2 | 5.98e-10 | 9.18e-7 |
| 0.04 | 4 | 1.17e-10 | 2.30e-7 |
| 0.02 | 8 | 5.50e-12 | 5.74e-8 |

The wavefunction self-convergence order is 1.9986. Signed-step reversal and exact
same-runtime replay pass. The finest trajectory's weighted wavefunction difference
from independent sparse-exponential grid evolution is about 1.46e-3. This remaining
projection error is much larger than the time-refinement difference: temporal
convergence does not establish basis completeness or continuum accuracy.

Regression tests additionally cover two-Gaussian propagation/replay, short analytic
singlet/triplet SOC propagation in coordinate-dependent gauges, correlated-width
serialization, modified receipts, re-signed false results, operator mismatch,
invalid settings, and a deliberately unconverged solve. Midpoint in a nonlinear
parameter chart is not exactly gauge covariant at finite step size; the short SOC
comparison uses an explicit numerical tolerance.

## Cumulative validation and historical evidence

The full CI suite now runs on pull requests as well as main and installs the pinned
PySCF scientific runtime. Three disjoint file groups run concurrently: core tests,
earlier cumulative campaigns, and later cumulative campaigns. No test file is
excluded. Development and cumulative jobs upload JUnit reports;
the development job recomputes the midpoint evidence.

Seventeen JSON files referenced by historical packaging tests were missing from
Git. Their unchanged bytes are packaged in `tests/fixtures/v0300_historical_evidence.tar.xz`,
from the verified v0.30.0 checkpoint archive. At pytest startup, missing files are
extracted after validating the exact member list and SHA-256 hashes. Existing local
result files are never replaced.
`results/v0302_historical_evidence_recovery.json` records the archive SHA-256 and
each file's SHA-256. These are historical fixtures, not current execution results.
Runtime identity and scientific compatibility are distinct; pinning dependencies
does not make a hosted runner's Python binary identical to the sealed release.
The cumulative regression assertions independently recompute both profiles, require
every scientific gate, and require the identity gate to agree with the measured
runtime. On a noncanonical binary, exactly that identity gate must fail and the
original campaign's `acceptance.passed` must remain false. No campaign implementation,
release threshold or identity policy is changed. JUnit success means these regression
contracts pass; it does not confer canonical release acceptance. Per-version CI
artifacts preserve that distinction. Negative controls reject any additional failed
gate, a forged identity report, or a false claim of canonical release success.

## Remaining scientific work

The [v0.30.3 campaign](V303_SOC_TRAJECTORIES.md) extends SOC validation to a
ten-times-longer analytic singlet–triplet trajectory with temporal and gauge
refinement. Larger and multidimensional packet campaigns,
adaptive topology, changing patches, robust stiff integration, and live molecular
provider accuracy remain unvalidated by this slice. The next data milestone is an
audited adapter between the package-level v0.29.4 records and the top-level
v0.29.5 certification harness, followed by physical symmetry re-evaluations and
trajectory-level split/leakage checks.
