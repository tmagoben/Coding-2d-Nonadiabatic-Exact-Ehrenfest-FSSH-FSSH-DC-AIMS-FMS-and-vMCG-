# v0.30.7 development: adaptive SOC convergence

v0.30.6 established the scalar newborn lifecycle. v0.30.7 exercises an unseeded
birth, coefficient growth, audited shape activation and subsequent shape motion
under the complete four-state singlet/triplet Hamiltonian. It separates two
numerical questions: convergence of the midpoint integration at a fixed
adaptation schedule, and sensitivity to how frequently activation is checked.
The stable release remains v0.27.0.

## Physical checkpoints and integration substeps

`molecular_checkpoints_v307.py` integrates to each requested physical time with a
fixed active-shape mask. At that endpoint it ages each packet once, audits
activation, and applies scheduled v0.30.5 topology only if all shapes are active.
Packet ages count checkpoints, not integration substeps. Refining the integration
mesh therefore preserves the lifecycle schedule. Changing checkpoint cadence is
a different policy experiment; it also changes the checkpoint age count.

Every accepted substep retains the v0.30.6 midpoint, metric, null-space and
transport receipts. No hidden timestep retry, coefficient seed, renormalization
or relaxed acceptance threshold is introduced. A failed substep does not mutate
the input state or publish an accepted checkpoint. Full serialization binds the
schedule, states, metadata, settings and nested receipts. Replay recomputes every
substep and lifecycle decision and requires the complete receipts to match.

## Independently checked analytic tangent

The restricted `AnalyticSOCFrameBridgeV307` accepts only the declared analytic
singlet/triplet model and its one-dimensional linear-phase gauge. Its connection
is constant, so segment transport is exactly

$$W(R\leftarrow q)=\exp[-A(R-q)].$$

For the Gaussian with width $\Gamma=\exp(E)$, chirp $B$, momentum $p$ and
$x=R-q$, the shape derivatives are

$$
\partial_q(gWc)=[(\Gamma-iB)x-ip]gWc+gWAc,\qquad
\partial_p(gWc)=ixgWc,
$$

$$
\partial_E(gWc)=(1/4-\Gamma x^2/2)gWc,\qquad
\partial_B(gWc)=ix^2gWc/2.
$$

The center-connection term is essential for gauge covariance. Coefficient
columns are $gW$ and $igW$; all columns include quadrature square-root weights.
The optional bridge Jacobian hook leaves generic molecular providers on central
differences and preserves their previous fingerprints and solver gates. The
specialized bridge binds its tangent backend, model and gauge to its fingerprint
and rejects changes to that configuration.

Tests compare exact transport against independent NAC-ODE integration and frame
endpoints, compare the full grid operator and one-packet velocity against the
generic numerical bridge, check second-order central-difference agreement for a
populated two-packet tangent, and test shape velocity and physical time-derivative
gauge covariance. Provider/configuration corruption must fail closed.

Finite-difference tangent noise caused irregular small-step midpoint failures in
the initial probe. The analytic tangent removes that source of derivative noise.
The fixed-point midpoint still needs explicit small steps after activation; this
iteration does not supply a general stiff nonlinear solver.

## Bounded campaign

| Quantity | Setting |
|---|---|
| Initial conditions $(q_0,p_0)$ | $(0.1,0.3)$ and $(-0.25,0.55)$ |
| Electronic model | One singlet plus three triplet components; SOC scale 0.2 |
| Mass, width, chirp | 2, 1.1, 0.08 |
| Spatial grid | 13 points on $[-3,3]$ |
| Packet budget | Two; one v0.30.5 unseeded admission at $t=0$ |
| Fixed adaptation times | 0.02, 0.04, 0.06, 0.08 |
| Integration steps per interval | 20, 40, 80 |
| Integration timesteps | 0.001, 0.0005, 0.00025 |
| Frames | Fixed spin-diabatic and coordinate-dependent unitary gauge |
| Separate checkpoint spacings | 0.02, 0.01, 0.005, 0.0025, each with integration step 0.0005 |
| Reference | Matrix-exponential evolution of the same complete discrete Hamiltonian |

Physical fields and triplet populations in the moving frame are transformed back
to the fixed spin basis before comparison. The nonlinear gauge transformation
includes moving Gaussian centers, so discrete midpoint trajectories need not
agree exactly between gauges; their physical difference must converge away.

The offline evidence auditor recomputes endpoint errors, observed orders, gauge
errors, sampled activation brackets, triplet populations and norm/energy drifts
from saved fields and diagnostics. It does not trust saved pass flags or summary
measurements. Five distinct initial admissions and five complete trajectories (one coarse
trajectory per initial-condition/gauge pair plus the SOC-off control) are
physically replayed during generation. Finer trajectories are freshly integrated
and compared, with their complete trace digests retained.

The committed evidence stores checkpoint states and every integration substep's
diagnostics. The driver can serialize complete replayable traces; the compact
campaign summaries are not themselves substitutes for those full traces.

## Recorded results

The committed evidence was generated from source
`5c4d23a621804bdae1d872cef9f0059b32dc25dc` with Python 3.12.14,
NumPy 2.5.2 and SciPy 1.18.0. All 21 acceptance gates pass.

| Initial condition / frame | Observed timestep order | Finest-mesh finite-grid reference error |
|---|---:|---:|
| central / fixed | 1.996932 | 0.00059847263 |
| central / moving | 1.996599 | 0.00059847254 |
| displaced / fixed | 1.997991 | 0.000541519 |
| displaced / moving | 2.002008 | 0.00054151907 |

The finest-mesh physical gauge differences are 2.04234e-8 (central) and
3.82919e-8 (displaced). Both decrease at approximately second order. All fixed
checkpoint meshes activate the child at t=0.04, and its geometry subsequently
moves. The SOC-off control has triplet population 1.01225e-35 and leaves the
child dormant throughout the interval.

| Checkpoint spacing | First accepted activation | Last rejected / first accepted bracket |
|---:|---:|---:|
| 0.02 | 0.04 | (0.02, 0.04] |
| 0.01 | 0.04 | (0.03, 0.04] |
| 0.005 | 0.035 | (0.03, 0.035] |
| 0.0025 | 0.0325 | (0.03, 0.0325] |

Successive endpoint differences in this cadence study are 5.83053e-13,
5.35746e-5 and 2.71819e-5. The first pair shares the same sampled activation
time; its near-zero difference is a plateau, not evidence for an event-time
convergence order. The largest accepted scaled midpoint residual is 0.998378
against a limit of 1. No solver threshold was relaxed.

## Reproduction and limits

After committing the implementation, run with the CI thread settings:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python examples/164_recompute_v0307_adaptive_soc_convergence.py
python -m pytest -q tests/test_v0307_*.py
```

The evidence records the full source commit and numerical runtime. CI independently
regenerates the campaign and uploads it alongside all existing validation jobs.
Exact replay is scoped to the same runtime. Independent cases run in two processes
by default; `--workers 1` selects serial execution without changing a solve.
For long local runs, `--checkpoint-directory PATH` atomically saves completed
cases, keyed by source commit, numerical runtime and case inputs. Corrupted
checkpoints fail closed; failed cases are not published. CI uses fresh computation
without this optional cache.

The study establishes timestep self-convergence on a fixed, small spatial grid
and a two-packet analytic fixture. The grid-reference error also includes finite
Gaussian-basis error and should not be called temporal error. Activation times
are first accepted sampled checkpoints, not accurately located event roots.
Checkpoint refinement gives narrower sampled false/true brackets; no second-order
event-time convergence is asserted. There is no merge/prune event in this campaign.

Spatial/domain convergence, larger adaptive bases, repeated SOC spawning,
merge/prune convergence, general molecular-provider trajectories and robust stiff
solves remain open. `dataset_certified`, `training_eligible`, and
`continuum_certified` remain false.
