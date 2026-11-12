# v0.30.6 development: dormant-shape propagation and audited activation

v0.30.5 could admit a zero-population child but could not evolve that adaptive
basis. v0.30.6 keeps every real and imaginary electronic coefficient active while
omitting dormant shape columns from the real TDVP metric. Accepted midpoint steps
copy dormant centers, momenta, width matrices and chirp matrices exactly. The
all-active path retains v0.30.2 behavior; this does not change the v0.27.0 release.

## Lifecycle contract

Each forward step propagates, increments packet ages, audits activation, then
evaluates scheduled topology. Topology is deferred while any shape is dormant:
the v0.30.5 residual oracle uses a fully active tangent, so it must not score the
constrained state. Once all shapes are active, v0.30.5's merge/prune/spawn ordering,
projection and full-Hamiltonian energy gates apply, with at most one topology
event per checkpoint. Newborn shapes begin dormant with age zero and no seeded
electronic amplitude. All coefficient components remain free to evolve.

Activation processes dormant blocks in stable packet-ID order. It requires:

1. Coefficient row norm squared at least the existing population threshold.
2. A compatible TDVP solve with the proposed shape block included.
3. Full rank and acceptable conditioning of the new shape tangent after projecting
   out the existing retained real tangent; singular values are squared into metric
   units before applying the existing tangent cutoffs.
4. An acceptable ratio of trial to base shape-coordinate velocity norms.

For the weighted complex Jacobian $J_A$ of the currently active parameters, form
$X_A=[\mathrm{Re}\,J_A;\mathrm{Im}\,J_A]$. If $U_r$ spans its retained left
singular vectors, the proposed shape columns are tested through

$$
\widetilde X=(I-U_rU_r^T)X_{\mathrm{new}},\qquad
\lambda_k=\sigma_k(\widetilde X)^2.
$$

The new block needs all $2d+d(d+1)$ shape directions: centers, momenta, symmetric
log widths and symmetric chirps. This projected spectrum measures independent
shape freedom after existing coefficient and shape directions are removed.

The velocity ratio excludes electronic coefficient derivatives because those
depend on the electronic gauge. The population is a coefficient norm diagnostic,
not a physical population assigned to one nonorthogonal Gaussian. Activation
changes only the mask, never the physical state. If all shapes are dormant and the
base shape speed is zero, the velocity gate conservatively blocks nonzero trial
shape motion. Ordinary newborns have an active parent.

There are no hidden retries, timestep changes, renormalization or amplitude seeds.
Failed solves raise without mutating the input. Small post-activation timesteps
are supplied explicitly by the caller. Receipts bind the complete states, masks,
IDs, ages, schedule, settings, provider identity, midpoint and metric residuals,
activation spectra and shape speeds, and topology decisions. Replay recomputes
all physical stages and compares exact receipts in the same numerical runtime;
fingerprints alone are not proof of correctness.
Only identified TDVP metric failures defer a proposed activation. Provider,
transport and input-validation failures propagate to the caller.

## Evidence and controls

Run from committed implementation source:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python examples/163_recompute_v0306_adaptive_propagation.py
```

The fresh campaign uses the v0.30.2 scalar zero-potential model, mass 2, thirteen
grid points on [-3, 3], and its full discrete Hamiltonian. An initial v0.30.5
admission at t=0 is followed by this explicit schedule:

| Count | dt (atomic units) | Cumulative time |
|---:|---:|---:|
| 12 | 0.02 | 0.24 |
| 4 | 0.001 | 0.244 |
| 4 | 0.004 | 0.26 |
| 2 | 0.0001 | 0.2602 |

The packet cap is three and topology is scheduled every step. All activation,
projection and solver gates retain their defaults. The campaign requires two
unseeded spawns (including the initial admission), actual coefficient growth,
an accepted activation, subsequent motion of that shape, exact dormancy, complete
metadata continuity, and physical replay of the initial admission and all 22
steps. Saved fields also compare the endpoint with direct finite-grid matrix
exponentiation and a fixed one-packet baseline on the identical time schedule.

Regression budgets at the propagation endpoints are 2e-7 for maximum finite-grid
norm drift and 1e-7 hartree for energy drift; event projection and energy changes
have their separate admission gates. The endpoint wavefunction error budget is
0.002 against the same
discrete Hamiltonian. These are bounded fixture checks, not time or spatial
convergence certificates. The evidence gate evaluator derives acceptance from
measurements, spectra and state bindings; corruption tests ensure that stored
success flags cannot override them.

Independent tests cover generalized Crank–Nicolson evolution with all shapes
frozen, including the complete four-component analytic SOC spinor in fixed and
coordinate-dependent electronic gauges; masked scalar-gauge evolution;
activation covariance; packet permutation; gate rejection; invalid inputs;
failed-step immutability; and rehashed false replay receipts.

The [committed evidence](results/v0306_adaptive_propagation_evidence.json) was
generated from source revision `b1d48fb4fdfc071a868fe52391c92aff3dd629c2` using
Python 3.12.14, NumPy 2.5.2 and SciPy 1.18.0. All 16 campaign gates pass. The initial
admission and all 22 adaptive steps replay exactly in that runtime. The first child
activates at t=0.26 with velocity amplification 99.3608; a second child is admitted
at that checkpoint, and the trajectory continues to t=0.2602.

| Measurement | Observed value |
|---|---:|
| Adaptive weighted wavefunction error against the finite-grid reference | 0.001221624 |
| Fixed one-packet error on the same timestep schedule | 0.002387399 |
| Maximum sampled norm drift | 5.99031e-8 |
| Maximum sampled energy drift (hartree) | 2.34115e-8 |
| Maximum scaled midpoint residual (acceptance limit 1) | 0.693201 |

The adaptive endpoint error is about 49% smaller in this fixture. This comparison
does not establish general adaptive accuracy or convergence through topology
events. There are 69 new tests, including 22 evidence checks; CI regenerates the
lifecycle campaign independently from its checked-out source.

## Remaining work

The long lifecycle campaign is scalar and one-dimensional. The SOC control tests
coefficient propagation, not a complete adaptive SOC trajectory. General
coordinate covariance of activation, repeated merge/prune trajectories,
multi-initial-condition adaptive SOC, timestep refinement through topology events,
spatial convergence, and live molecular-provider validation remain open. Exact
replay is runtime-specific. Dataset certification and training eligibility remain
false; v0.30.5 export labels do not imply adaptive-trajectory certification.

The next scientific milestone is a convergence-controlled adaptive SOC campaign,
including event-time sensitivity and a broader trajectory set. A robust nonlinear
solver or explicit retry policy may be needed for stiff newly activated shapes.
