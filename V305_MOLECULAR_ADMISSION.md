# v0.30.5 — molecular topology admission and full-spinor dataset adapter

This development increment follows merged v0.30.4. The stable distribution remains
v0.27.0. It adds an explicit single-checkpoint topology controller and a separate
dataset schema; it does not propagate an adaptive trajectory or train a surrogate.

## Final admission

`evaluate_molecular_admission_v305` first recomputes every v0.30.4 residual candidate
from the complete molecular grid Hamiltonian. It then applies at most one topology
event, in the inherited order: merge, prune, spawn, no event. Packet IDs, ages,
shape-activity flags and the next unused packet serial are validated and recorded.
The caller chooses when to invoke the checkpoint; this is not a timestep scheduler.

Merge and prune use complete curved-wavefunction SVD projection. Spawning selects
the largest passing residual score, with stable packet IDs and candidate geometry
as deterministic tie breakers. It attempts only that candidate. Every attempted
event must pass coefficient-basis rank/condition, projection loss, and the absolute
full-Hamiltonian energy-jump limit. The older curved controller omitted the spawn
energy check; the new controller evaluates it explicitly. Historical controllers
and the sealed release contracts are unchanged.

| Candidate outcome | Meaning |
| --- | --- |
| `residual_rejected` | The complete-complex residual oracle rejected it. |
| `not_selected` | Another candidate ranked higher; projection was not attempted. |
| `earlier_topology_event` | A merge or prune consumed this checkpoint's event. |
| `dormant_shape_pending` | An existing packet still has dormant shape parameters. |
| `maximum_packet_count` | The packet budget blocks spawning. |
| `projection_rejected` | The selected candidate failed rank, condition or loss gates. |
| `energy_rejected` | Projection completed but the energy jump exceeded its limit. |
| `admitted` | The selected candidate passed all admission gates and was appended. |

Only `admitted` has `final_admitted=true`. A nonselected or blocked candidate is not
an intrinsic negative example of physical admissibility. Residual acceptance stays
separate from final admission in every record. The candidate spans the complete
complex electronic spinor; `target_state=null` is intentional.

An admitted packet starts with age zero and dormant shape parameters. Projection
preserves the existing physical wavefunction up to its measured error, so newborn
coefficients may initially be numerically zero. This module neither seeds an
artificial population nor claims that newborn activation or propagation has passed.
A no-event decision returns the original state object without numerical rewriting.
No packet ages advance at the instantaneous checkpoint.

## Exact replay and export

`replay_molecular_admission_v305` reconstructs the source state and lifecycle
metadata, recomputes all residual labels and attempted events, and requires an exact
match of the complete receipt and final state. Rehashing a changed decision, energy,
activation mask or candidate list cannot make it replay. Exactness is a same-runtime
contract, distinct from tolerance-based physical symmetry comparisons.

`materialize_admission_dataset_v305` accepts serialized checkpoints with their
actual bridges and Hamiltonians. It physically replays every checkpoint before
export. Rows separate predictor features, oracle targets and source provenance.
The allowlist contains Gaussian geometry scalars and the source coefficient-row
norm; it excludes parameter-metric condition numbers and residual labels from the
predictor inputs. A coefficient-row norm is not labeled a physical state population.

The new schema is `gnd-full-spinor-admission-dataset-v0.30.5-dev`. The frozen top-level
v0.29.4 state-targeted schema and v0.29.5 certification API are unchanged. Their rows
must not be mixed with this export without a future explicit migration.

All steps, candidates, control settings and electronic-frame variants of a complete
`(family_id, trajectory_id)` group share one deterministic split. Group identities
use tuple encoding rather than delimiter concatenation. Duplicate checkpoint/variant
identities fail. Canonical output order is independent of input order. The manifest
records a full source commit, checkpoint hashes, row hashes and split counts; retain
it alongside the JSONL rows. Provenance alone is not a physics certificate.

Exports explicitly set `exact_replay_verified=true`, `dataset_certified=false`,
`physical_invariance_certified=false` and `training_eligible=false`.

## Validation

From a committed implementation, using the single-thread environment in CI:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python examples/162_recompute_v0305_molecular_admission.py
```

The command refuses to attribute dirty implementation code to `HEAD`. An explicit
`--generation-commit` is available for externally verified source trees. The output
contains the complete admission evidence and a JSONL companion. Hosted CI recomputes
the campaign and uploads both files from the checked-out commit.

The campaign uses two distinct analytic SOC initial conditions on the existing
13-point grid. Both are evaluated in fixed and coordinate-dependent electronic
frames. It compares residual scores, final decisions, exported scalar features and
the physical wavefunction after admission. Projection-conditioning, packet-budget
and residual-threshold controls preserve rejected candidates. All seven checkpoints
are physically replayed before the 28 candidate rows are materialized.

The committed [evidence](results/v0305_molecular_admission_evidence.json) and
[JSONL rows](results/v0305_molecular_admission_evidence.jsonl) were generated from
implementation commit `29852181b0016489dcd6e828fc79be0d19e53176`. All **18 gates pass**;
the largest post-admission physical gauge difference is **2.37e-14**. The related
local suite passes 137 tests, with the additional actual-newborn checkpoint test
and 16 evidence-corruption controls passing separately: **154 checks in total**.
The wheel builds and all new modules import from an isolated installation.

Focused tests also cover merge/prune priority, packet-order independence of removal,
real pruning energy rejection, an injected spawn energy jump, input validation,
exact no-event behavior, complete-complex labels, rehashed corruption, provenance,
split controls and duplicate exports. CI retains all prior cumulative, integration,
midpoint, SOC and spawning campaigns and adds a fresh admission job.

## Remaining boundary

This is a fixed-grid, single-event admission audit. It does not establish continuum
convergence, Gaussian-basis completeness, live ab-initio SOC accuracy, adaptive
trajectory stability, newborn activation, or certified training data. A complete
coordinate/permutation feature audit, multi-step adaptive propagation, and a
versioned campaign-certification layer are the next acceptance steps. Gauge variants
must remain in their parent trajectory's split throughout that work.
