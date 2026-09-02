# v0.11 Basis-Completeness Benchmark Campaign

## 1. Run the compact release campaign

```bash
python examples/39_v011_basis_ladder.py
```

or in Python:

```python
from gaussian_dynamics import run_v011_release_benchmark

result = run_v011_release_benchmark()
```

The campaign is intentionally more expensive than earlier toy examples because it
repeats the strong CI passage over multiple basis sizes and ablations.

---

## 2. Exact reference

The benchmark uses the same v0.10 exact reference:

```text
grid_n = 64
dt     = 0.0025
```

with the global-diabatic reduced electronic density as the primary observable.

The exact reference should still be interpreted together with the v0.10 exact
grid/timestep surface.

---

## 3. v0.10 baseline

The inherited v0.10 managed settings are rerun on the same initial condition.

This is important: the v0.11 comparison is not against a remembered number copied from
documentation.

The old algorithm and new algorithm are executed in the same repository/environment.

---

## 4. v0.11 reference settings

The compact v0.11 reference uses

```text
dt                         = 0.005
SPA                        = 1
integrated spawn action    = 1e-4
max basis                  = 10
max generation             = 5
children per event         = 2
position shifts            = 0, +0.05, -0.05
width scales               = 0.65, 1.0, 1.55
momentum directions        = NAC, parent momentum
existing overlap blocker   = 0.9999
sibling overlap blocker    = 0.995
```

The search is deterministic.

---

## 5. Basis ladder

Run

```text
Nmax = 2, 4, 6, 8, 10
```

without changing the remaining v0.11 parameters.

For every row inspect:

```text
population error
full density error
purity
norm drift
condition number
actual basis size
```

Do not infer convergence from basis size alone.

A larger but ill-conditioned basis may be less reliable.

---

## 6. Position-search ablation

Set

```python
position_shifts=(0.0,)
```

This reduces the child search to same-position spawning.

The difference from the reference quantifies the importance of allowing a small
nonclassical position adjustment.

---

## 7. Width-bank ablation

Set

```python
width_scales=(1.0,)
```

All children inherit the parent width.

Compare this with

```python
width_scales=(0.65,1.0,1.55)
```

to determine whether shape diversity adds genuine basis completeness or merely
numerical complexity.

---

## 8. Child multiplicity ablation

Set

```python
children_per_event=1
```

and compare with two nonredundant children per event.

This tests whether one local phase-space branch is enough to represent transferred
wavepacket structure.

---

## 9. Acceptance interpretation

The default v0.11 acceptance checks are

```text
population L2 error      <= 0.05
density Frobenius error  <= 0.10
purity error             <= 0.05
maximum norm drift       <= 0.01
condition number         <= 1e6
```

These are release-regression thresholds for this analytic benchmark.

They are **not universal chemical-accuracy standards**.

The full density criterion is intentionally stricter than population-only matching.

---

## 10. What to do if populations converge but density does not

Suppose

```text
population error   small
purity error       small
full density error still significant
```

Then the remaining problem is likely electronic coherence/phase structure rather than
gross branching probability.

The next development should then target:

- phase accuracy of coupled Gaussian amplitudes;
- more accurate off-diagonal Hamiltonian matrix elements;
- spawn timing through the coupling window;
- or variational Gaussian motion.

Do not simply add more TBFs until one scalar population happens to match.


## 11. Campaign execution and memory use

The saved release result is assembled from independent calculations.

The base call

```python
run_v011_release_benchmark(include_ablations=False)
```

computes the exact reference, v0.10 baseline, and v0.11 basis ladder.

The release ablations are run independently with `run_v011_case(...)` and merged into

```text
results/v011_basis_completeness_campaign.json
```

This prevents several large dynamic gauge graphs from remaining resident at the same
time. Examples 39 and 40 read the saved release campaign and therefore execute
quickly. Example 41 deliberately recomputes the base campaign.
