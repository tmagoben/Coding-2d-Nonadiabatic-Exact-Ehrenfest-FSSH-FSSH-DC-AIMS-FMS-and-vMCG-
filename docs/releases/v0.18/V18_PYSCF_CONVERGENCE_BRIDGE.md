# v0.18 PySCF / Molecular Convergence Bridge

v0.18 still uses the analytic two-state LVC benchmark for release acceptance.

However, the convergence architecture was chosen specifically to make a future
molecular direct-dynamics benchmark more defensible.

## 1. Full-wavefunction metrics do not transfer directly to arbitrary molecules

The analytic model has a globally defined diabatic two-component grid wavefunction.

For a PySCF molecular trajectory, one generally does not have an independent exact
high-dimensional grid wavefunction.

Therefore the v0.18 full-wavefunction metrics become a **small-model validation
standard**, not a directly available molecular observable.

A molecular release will need surrogate convergence quantities such as:

```text
projected TDSE residual
basis-to-basis self convergence
reduced electronic density convergence
nuclear moments
trajectory observables
local operator/matrix audits
```

## 2. Physical-time normalized controls do transfer

The following v0.18 change should carry directly into molecular calculations:

```text
defect checks specified in physical time
minimum adaptation separation specified in physical time
prune age specified in physical time
audit cadence specified in physical time
cost horizon specified in physical time
```

Without this, changing `dt` changes both the integrator and the adaptive algorithm,
making a convergence study ambiguous.

## 3. Batched residual ranking is directly useful

Electronic-structure direct dynamics can have expensive candidate dictionaries.

The v0.18 batching strategy bounds candidate-grid memory without changing the residual
projection algebra.

A future molecular implementation can also batch:

```text
electronic cache queries
candidate operator estimates
local residual contractions
```

## 4. Sampled audits are the right direction, but not yet sufficient

A full dense molecular S/H/T rebuild every audit interval would be impractical.

v0.18 replaces normal dense audits with sampled omitted-edge checks, while retaining
dense sentinels only because the analytic benchmark is small.

A molecular successor should calibrate a hierarchy such as:

```text
frequent:
    local score budget
    residual diagnostics

medium cadence:
    priority + random omitted-edge samples

rare:
    neighborhood-dense audits

validation-only:
    full dense audit on small molecular models
```

## 5. Electronic gauge consistency remains mandatory

Any molecular $H_{ij}$, derivative coupling, or future SOC matrix entering an edge
score must be represented in consistent electronic frames.

Thus the v0.6-v0.8 state-tracking and gauge-graph infrastructure remains a prerequisite
for the later molecular/SOC releases.

## 6. Why SOC is still deferred

v0.18 exposes a nontrivial remaining full-wavefunction error even though reduced
electronic errors are already small.

Adding SOC now would make it harder to distinguish:

```text
basis incompleteness
sparse truncation
electronic gauge error
SOC physics
```

The cleaner roadmap is still:

```text
v0.18  convergence completeness
v0.19  molecular/direct-dynamics integration
v0.20  controlled analytic SOC dynamics
v0.21  ab-initio SOC bridge
```
