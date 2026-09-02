# v0.10 PySCF Benchmark Protocol

The PySCF backend inherited from v0.5-v0.8 should be introduced only after the analytic
benchmark campaign is understood.

## Stage 1 — single-point electronic validation

At representative geometries verify:

- SCF convergence;
- SA-CASSCF convergence;
- active orbitals;
- active-electron count;
- state-average weights;
- state character;
- analytical gradients;
- analytical NAC convention.

## Stage 2 — geometry-step tracking refinement

For a path sampled with step $\Delta q$, repeat with

$$
\Delta q,\quad
\Delta q/2,\quad
\Delta q/4.
$$

Compare:

- many-electron overlap matrices;
- assigned root permutations;
- assigned overlap magnitudes;
- overlap unitarity defects;
- analytic projected NAC;
- overlap-derived directional NAC.

## Stage 3 — state-manifold convergence

Increase the number of averaged/tracked states.

A two-state model is inadequate if substantial overlap leaks into a third state.

The state-overlap singular values and unitarity defect provide direct diagnostics.

## Stage 4 — active-space convergence

Compare at least two scientifically defensible active spaces.

Track whether:

- state ordering changes;
- state character changes;
- gradients change;
- NAC location/magnitude changes.

A dynamics difference caused by a different active space is an electronic-structure
difference, not a Gaussian convergence error.

## Stage 5 — basis-set convergence

Repeat representative electronic points with the next practical basis level.

The goal is not to make every trajectory prohibitively expensive.

The goal is to establish that the selected cheaper level does not qualitatively move
the crossing or coupling region.

## Stage 6 — cache provenance

Every cached point used in a benchmark campaign should identify:

```text
geometry
basis
charge/spin
SCF reference
ncas
nelecas
number of states
state weights
SCF tolerance
CASSCF tolerance
CASSCF gradient tolerance
ETF convention
PySCF version
tracking/gauge information
```

Do not merge cached points generated under different electronic contracts into one
convergence table.

## Stage 7 — nuclear dynamics convergence

Only after Stages 1-6:

- vary Gaussian timestep;
- vary SPA order;
- vary spawning action;
- vary overlap blocker;
- vary maximum basis;
- vary pruning tolerances;
- increase initial-condition count.

The final uncertainty statement should distinguish electronic-structure sensitivity
from nuclear/Gaussian sensitivity.
