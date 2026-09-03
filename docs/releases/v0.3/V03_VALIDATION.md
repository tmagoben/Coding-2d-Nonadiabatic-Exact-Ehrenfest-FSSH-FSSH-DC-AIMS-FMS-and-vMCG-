# v0.3 Validation Contract

Core provider tests are backend independent and run in CI.

The optional PySCF adapter has an additional validation requirement before research
use:

1. verify RHF and SA-CASSCF convergence;
2. verify the requested active space/state manifold;
3. compare state gradients against finite differences for a small geometry set;
4. verify NAC antisymmetry/convention;
5. compare the projected 1D gradient and NAC against explicit finite differences along
   the chosen coordinate;
6. record PySCF version, basis, active space, weights, convergence thresholds, and
   coordinate definition;
7. cache raw Cartesian gradients/NACs as well as their projected scalar values.

A backend passing software tests is not automatically chemically adequate. Electronic
structure level, active space, and state selection remain scientific convergence
questions.
