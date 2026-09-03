# v0.8 Validation Contract

Version 0.8 must satisfy all earlier v0.1-v0.7 regressions and additionally:

- overlap/local-diabatic electronic propagation is covariant under independent endpoint $U(N)$ gauges;
- explicit-NAC and overlap propagation converge to the same analytic CI result;
- both electronic propagators preserve norm;
- incremental temporal graph links are unitary after polar projection;
- dynamic center/centroid graph growth creates the expected cycles;
- the moving-basis connection satisfies $\dot S=T+T^\dagger$ numerically;
- metric correction preserves the seed anti-Hermitian part;
- the time-dependent generalized coefficient step preserves $C^\dagger SC$ to integration tolerance;
- dynamic spawning inserts exactly zero child amplitude at birth;
- later coupled propagation produces nonzero child amplitude;
- graph node count grows monotonically in the dynamic benchmark;
- the incremental many-electron snapshot graph computes only requested new overlap edges;
- raw snapshot-edge singular values and unitarity defects remain inspectable.

Passing these tests establishes numerical and gauge consistency of the v0.8 benchmark. It does not establish chemical convergence for a real molecule.
