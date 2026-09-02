# v0.10 Build Validation Report

Validated on 2026-08-12.

## Source syntax validation

All 135 Python files under the repository parse successfully with
Python's `ast` parser.

## Automated regression suite

```text
[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m [ 64%]
[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m                                  [100%][0m
[32m[32m[1m111 passed[0m[32m in 3.09s[0m[0m
```

The v0.10 tests add coverage for:

- Gaussian Wigner covariance and deterministic sampling;
- ensemble statistics;
- exact-grid campaign construction;
- managed multidimensional convergence-surface construction;
- explicit benchmark acceptance criteria;
- non-additive sensitivity/error budgets;
- reduced electronic density matrices;
- unitary density-matrix frame transformations;
- global-diabatic analytic CI observables;
- optional repeated spawning capable of growing beyond one child.

All inherited v0.1-v0.9 tests remain passing.

## Representative executable examples

### 31 — common observable / reduced-density benchmark

```text
v0.10 exact versus graph-Gaussian reduced electronic density
-----------------------------------------------------------
exact diabatic populations:   [0.22600611 0.77399389]
managed diabatic populations: [0.96114092 0.03885908]
population L2 error:           1.0396376243628043
density Frobenius error:       1.0456430442559244
exact purity:                  0.6762081969769371
managed purity:                0.9999999963862384
exact linear entropy:          0.3237918030230629
managed linear entropy:        3.613761623277867e-09

A norm-conserving Gaussian run may still fail this observable test.
```

### 32 — exact grid x timestep surface

```text
Exact grid x timestep surface
-----------------------------
N= 32  dt=  0.0025  P=[0.40202485 0.59797515]  |P-P_finest|=5.549284e-04  norm=1.000000000000005
N= 32  dt=   0.005  P=[0.40202485 0.59797515]  |P-P_finest|=5.549256e-04  norm=1.000000000000002
N= 32  dt=    0.01  P=[0.40202486 0.59797514]  |P-P_finest|=5.549142e-04  norm=0.999999999999998
N= 48  dt=  0.0025  P=[0.40349821 0.59650179]  |P-P_finest|=1.528719e-03  norm=0.999999999999804
N= 48  dt=   0.005  P=[0.40349821 0.59650179]  |P-P_finest|=1.528721e-03  norm=0.999999999999908
N= 48  dt=    0.01  P=[0.40349822 0.59650178]  |P-P_finest|=1.528730e-03  norm=0.999999999999953
N= 64  dt=  0.0025  P=[0.40241724 0.59758276]  |P-P_finest|=0.000000e+00  norm=1.000000000000046
N= 64  dt=   0.005  P=[0.40241724 0.59758276]  |P-P_finest|=3.935559e-09  norm=1.000000000000022
N= 64  dt=    0.01  P=[0.40241726 0.59758274]  |P-P_finest|=1.969305e-08  norm=1.000000000000010

Finest candidate: {'grid_n': 64, 'dt': 0.0025, 'norm': 1.0000000000000464, 'populations': array([0.40241724, 0.59758276])}
The finest row is only a candidate reference; inspect its neighboring refinement differences before calling it converged.
```

### 33 — managed convergence surface

```text
Managed graph-Gaussian convergence surface
------------------------------------------
dt= 0.0100  SPA=0  Nmax=2  block=0.9000  basis_used=2  spawns=1  cond_max=1.0028e+00  norm_err=2.285e-07
dt= 0.0100  SPA=0  Nmax=2  block=0.9999  basis_used=2  spawns=1  cond_max=1.0028e+00  norm_err=2.285e-07
dt= 0.0100  SPA=0  Nmax=4  block=0.9000  basis_used=2  spawns=1  cond_max=1.0028e+00  norm_err=2.285e-07
dt= 0.0100  SPA=0  Nmax=4  block=0.9999  basis_used=4  spawns=3  cond_max=1.9602e+04  norm_err=2.265e-07
dt= 0.0100  SPA=1  Nmax=2  block=0.9000  basis_used=2  spawns=1  cond_max=1.0028e+00  norm_err=2.285e-07
dt= 0.0100  SPA=1  Nmax=2  block=0.9999  basis_used=2  spawns=1  cond_max=1.0028e+00  norm_err=2.285e-07
dt= 0.0100  SPA=1  Nmax=4  block=0.9000  basis_used=2  spawns=1  cond_max=1.0028e+00  norm_err=2.285e-07
dt= 0.0100  SPA=1  Nmax=4  block=0.9999  basis_used=4  spawns=3  cond_max=1.9602e+04  norm_err=2.265e-07
dt= 0.0050  SPA=0  Nmax=2  block=0.9000  basis_used=2  spawns=1  cond_max=1.0029e+00  norm_err=7.147e-09
dt= 0.0050  SPA=0  Nmax=2  block=0.9999  basis_used=2  spawns=1  cond_max=1.0029e+00  norm_err=7.147e-09
dt= 0.0050  SPA=0  Nmax=4  block=0.9000  basis_used=2  spawns=1  cond_max=1.0029e+00  norm_err=7.147e-09
dt= 0.0050  SPA=0  Nmax=4  block=0.9999  basis_used=4  spawns=3  cond_max=1.9960e+04  norm_err=6.408e-09
dt= 0.0050  SPA=1  Nmax=2  block=0.9000  basis_used=2  spawns=1  cond_max=1.0029e+00  norm_err=7.147e-09
dt= 0.0050  SPA=1  Nmax=2  block=0.9999  basis_used=2  spawns=1  cond_max=1.0029e+00  norm_err=7.147e-09
dt= 0.0050  SPA=1  Nmax=4  block=0.9000  basis_used=2  spawns=1  cond_max=1.0029e+00  norm_err=7.147e-09
dt= 0.0050  SPA=1  Nmax=4  block=0.9999  basis_used=4  spawns=3  cond_max=1.9960e+04  norm_err=6.429e-09

The stored state-label populations are an internal convergence proxy. Use example 31 for a rigorous exact/Gaussian population comparison.
```

### 34 — Wigner sampling

```text
Gaussian Wigner initial-condition sampling
------------------------------------------
target q mean: [-0.6   0.25]
sample q mean: [-0.5908186   0.25081154]

target p mean: [10.  0.]
sample p mean: [ 1.00041281e+01 -1.33655719e-03]

target Cov(q):
 [[0.35714286 0.        ]
 [0.         0.35714286]]
sample Cov(q):
 [[ 0.35395772 -0.00587429]
 [-0.00587429  0.36179946]]

target Cov(p):
 [[0.7 0. ]
 [0.  0.7]]
sample Cov(p):
 [[ 0.68301989 -0.00752692]
 [-0.00752692  0.70163067]]
```

### 35 — compact error budget

```text
v0.10 compact sensitivity/error budget
----------------------------------------
total_vs_exact: 1.0242756986247297
exact_discretization_proxy: 5.760829006233664e-05
managed_timestep_proxy: 0.00038484812098345305
spa_truncation_proxy: 0.0008869391506525756
spawn_threshold_proxy: 0.0
basis_size_proxy: 0.5346564941377776
dominant_proxy: basis_size_proxy

Acceptance result:
{'passed': False, 'checks': {'norm': True, 'population_sum': True, 'conditioning': True, 'pruning_loss': True, 'exact_reference_population': False}, 'metrics': {'final_populations': [0.01298193428853648, 1.0194141115734356], 'final_norm': 0.9999999194441826, 'max_norm_error': 8.055581735000317e-08, 'max_condition_number': 19927.274562082956, 'max_basis_size': 4, 'spawn_count': 3, 'prune_count': 0, 'total_pruning_loss': 0.0, 'max_spa1_relative_correction': 2.2487045511537406e-05, 'population_l2_vs_reference': None, 'observed_populations': [0.9502784027694906, 0.049721597230509486], 'observed_population_l2_vs_reference': 1.0242756986247297}}

The sensitivity terms are not independent statistical errors and are not added in quadrature.
```

## Scientific status

The compact near-CI campaign intentionally **does not pass** the configured
exact-reference population criterion.

The current small graph-Gaussian basis remains much too pure compared with the exact
wavepacket and does not reproduce the exact reduced electronic density for this
stronger passage.

The build itself is validated; the demanding physical benchmark exposes a method
limitation that should guide the next release.

## PySCF status

The explicit PySCF backend, CASSCF many-electron state tracking, graph gauge transport,
and incremental graph layers from earlier releases remain included.

PySCF is not installed in the build environment, so real PySCF binary calculations
are not executed here.  `V10_PYSCF_BENCHMARK_PROTOCOL.md` specifies the convergence
sequence for a PySCF-enabled machine.
