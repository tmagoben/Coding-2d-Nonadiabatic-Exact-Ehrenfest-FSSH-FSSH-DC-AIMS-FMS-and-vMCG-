# v0.5 Validation Contract

## Backend-independent data contract

A Cartesian electronic-structure point must satisfy:

- finite state energies;
- one Cartesian gradient per state;
- one Cartesian NAC vector per ordered state pair;
- exact zero diagonal NAC;
- antisymmetric real NAC tensor;
- atom count matching the requested geometry;
- atomic mass count matching the geometry.

## Coordinate projection

For a constant linear map $R=R_0+Jq$:

- projected gradients equal $J^Tg_R$;
- projected NACs equal $J^Td_R$;
- generalized mass matrix equals $J^TM_RJ$;
- $M_q$ is symmetric positive definite for independent modes.

## Gaussian local matrix elements

- analytic overlap matches numerical 2D quadrature;
- analytic gradient matrix element matches quadrature;
- analytic kinetic matrix element matches quadrature;
- local $S$ and $H$ are Hermitian;
- moving-basis identity satisfies
  $\dot S=T^{basis}+T^{basis\dagger}$.

## Direct spawned dynamics

- $C^\dagger SC$ is the monitored norm;
- zero-amplitude basis insertion is continuous;
- child momentum satisfies the generalized-mass energy equation;
- deterministic inputs produce deterministic spawn decisions;
- an inserted child receives amplitude only through coupled propagation.

## Explicit PySCF backend

The software test suite verifies the backend call contract with a deterministic fake
PySCF implementation even when PySCF is not installed:

- molecule built with `unit="Bohr"`;
- requested RHF/ROHF reference is used;
- SCF tolerances and macro settings are forwarded;
- `state_average_(weights)` is called;
- state-specific gradient calls are made;
- internal NAC convention requests `state=(J,I)`;
- dynamics NAC always uses `mult_ediff=False`;
- configured `use_etfs` is forwarded;
- optional scaled NAC calls use `mult_ediff=True`;
- convergence failures raise rather than silently return data.

## Real PySCF acceptance before research use

A machine with PySCF installed must additionally run the backend-specific validation
example and check:

1. PySCF version;
2. RHF/ROHF convergence;
3. SA-CASSCF convergence;
4. active-space choice;
5. state ordering/character;
6. finite-difference state gradients at selected geometries;
7. NAC antisymmetry and convention;
8. translation/ETF choice;
9. continuity along the intended trajectory coordinate.

The build environment used to create v0.5 does **not** contain PySCF, so the real
backend is source/API validated and fake-backend integration tested here, but not
numerically executed against the PySCF binary in this environment.
