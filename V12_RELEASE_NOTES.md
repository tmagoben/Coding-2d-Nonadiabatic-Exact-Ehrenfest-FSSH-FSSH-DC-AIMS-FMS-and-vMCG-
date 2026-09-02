# v0.12 Release Notes

Version 0.12 is the **representation-consistency and coherence-validation release**.

## Main correction

The strong-CI exact benchmark starts from

$$
g(R)\Phi_a(R),
$$

where the electronic eigenvector varies across the nuclear coordinate.

The earlier center-based reduced-density diagnostic approximated a TBF as

$$
g(R)\Phi_a(q_i).
$$

v0.12 shows that these are materially different for the broad strong-CI packet.

For the default initial state:

```text
exact coordinate-dependent initial purity: 0.6764597760317345
center-frozen initial purity:              0.9999999999999998
center-frozen initial density error:       0.28703562527170995
```

Therefore v0.10-v0.11's full-density discrepancy should not be attributed solely to
time propagation or spawning.

## New modules

```text
lvc_exact_gaussian.py
moving_basis_v12.py
coherence_metrics.py
local_diabatic_tbf_v12.py
coherent_lvc_dynamics_v12.py
spinor_complete_lvc_v12.py
paired_basis_management_v12.py
spinor_complete_dynamics_v12.py
born_huang_grid_v12.py
born_huang_dynamics_v12.py
initial_projection_v12.py
v12_benchmark.py
```

## Main release benchmark

The nine-Gaussian projected reference gives:

```text
initial projection fidelity:             0.832276023595292
initial reduced-density error:           0.03545457994295867
projected-state dynamics density error:  0.00029022869338069174
original-target density error:           0.03500028070905269
population error:                        0.02810899300694737
coherence phase error / rad:             0.0019607485027196615
purity error:                            0.013826187810993096
norm drift:                              1.3083560634896685e-06
condition number:                        2235.290713199147
```

All configured release criteria pass.

## Interpretation

The small projected-state dynamics error

```text
0.00029022869338069174
```

shows that the exact analytic LVC Gaussian coefficient propagation is highly accurate
for the state actually represented by the nine-Gaussian bank.

The larger error against the original target

```text
0.03500028070905269
```

is dominated by finite initial-state representation.

## v0.11 comparison

v0.11 remains stronger on one scalar diagonal-population error, but its center-frozen
full-density/coherence comparison is much poorer:

```text
v0.11 population error:              0.012877374121210683
v0.12 population error:              0.02810899300694737

v0.11 full-density error:            0.15991833275047374
v0.12 full-density error:            0.03500028070905269

v0.11 coherence phase error / rad:   1.367544547628621
v0.12 coherence phase error / rad:   0.0019607485027196615
```

v0.12 therefore advances the validation from population matching to
representation-consistent density-matrix/coherence matching.

## Terminology

The spinor-complete LVC path is a classically guided vector-Gaussian benchmark model.

The Born-Huang grid path is a projected benchmark/reference implementation.

Neither should be called a production molecular AIMS implementation.


## Automated validation

```text
149 passed
```
