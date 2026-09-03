# v0.21.4 Algorithmic Complexity

Let N be the Gaussian count, s the electronic dimension, d the nuclear dimension, E the
active sparse-edge count, and $C_{\mathrm{provider}}$ the cost of one provider snapshot.

## Differential provider audit

One audit geometry requires $2d+1$ provider evaluations and $2d$ cross-geometry
overlaps. The non-provider algebra performs two polar factorizations and Hamiltonian
transports per coordinate:

$$
O\!\left(d\,s^3\right).
$$

The total cost is therefore approximately

$$
O\!\left((2d+1)C_{\mathrm{provider}}+d\,C_{\mathrm{overlap}}+d\,s^3\right).
$$

The audit is a certification operation, not a per-timestep requirement.

## Checkpoint storage and I/O

The numerical state uses

$$
O\!\left(Nd^2+Ns+Ns^2+E\right)
$$

scalars: widths dominate the nuclear state for general d, guide densities contribute
$Ns^2$, and active graph state contributes E UID pairs. Digest construction, validation,
and compressed I/O are linear in the serialized byte count, aside from validating N
positive-definite width matrices at $O(Nd^3)$.

## Resume overhead

Restart reconstructs and validates N provider snapshots, restores $O(Ns^2)$ guidance
state and $O(E)$ sparse state, then resumes the unchanged propagation algorithm. Thus
the one-time non-provider restart overhead is

$$
O\!\left(Nd^3+Ns^2+E\right),
$$

plus $NC_{\mathrm{provider}}$. Checkpointing does not alter the asymptotic cost of subsequent
dense or sparse propagation steps.

## Zero-SOC rehearsal

Composing explicit zero H_SOC and K_SOC arrays copies or adds $O(ds^2)$ electronic data
per provider point. It does not change Gaussian graph scaling and is used for release
certification rather than to claim physical spin dynamics.

