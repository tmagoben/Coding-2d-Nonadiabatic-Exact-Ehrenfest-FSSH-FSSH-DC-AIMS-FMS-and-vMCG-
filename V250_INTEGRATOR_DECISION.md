# v0.25.0 integrator decision

## Decision summary

| Question | v0.25.0 decision | Boundary |
|---|---|---|
| Use a time-dependent variational procedure? | Yes: validate the single-packet, complete-spinor TDVP restriction first. | Full coupled multi-Gaussian TDVP remains closed. |
| Use velocity Verlet for nuclei? | Yes, for constant-mass canonical `(q,p)` only. | Reject coordinate-dependent mass and general noncanonical TDVP coordinates. |
| Polar decomposition or SVD? | Both: use the polar factor physically and compute/certify it by SVD. | Never propagate amplitudes with the raw contractive overlap. |
| What for full TDVP later? | Implicit midpoint/discrete variational residual solve. | Must add metric, constraint, null-space, and nonlinear-solver evidence. |

## Rationale

Velocity Verlet is symmetric, second order, inexpensive, and appropriate for the
released separable kinetic energy `p^T M^-1 p/2`. Its assumptions are explicit and
runtime checked. Calling the same update a solver for a general Gaussian TDVP metric
would overstate the mathematics.

The phrase “polar versus SVD” mixes a mathematical object with an algorithm. The
unitary polar factor is the required gauge/frame transport. SVD supplies that factor
robustly while also exposing the information that determines whether the retained
electronic manifold is trustworthy. v0.25.0 therefore freezes

$$O=U\Sigma V^\dagger,\qquad W=UV^\dagger.$$

No inverse singular values enter the transport itself; small singular values instead
close the trajectory gate.
