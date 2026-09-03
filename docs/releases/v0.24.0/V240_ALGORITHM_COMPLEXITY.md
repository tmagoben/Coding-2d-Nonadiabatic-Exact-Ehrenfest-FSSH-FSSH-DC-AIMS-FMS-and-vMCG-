# v0.24.0 algorithm complexity

Let `A` be the number of atoms, `Q=3A` the Cartesian dimension, `S` the number of
spin components, `L` the displacement ladder length, and `R=1+2QL` the record count.

- Artifact hashing is linear in total byte count.
- Strict parsing stores `O(R S^2 + R Q)` numeric data.
- Polar transport and complete-manifold SVDs cost `O(R S^3)`.
- Centered component derivatives cost `O(Q L S^3)` including transports.
- Static electronic propagation precomputes one `S x S` exponential in `O(S^3)` and
  advances each step in `O(S^2)`.

For the frozen H2O protocol, `A=3`, `Q=9`, `S=4`, `L=3`, and `R=55`. The campaign is
therefore dominated by the inherited real PySCF validation, not these small matrix
audits. Production OpenMolcas wall time and storage are backend costs and are not
estimated by this framework complexity statement.

