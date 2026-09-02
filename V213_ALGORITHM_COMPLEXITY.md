# v0.21.3 Algorithmic Complexity

Let $N$ be the Gaussian count, $s$ the electronic model-space dimension, $d$ the
nuclear dimension, $E$ the active Gaussian-edge count, and $K$ the number of nuclear
grid points used for initialization.

## Strict invariants

Checking H costs $O(s^2)$. Checking all physical derivative and connection matrices
costs $O(ds^2)$. A frame Gram matrix costs $O(rs^2)$ for a representation dimension
$r$ and reduces to $O(s^3)$ for a square frame. These checks do not change the
propagation asymptotics.

## Density-matrix guidance

Evaluating

$$
F_{i,a}=-\operatorname{Tr}(\rho_iK_a)
$$

for one Gaussian costs $O(ds^2)$. All local forces cost

$$
O(Nds^2)
$$

after provider evaluation. Transport requires one nearest-unitary factorization per
moved guide density, $O(s^3)$ each, for a worst-case step cost

$$
O\left[N(ds^2+s^3)\right].
$$

Guide-density storage is $O(Ns^2)$. This replaces an $O(s^3)$ eigenvector selection
whose result was not physically defined at degeneracy.

## Model-space and provenance validation

State and multiplet validation is $O(s)$ apart from metadata serialization. Fingerprint
generation is linear in the canonical provenance payload size and is normally
negligible relative to electronic-structure evaluation.

## Grid projection

Building all analytic nuclear overlaps costs approximately

$$
O(N^2d^3)
$$

for general dense width matrices. Evaluating all Gaussian/grid right-hand sides and
reconstructing the wavefunction costs

$$
O(KNs).
$$

The current robust dense solve on the $Ns$ block metric is formally

$$
O((Ns)^3)
$$

with $O((Ns)^2)$ storage. The fixed-frame Kronecker structure can support a later
nuclear-only factorization, but v0.21.3 prioritizes a transparent reference
implementation.

## Fixed-frame complex cache

Each entry stores $O(ds^2+s^2+d^2)$ numerical data plus metadata. Hashing and JSON
work are linear in coordinate/provenance size; compressed array I/O is linear in entry
size. Cache hits avoid the provider's electronic-structure cost.

## Impact of later SOC

SOC changes the content of H and K and may increase $s$ by expanding full multiplets.
It does not change the Gaussian graph or propagation formulas. The important scaling
effect is therefore the existing block dependence on $s$, especially $s^2$ storage and
$s^3$ dense electronic factorizations/products. Spin-free mode retains its original
model-space size.
