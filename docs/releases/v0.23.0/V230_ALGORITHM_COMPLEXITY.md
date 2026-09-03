# v0.23.0 algorithmic complexity

Let $R$ be the number of replay records, $d$ the nuclear-coordinate dimension, $s$ the
electronic-state dimension, and $P$ the number of physical projectors.

## Capture

Capturing point operators requires $R$ electronic evaluations. Capturing all ordered
cross-record overlaps requires $R^2$ overlap evaluations. Excluding backend cost, array
assembly uses

$$
O\!\left(Rd s^2 + R^2s^2 + Rd^2 + Ps^2\right)
$$

memory and comparable serialization work. The $R^2s^2$ overlap table dominates large
replays. This quadratic design is intentional in v0.23.0: every differential and
tracking relation is reproducible without calling the source backend.

## Load and validation

Array hashing and finite-value scans are linear in stored bytes. Checking all overlap
isometries performs $R^2$ dense $s\times s$ matrix products, costing

$$
O(R^2s^3)
$$

time and $O(R^2s^2)$ memory. Structural checks on the component derivatives cost
$O(Rds^2)$ apart from the inherited matrix audits.

## Exact lookup

The file-backed provider constructs a hash map of rounded coordinate tuples in
$O(Rd)$ time and memory. An exact-coordinate lookup is $O(d)$ on average, after which
returning an electronic record copies $O(ds^2+d^2)$ values. No interpolation system is
built.

## Admission audit

For a center geometry and $h$ finite-difference scales, the inherited component audit
requests $O(dh)$ neighboring records. Dense matrix residuals cost $O(dhs^2)$; overlap
isometry checks cost $O(dhs^3)$. The evidence gates themselves are linear in the basis
and method ladder lengths.

## Production implication

For a future large external dataset, the all-pairs overlap table is the first scaling
limit. A later format may store a fingerprinted sparse overlap graph, but that change
must preserve the exact tracking paths required by every admitted trajectory and must
introduce a new format version. v0.23.0 does not silently sparsify replay evidence.
