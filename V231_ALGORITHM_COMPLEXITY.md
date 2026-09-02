# v0.23.1 algorithmic complexity

Let $R$ be the number of replay records, $s$ the electronic-state dimension, $d$ the
nuclear-coordinate dimension, $E_t$ the number of tracking edges, $G$ the number of
physical manifolds, and $L_b,L_m$ the basis and method ladder lengths.

## Receipt and artifact validation

The dossier contains $R+L_b+L_m+3$ calculation receipts and two raw artifacts per
receipt, plus fixed template, environment, reference, and optional runtime artifacts.
Structural receipt validation is linear in the number of receipts. Hash verification is

$$
O(B_{\mathrm{raw}}),
$$

where $B_{\mathrm{raw}}$ is the total raw-artifact byte count. Memory use is linear in the
manifest and one hashing block; artifacts are streamed rather than loaded together.

## Evidence derivation

Reference and frame metrics are linear in their stored observable sizes. Basis and
method derivation costs

$$
O\!\left((L_b+L_m)N_o\right),
$$

for $N_o$ stored observable components. No convergence Boolean is trusted.

## Manifold tracking

For each tracking edge and manifold, the assigned overlap block is decomposed by SVD
and competing blocks receive spectral norms. With manifold sizes $s_g$, the dominant
dense cost is approximately

$$
O\!\left(E_t\sum_{g=1}^{G}s_g^3
+E_t\sum_{g\ne h}\min(s_g,s_h)^2\max(s_g,s_h)\right).
$$

The canonical connected path uses $E_t=R-1$. The full v0.23.0 replay still stores all
$R^2$ overlap blocks, so its $O(R^2s^2)$ storage and $O(R^2s^3)$ integrity validation
remain the larger asymptotic costs.

## Admission

Canonical JSON fingerprinting is linear in dossier size. The inherited physical audit
retains its finite-difference and dense-matrix costs. A real backend-specific parser is
method dependent and can dominate all framework-side work, but it must execute only at
admission or replay creation—not during ordinary file-backed propagation.
