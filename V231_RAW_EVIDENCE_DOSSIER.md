# v0.23.1 raw-evidence admission dossier

## Purpose

v0.23.0 defined which evidence a molecular SOC backend needs. v0.23.1 records where
that evidence came from and recomputes every numerical summary. A passing number cannot
be inserted without changing the fingerprinted dossier, and real admission additionally
requires an executable backend-specific artifact parser.

## Bundle structure

Each canonical bundle contains:

```text
molecular_soc_admission_dossier_v231.json
replay/
  molecular_soc_manifest_v230.json
  molecular_soc_arrays_v230.npz
raw/
  calculation_template.json
  environment_lock.json
  independent_reference.json
  receipts/*.input.json
  receipts/*.output.json
```

The canonical fixtures each contain 17 receipts: nine exact trajectory records, three
basis-ladder calculations, two method-ladder calculations, and three frame calculations.
Their 37 raw artifacts include one template, one environment lock, one independent
reference, and distinct input/output artifacts for every receipt.

## Receipt contract

A receipt freezes:

- a unique record ID and role;
- the exact generalized coordinate in bohr;
- backend name/version and source classification;
- electronic method, basis, SOC operator, and derivative method;
- distinct input and output artifact names;
- SCF, correlated-state, SOC, derivative, and overlap convergence.

Every replay record must have exactly one trajectory receipt at the identical
coordinate. Receipt inputs and outputs cannot be reused across calculations. Basis and
method labels must agree with the corresponding raw convergence ladder.

## Derived evidence

For stored observable arrays $x$ and $x_{\rm ref}$, the reference error is computed as
either

$$
\epsilon_{\rm ref}=\max_i|x_i-x_{{\rm ref},i}|
$$

or the root-mean-square difference. A basis or method ladder derives one change between
every adjacent pair:

$$
\Delta_k=\operatorname{metric}(x_{k+1},x_k).
$$

Convergence is determined by the final adjacent change, not by a stored Boolean.
Translation and rotation residuals are recomputed from the base, transformed, and
expected-transformed observations. The rotation matrix must be finite, orthogonal, and
proper with determinant $+1$.

## Degenerate-manifold tracking

Statewise overlaps are not stable inside a degenerate multiplet. The tracking dossier
therefore partitions the full state space into named physical manifolds. For each edge
and manifold $g$, it derives

$$
s_g=\sigma_{\min}(O_{gg}),
$$

and the largest leakage to a competing manifold,

$$
c_g=\max_{h\ne g}\|O_{gh}\|_2.
$$

The assignment margin is $m_g=s_g-c_g$. Admission uses the minimum $s_g$ and $m_g$
over a connected record graph. The manifold groups must be disjoint and partition all
electronic states.

## Integrity and trust boundary

The dossier stores a relative path, role, byte size, and SHA-256 for every raw artifact.
Absolute paths and `..` traversal are forbidden. The dossier itself has a canonical
JSON fingerprint and is bound to the replay dataset fingerprint and provider identity.

Hashes establish byte identity, not scientific meaning. External or live admission
therefore also requires a method-specific validator whose name/version match the runtime
attestation and whose raw-artifact parser executes successfully. A manually constructed
attestation is insufficient.
