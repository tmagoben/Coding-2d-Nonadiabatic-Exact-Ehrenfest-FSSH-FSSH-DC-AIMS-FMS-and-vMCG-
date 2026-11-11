# v0.29.3 Build Validation

## Classification

Development checkpoint only. Package metadata remains tied to the certified v0.27.0 release. The retained v0.29.2 development layer is preserved and extended; no public-release certification is claimed.

## Focused test receipts

- v0.29.3 new geometry/adaptation tests: **10/10 passed**.
- v0.29.3 + v0.29.2 + v0.27 correlated core replay: **37/37 passed**.
- additional v0.27 correlated validation replay: **4/4 passed**.
- v0.29.2 evidence script recomputed successfully with its prior scientific values.

## v0.29.3 deterministic evidence

Fingerprint:

`3d70cf7a13515e376df639b089d84365f382d3ea1116626308b4753a4d3f3a7c`

Selected results:

| quantity | value |
|---|---:|
| constant-connection reduction max abs | `0.0` |
| Wilson gauge covariance error | `0.0` |
| straight-vs-dogleg transporter separation | `1.3398687483e-03` |
| reverse-path inverse error | `2.5064762736e-16` |
| small-plaquette curvature relative error | `4.5429329599e-04` |
| patch local connection relation error | `0.0` |
| same-path two-patch transport error | `1.7608938669e-16` |
| spawn gain prediction relative error | `1.0988915010e-15` |
| spawn projection relative loss | `1.7916384190e-31` |
| spawn fidelity | `1.0` |
| newborn coefficient norm | `1.7735003729e-16` |
| newborn shape dormant | `True` |

## Interpretation

The checkpoint demonstrates that the path/patch geometry is internally consistent and that the v0.29.2 residual-driven spawning oracle survives replacement of a constant connection by a genuinely coordinate-dependent non-Abelian field.

The nonzero straight-vs-dogleg separation is an expected positive result: when curvature is nonzero, two geometrically different paths need not produce the same Wilson transporter.

## Honest boundaries

This is not a complete historical release replay. It does not certify a live molecular NAC/SOC trajectory under coordinate-dependent patch geometry, nor automatic atlas/path optimization.
