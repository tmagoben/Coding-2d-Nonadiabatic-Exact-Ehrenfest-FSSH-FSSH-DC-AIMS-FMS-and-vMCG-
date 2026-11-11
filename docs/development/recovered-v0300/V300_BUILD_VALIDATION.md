# v0.30.0 Build Validation

## Baseline audit

The recovered v0.29.4 checkpoint compiled successfully. Its unrestricted inherited
suite produced 584 passes and 16 optional-runtime skips; 11 cumulative release tests
failed only because PySCF was initially absent. After installing the pinned PySCF
2.13.1 runtime and the locked NumPy/SciPy versions, all directly targeted molecular,
differential, v0.29.2, v0.29.3, and v0.29.4 tests passed.

The only remaining cumulative red gate was
`release_locked_profile_passes_canonical_runtime`. This is an intentional
byte-identical-machine gate: the audit runtime is scientifically compatible but is not
the exact Python executable used to seal v0.27.0. No physics or numerical assertion
failed, so changing the release identity would have been incorrect.

The full portable scientific suite, excluding only the ten downstream cumulative
meta-tests that require that byte-identical executable, passed **611/611 tests** in the
locked PySCF/NumPy/SciPy environment. The underlying v0.23.2 real-PySCF release gate
was included in that run.

## v0.30.0 additions

- `gaussian_dynamics/molecular_frame_bridge_v300.py`
- `gaussian_dynamics/v300_benchmark.py`
- `tests/test_v0300_molecular_frame_bridge.py`
- `results/v0300_molecular_frame_bridge_evidence.json`
- `V300_MOLECULAR_FRAME_BRIDGE.md`
- `V300_PROGRAM_ARCHITECTURE.md`

## Focused acceptance

- v0.30.0 Python tests: **10/10 passed**;
- combined inherited v0.27 correlated-width, v0.29.2-v0.29.4 curved/pre-ML,
  and v0.30.0 integration suite: **58/58 passed**;
- deterministic scientific/data-contract gates: **12/12 passed**;
- maximum NAC-Wilson versus overlap-polar residual:
  `3.70e-16` (well below `3e-9`);
- maximum Wilson unitarity residual: `<4.0e-16`;
- maximum overlap-polar unitarity residual: `<6.7e-16`;
- Hamiltonian, physical derivative, and connection covariance residuals: `0.0`
  for the analytic moving-frame SOC fixture.

## Adversarial controls

The focused suite confirms that spectrally expansive raw overlaps,
non-anti-Hermitian derivative connections, missing production numerical identity,
state/path fingerprint mismatches, and attempts to select an unspecified transport
mode all fail closed.

## Claim boundary

v0.30.0 certifies the molecular-frame **interface and local transport agreement**, not
a full molecular trajectory or ab-initio predictive result.
