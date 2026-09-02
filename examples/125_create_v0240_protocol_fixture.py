"""Create a conspicuously synthetic v0.24.0 intake fixture."""

from pathlib import Path

from gaussian_dynamics import build_v240_protocol_fixture


root = Path(__file__).resolve().parents[1]
target = root / "results/v0240_openmolcas_protocol_fixture"
if target.exists():
    raise SystemExit(f"Refusing to overwrite existing fixture: {target}")
print("Created:", build_v240_protocol_fixture(target))
print("This is a protocol fixture, not OpenMolcas or ab-initio evidence.")

