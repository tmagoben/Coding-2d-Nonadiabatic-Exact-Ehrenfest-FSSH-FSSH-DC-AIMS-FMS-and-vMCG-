"""Recompute molecular residual-oracle scores, gauge comparisons and exact replay."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gaussian_dynamics.v304_benchmark import run_v304_spawn_campaign

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/v0304_molecular_spawning_evidence.json")
    args = parser.parse_args()
    evidence = run_v304_spawn_campaign(verbose=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(evidence, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"passed": evidence["passed"], "gates": evidence["gates"]}, indent=2))
    raise SystemExit(0 if evidence["passed"] else 1)
