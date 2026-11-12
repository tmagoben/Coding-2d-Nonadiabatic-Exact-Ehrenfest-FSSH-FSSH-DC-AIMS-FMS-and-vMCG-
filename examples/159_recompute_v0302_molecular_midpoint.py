"""Recompute v0.30.2 midpoint convergence, reversal and replay evidence."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gaussian_dynamics.v302_benchmark import run_v302_midpoint_campaign

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/v0302_molecular_midpoint_evidence.json")
    args = parser.parse_args()
    result = run_v302_midpoint_campaign()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({k: result[k] for k in ("passed", "gates", "observed_time_order", "refinement")}, indent=2))
    raise SystemExit(0 if result["passed"] else 1)
