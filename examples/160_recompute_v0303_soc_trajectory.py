"""Recompute fixed/moving-frame SOC trajectory convergence and controls."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gaussian_dynamics.v303_benchmark import run_v303_soc_campaign

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/v0303_soc_trajectory_evidence.json")
    parser.add_argument("--workers", type=int, choices=(1, 2, 3), default=2)
    args = parser.parse_args()
    evidence = run_v303_soc_campaign(workers=args.workers, verbose=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(evidence, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"passed": evidence["passed"], "gates": evidence["gates"],
        "time_orders": {key: value["observed_time_order"] for key, value in evidence["frames"].items()},
        "gauge_wavefunction_errors": evidence["gauge_wavefunction_errors"],
        "soc_disabled": evidence["soc_disabled"]}, indent=2))
    raise SystemExit(0 if evidence["passed"] else 1)
