"""Recompute unseeded molecular packet activation and adaptive propagation."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gaussian_dynamics.v306_benchmark import run_v306_adaptive_campaign

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/v0306_adaptive_propagation_evidence.json")
    args = parser.parse_args()
    paths = ["gaussian_dynamics", "examples", "pyproject.toml"]
    dirty = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", *paths], cwd=ROOT).returncode
    untracked = subprocess.check_output(["git", "ls-files", "--others", "--exclude-standard", "--", *paths], cwd=ROOT, text=True)
    if dirty or untracked:
        parser.error("commit the implementation before recording generation provenance")
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    evidence = run_v306_adaptive_campaign(generation_commit=commit, verbose=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(evidence, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"passed": evidence["passed"], "gates": evidence["gates"], "measurements": evidence["measurements"]}, indent=2))
    raise SystemExit(0 if evidence["passed"] else 1)
