"""Recompute final admission, gauge audit and replay-verified dataset export."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gaussian_dynamics.molecular_dataset_v305 import admission_dataset_jsonl_v305
from gaussian_dynamics.v305_benchmark import run_v305_admission_campaign

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/v0305_molecular_admission_evidence.json")
    parser.add_argument("--generation-commit", default=None)
    args = parser.parse_args()
    if args.generation_commit is None:
        paths = ["gaussian_dynamics", "examples", "pyproject.toml"]
        dirty = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", *paths], cwd=ROOT).returncode
        untracked = subprocess.check_output(["git", "ls-files", "--others", "--exclude-standard", "--", *paths], cwd=ROOT, text=True)
        if dirty or untracked:
            parser.error("commit the implementation before recording generation provenance; use --generation-commit only for an externally verified source revision")
    commit = args.generation_commit or subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    evidence = run_v305_admission_campaign(generation_commit=commit, verbose=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(evidence, indent=2, sort_keys=True, allow_nan=False) + "\n")
    args.output.with_suffix(".jsonl").write_text(admission_dataset_jsonl_v305(evidence["dataset"]))
    print(json.dumps({"passed": evidence["passed"], "gates": evidence["gates"]}, indent=2))
    raise SystemExit(0 if evidence["passed"] else 1)
