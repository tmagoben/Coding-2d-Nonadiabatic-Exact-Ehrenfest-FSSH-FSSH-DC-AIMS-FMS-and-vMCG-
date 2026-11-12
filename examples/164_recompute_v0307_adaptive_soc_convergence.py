"""Recompute adaptive SOC convergence with fixed physical adaptation checkpoints."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from gaussian_dynamics.v307_benchmark import run_v307_convergence_campaign

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT / "results/v0307_adaptive_soc_convergence_evidence.json")
    parser.add_argument("--workers",type=int,default=2,help="independent campaign workers; use 1 for serial execution")
    parser.add_argument("--checkpoint-directory",type=Path,help="optional source-bound cache of completed independent cases")
    args = parser.parse_args()
    paths = ["gaussian_dynamics","examples","pyproject.toml"]
    dirty = subprocess.run(["git","diff","--quiet","HEAD","--",*paths],cwd=ROOT).returncode
    untracked = subprocess.check_output(["git","ls-files","--others","--exclude-standard","--",*paths],cwd=ROOT,text=True)
    if dirty or untracked:
        parser.error("commit the implementation before recording generation provenance")
    commit = subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip()
    evidence = run_v307_convergence_campaign(generation_commit=commit,verbose=True,workers=args.workers,checkpoint_directory=args.checkpoint_directory)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(evidence,sort_keys=True,separators=(",",":"),allow_nan=False)+"\n")
    print(json.dumps({"passed":evidence["passed"],"gates":evidence["gates"],"measurements":evidence["measurements"]},indent=2))
    raise SystemExit(0 if evidence["passed"] else 1)
