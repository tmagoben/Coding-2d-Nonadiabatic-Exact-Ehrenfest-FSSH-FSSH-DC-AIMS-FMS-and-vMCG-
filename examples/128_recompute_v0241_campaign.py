"""Recompute the canonical v0.24.1 PySCF static-SOC release campaign."""

import os
import sys
from pathlib import Path


for _thread_key in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_thread_key] = "1"

root = Path(__file__).resolve().parents[1]
if str(root) not in sys.path:
    sys.path.insert(0, str(root))

from gaussian_dynamics import run_v0241_release_benchmark
from gaussian_dynamics.campaign_io import save_campaign_json


output = run_v0241_release_benchmark()
campaign_path = save_campaign_json(
    root / "results/v0241_pyscf_static_soc_campaign.json", output
)
runtime_path = save_campaign_json(
    root / "results/v0241_pyscf_static_soc_evidence.json",
    output["pyscf_static_soc_runtime_evidence"],
)
print("Saved campaign:", campaign_path)
print("Saved runtime evidence:", runtime_path)
print("Acceptance:", output["acceptance"])
print("Claims:", output["claims"])
