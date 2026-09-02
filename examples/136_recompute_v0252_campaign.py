"""Recompute the cumulative v0.25.2 adaptive TDVP campaign."""

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

from gaussian_dynamics import run_v0252_release_benchmark
from gaussian_dynamics.campaign_io import save_campaign_json


output = run_v0252_release_benchmark()
campaign_path = save_campaign_json(
    root / "results/v0252_adaptive_multigaussian_tdvp_campaign.json",
    output,
)
evidence_path = save_campaign_json(
    root / "results/v0252_adaptive_multigaussian_tdvp_evidence.json",
    output["adaptive_multigaussian_validation_evidence"],
)
print("Saved campaign:", campaign_path)
print("Saved adaptive multi-Gaussian TDVP evidence:", evidence_path)
print("Acceptance:", output["acceptance"])
print("Claims:", output["claims"])

