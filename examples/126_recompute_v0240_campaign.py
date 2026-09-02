"""Recompute the canonical v0.24.0 external-SOC intake campaign."""

import os


for _thread_key in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_thread_key] = "1"

from pathlib import Path

from gaussian_dynamics import run_v0240_release_benchmark
from gaussian_dynamics.campaign_io import save_campaign_json


root = Path(__file__).resolve().parents[1]
output = run_v0240_release_benchmark()
campaign_path = save_campaign_json(
    root / "results/v0240_external_soc_intake_campaign.json", output
)
print("Saved campaign:", campaign_path)
print("Acceptance:", output["acceptance"])
print("Claims:", output["claims"])
