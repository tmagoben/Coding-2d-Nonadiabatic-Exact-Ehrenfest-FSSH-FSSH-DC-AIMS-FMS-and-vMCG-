"""Recompute the canonical v0.23.2 real-PySCF release campaign."""

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

from gaussian_dynamics import run_v0232_release_benchmark
from gaussian_dynamics.campaign_io import save_campaign_json


root = Path(__file__).resolve().parents[1]
output = run_v0232_release_benchmark()
campaign_path = save_campaign_json(
    root / "results/v0232_pyscf_runtime_admission_campaign.json", output
)
runtime_path = save_campaign_json(
    root / "results/v0232_pyscf_runtime_evidence.json",
    output["pyscf_runtime_evidence"],
)
print("Saved campaign:", campaign_path)
print("Saved runtime evidence:", runtime_path)
print("Runtime evidence SHA-256:", output["pyscf_runtime_evidence"]["evidence_sha256"])
print("Acceptance:", output["acceptance"])
print("Claims:", output["claims"])
