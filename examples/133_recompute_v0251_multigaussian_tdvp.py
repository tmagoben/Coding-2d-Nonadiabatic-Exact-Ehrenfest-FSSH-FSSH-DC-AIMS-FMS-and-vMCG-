"""Recompute deterministic v0.25.1 multi-Gaussian TDVP evidence."""

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

from gaussian_dynamics import run_multigaussian_tdvp_validation_evidence_v251
from gaussian_dynamics.campaign_io import save_campaign_json


evidence = run_multigaussian_tdvp_validation_evidence_v251()
evidence_path = save_campaign_json(
    root / "results/v0251_multigaussian_tdvp_evidence.json",
    evidence.as_dict(),
)
print("Saved multi-Gaussian TDVP evidence:", evidence_path)
print("Fingerprint:", evidence.fingerprint())
print("Acceptance:", evidence.audit.as_dict())
print("Claims:", evidence.claims)

