"""Recompute the deterministic v0.26.0 multidimensional validation evidence."""

import os
import sys
from pathlib import Path


for _thread_key in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_thread_key] = "1"

root = Path(__file__).resolve().parents[1]
if str(root) not in sys.path:
    sys.path.insert(0, str(root))

from gaussian_dynamics import save_multidimensional_validation_evidence_v260


if __name__ == "__main__":
    target = root / "results/v0260_multidimensional_evidence.json"
    evidence = save_multidimensional_validation_evidence_v260(target)
    print(f"saved {target}")
    print(f"passed {evidence.check_count}/{evidence.check_count} gates")
    print(f"fingerprint {evidence.fingerprint()}")
