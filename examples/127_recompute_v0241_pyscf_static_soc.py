"""Recompute the pinned OH PySCF BP-SOMF static-SOC evidence."""

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

from gaussian_dynamics import (
    run_pyscf_oh_static_soc_evidence_v241,
    save_pyscf_oh_static_soc_evidence_v241,
)


evidence = run_pyscf_oh_static_soc_evidence_v241()
path = save_pyscf_oh_static_soc_evidence_v241(
    root / "results/v0241_pyscf_static_soc_evidence.json", evidence
)
print("Saved evidence:", path)
print("Evidence fingerprint:", evidence.fingerprint())
print("Acceptance:", evidence.audit.passed)
print("Metrics:", evidence.audit.metrics)
print("Claims:", evidence.claims)
