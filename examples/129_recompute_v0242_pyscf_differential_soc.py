"""Recompute the pinned OH connected-geometry PySCF SOC evidence."""

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
    run_pyscf_oh_bond_differential_evidence_v242,
    save_pyscf_oh_bond_differential_evidence_v242,
)


evidence = run_pyscf_oh_bond_differential_evidence_v242()
path = save_pyscf_oh_bond_differential_evidence_v242(
    root / "results/v0242_pyscf_differential_soc_evidence.json",
    evidence,
)
print("Saved evidence:", path)
print("Evidence fingerprint:", evidence.fingerprint())
print("Acceptance:", evidence.audit.passed)
print("Runtime gates:", len(evidence.audit.checks))
print("Metrics:", evidence.audit.metrics)
print("Claims:", evidence.claims)
