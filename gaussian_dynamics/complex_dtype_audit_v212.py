from dataclasses import dataclass
from pathlib import Path
import re


@dataclass(frozen=True)
class ComplexDtypeAuditResultV212:
    files_scanned: int
    suspicious_casts: tuple
    intentional_legacy_casts: tuple

    @property
    def passed(self):
        return len(self.suspicious_casts)==0

    def as_dict(self):
        return {
            "files_scanned":int(self.files_scanned),
            "suspicious_casts":list(self.suspicious_casts),
            "intentional_legacy_casts":list(self.intentional_legacy_casts),
            "passed":bool(self.passed),
        }


CORE_FILES=(
    "electronic_operator_v21.py",
    "complex_gauge_v21.py",
    "block_sparse_molecular_v21.py",
    "block_dynamics_v21.py",
    "synthetic_operator_provider_v21.py",
    "self_consistent_block_v212.py",
    "electronic_observables_v212.py",
    "subspace_provider_v212.py",
    "block_basis_lifecycle_v212.py",
)

# This cast is intentional: the adapter's *source* contract is the legacy real,
# spin-free adiabatic provider.  The converted v0.21 operator matrices are complex.
INTENTIONAL=(
    ("electronic_operator_v21.py","nac = np.asarray(point.nac_q, dtype=float)"),
    ("electronic_operator_v21.py","E = np.asarray(point.energies, dtype=float)"),
    ("electronic_operator_v21.py","grad = np.asarray(point.gradients_q, dtype=float)"),
)

RISK_WORDS=(
    "H", "dH", "connection", "state_vectors", "coeff", "overlap",
    "observable", "nac", "wavefunction",
)


def audit_pre_soc_complex_core_v212(package_dir):
    package_dir=Path(package_dir)
    suspicious=[]; intentional=[]
    for name in CORE_FILES:
        path=package_dir/name
        text=path.read_text(encoding="utf-8")
        for lineno,line in enumerate(text.splitlines(),start=1):
            compact=line.strip()
            if not (
                "dtype=float" in compact
                or re.search(r"np\.asarray\([^\n]*,\s*float\)",compact)
            ):
                continue
            if any(file==name and snippet in compact for file,snippet in INTENTIONAL):
                intentional.append(f"{name}:{lineno}: {compact}")
                continue
            if any(word in compact for word in RISK_WORDS):
                # q, p, masses, forces and scalar diagnostics are physically real.
                real_context=(" q" in compact or "q)" in compact or "q," in compact or
                              "mass" in compact.lower() or "force" in compact.lower() or
                              "qdots" in compact or "pdots" in compact or "veloc" in compact)
                if not real_context:
                    suspicious.append(f"{name}:{lineno}: {compact}")
    return ComplexDtypeAuditResultV212(
        files_scanned=len(CORE_FILES),
        suspicious_casts=tuple(suspicious),
        intentional_legacy_casts=tuple(intentional),
    )
