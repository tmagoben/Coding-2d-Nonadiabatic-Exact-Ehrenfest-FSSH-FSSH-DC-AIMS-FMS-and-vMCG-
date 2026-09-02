from dataclasses import dataclass

from .sparse_complexity_v16 import (
    SparseComplexityLedgerV16,
    sparse_complexity_model_v16,
)


@dataclass
class SparseErrorComplexityLedgerV17(SparseComplexityLedgerV16):
    """v0.16 sparse ledger plus explicit error-control/audit work."""

    dense_audits: int = 0
    audit_pair_factorizations: int = 0
    score_relaxations: int = 0
    search_floor_relaxations: int = 0
    audit_seconds: float = 0.0

    max_audited_S_error: float = 0.0
    max_audited_H_error: float = 0.0
    max_audited_Snuc_error: float = 0.0

    def record_audit(self,audit):
        self.dense_audits+=1
        self.audit_pair_factorizations+=int(
            audit.get("dense_pair_factorizations",0)
        )
        self.max_audited_S_error=max(
            self.max_audited_S_error,
            float(audit["relative_S_frobenius_error"]),
        )
        self.max_audited_H_error=max(
            self.max_audited_H_error,
            float(audit["relative_H_frobenius_error"]),
        )
        self.max_audited_Snuc_error=max(
            self.max_audited_Snuc_error,
            float(audit["relative_Snuc_frobenius_error"]),
        )

    def as_dict(self):
        out=super().as_dict()
        out.update({
            "dense_audits":int(self.dense_audits),
            "audit_pair_factorizations":
                int(self.audit_pair_factorizations),
            "score_relaxations":
                int(self.score_relaxations),
            "search_floor_relaxations":
                int(self.search_floor_relaxations),
            "audit_seconds":float(self.audit_seconds),
            "max_audited_S_error":
                float(self.max_audited_S_error),
            "max_audited_H_error":
                float(self.max_audited_H_error),
            "max_audited_Snuc_error":
                float(self.max_audited_Snuc_error),
            "v17_complexity_note":(
                "Online dense audits intentionally reintroduce periodic O(N^2 d^3) "
                "calibration work. They are a correctness bridge for v0.17, not the "
                "final large-N sparse strategy."
            ),
            "asymptotic":{
                **sparse_complexity_model_v16().__dict__,
                "periodic_dense_audit":(
                    "O(N^2 d^3 + s^2 N^2) per audit checkpoint; amortized by the "
                    "chosen audit interval"
                ),
                "edge_importance":(
                    "O(M d^3) exact S/H/T scoring on KD-tree candidate pairs M"
                ),
            },
        })
        return out
