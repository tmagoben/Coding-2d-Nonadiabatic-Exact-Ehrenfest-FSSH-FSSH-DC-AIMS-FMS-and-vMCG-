from dataclasses import dataclass

from .sparse_error_complexity_v17 import (
    SparseErrorComplexityLedgerV17,
)


@dataclass
class ConvergenceComplexityLedgerV18(SparseErrorComplexityLedgerV17):
    """Complexity ledger for sampled-audit + batched-candidate v0.18."""

    sampled_audits: int = 0
    sampled_pairs_scored: int = 0
    sampled_audit_failures: int = 0
    sampled_audit_seconds: float = 0.0

    sentinel_dense_audits: int = 0
    sentinel_pair_factorizations: int = 0

    candidate_batches: int = 0
    candidate_dense_grid_elements: int = 0
    candidate_max_dense_grid_elements: int = 0
    candidate_peak_grid_elements: int = 0
    candidate_grid_elements_processed: int = 0

    wavefunction_metric_evaluations: int = 0
    wavefunction_metric_seconds: float = 0.0

    def record_sampled_audit(self,audit):
        self.sampled_audits+=1
        self.sampled_pairs_scored+=len(
            audit.sampled_pairs
        )
        if not audit.passed:
            self.sampled_audit_failures+=1

    def record_sentinel_audit(self,audit):
        self.sentinel_dense_audits+=1
        pairs=int(
            audit.get("dense_pair_factorizations",0)
        )
        self.sentinel_pair_factorizations+=pairs
        # Retain inherited dense-audit accounting for total correctness cost.
        self.record_audit(audit)

    def record_candidate_batching(self,diag):
        self.candidate_batches+=int(diag.batches)
        self.candidate_dense_grid_elements+=int(
            diag.dense_candidate_grid_elements
        )
        self.candidate_max_dense_grid_elements=max(
            self.candidate_max_dense_grid_elements,
            int(diag.dense_candidate_grid_elements),
        )
        self.candidate_peak_grid_elements=max(
            self.candidate_peak_grid_elements,
            int(diag.peak_candidate_grid_elements),
        )
        # Each batch is generated exactly once, so processed elements equal K*G even
        # though peak storage is bounded by B*G.
        self.candidate_grid_elements_processed+=int(
            diag.dense_candidate_grid_elements
        )

    @property
    def candidate_peak_memory_reduction_fraction(self):
        return float(
            1.0
            -self.candidate_peak_grid_elements
            /max(self.candidate_max_dense_grid_elements,1)
        )

    def as_dict(self):
        out=super().as_dict()
        out.update({
            "sampled_audits":int(self.sampled_audits),
            "sampled_pairs_scored":
                int(self.sampled_pairs_scored),
            "sampled_audit_failures":
                int(self.sampled_audit_failures),
            "sampled_audit_seconds":
                float(self.sampled_audit_seconds),
            "sentinel_dense_audits":
                int(self.sentinel_dense_audits),
            "sentinel_pair_factorizations":
                int(self.sentinel_pair_factorizations),
            "candidate_batches":
                int(self.candidate_batches),
            "candidate_dense_grid_elements":
                int(self.candidate_dense_grid_elements),
            "candidate_max_dense_grid_elements":
                int(self.candidate_max_dense_grid_elements),
            "candidate_peak_grid_elements":
                int(self.candidate_peak_grid_elements),
            "candidate_grid_elements_processed":
                int(self.candidate_grid_elements_processed),
            "candidate_peak_memory_reduction_fraction":
                self.candidate_peak_memory_reduction_fraction,
            "wavefunction_metric_evaluations":
                int(self.wavefunction_metric_evaluations),
            "wavefunction_metric_seconds":
                float(self.wavefunction_metric_seconds),
            "v18_complexity_note":(
                "Normal sparse audits are sampled. Full O(N^2) dense matrix audits are "
                "reserved for initial/final sentinels and small validation campaigns."
            ),
        })
        out["asymptotic"].update({
            "sampled_sparse_audit":(
                "O(J d^3) exact omitted-edge scoring for fixed sample size J, plus "
                "local KD-tree priority search; no full dense matrix rebuild"
            ),
            "batched_candidate_memory":(
                "O((N+B)G) complex grid storage instead of O((N+K)G), with batch "
                "size B << candidate count K"
            ),
            "sentinel_dense_audit":(
                "O(N^2 d^3 + s^2 N^2), used only at initial/final release sentinels"
            ),
        })
        return out
