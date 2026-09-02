from dataclasses import dataclass, field
from contextlib import contextmanager
import time


@dataclass(frozen=True)
class SparseComplexityModelV16:
    locality_screen: str
    active_pair_algebra: str
    sparse_matrix_storage: str
    sparse_cayley: str
    candidate_ranking: str
    electronic_cost: str


def sparse_complexity_model_v16():
    return SparseComplexityModelV16(
        locality_screen=(
            "O(N d^3 + N log N + M d) for width-eigenvalue preprocessing plus a "
            "safe cKDTree global-radius query and M pair-specific bound checks; "
            "worst case degrades to O(N^2 d) when the global radius spans the basis"
        ),
        active_pair_algebra=(
            "O(E d^3) pair algebra per snapshot, where E includes active canonical "
            "offdiagonal edges plus N diagonals"
        ),
        sparse_matrix_storage=(
            "O(s^2(N+2E_off)) nonzeros for block-sparse S/H/T before sparse-solver fill"
        ),
        sparse_cayley=(
            "sparse direct solve cost depends on graph fill/order; reported with nnz and "
            "wall time rather than falsely assigning a universal O(N^p) exponent"
        ),
        candidate_ranking=(
            "O(KG(N+s)+N^2K+N^3) residual shortlist plus O(K_short N d) local-degree "
            "and electronic-cost reranking"
        ),
        electronic_cost=(
            "provider-dependent; v0.16 records explicit cache-hit/new-geometry cost units "
            "instead of folding ab-initio cost into Gaussian linear algebra"
        ),
    )


@dataclass
class SparseComplexityLedgerV16:
    endpoint_graph_updates: int = 0
    midpoint_graph_updates: int = 0
    exact_pair_checks: int = 0
    screened_pairs: int = 0
    spatial_candidate_pairs: int = 0
    globally_screened_pairs: int = 0
    pair_factorizations: int = 0
    propagation_pair_factorizations: int = 0
    candidate_pair_factorizations: int = 0
    pair_requests: int = 0

    sparse_matrix_builds: int = 0
    sparse_time_builds: int = 0
    sparse_cayley_solves: int = 0
    sparse_defect_solves: int = 0

    candidate_searches: int = 0
    candidates_scored: int = 0
    cost_reranks: int = 0
    enrichments: int = 0
    pruning_events: int = 0

    electronic_cache_hits: int = 0
    electronic_cache_misses: int = 0
    electronic_cost_units: float = 0.0

    peak_basis_size: int = 0
    peak_active_edges: int = 0
    peak_S_nnz: int = 0
    peak_H_nnz: int = 0
    sum_active_edges: int = 0
    sum_total_edges: int = 0
    graph_samples: int = 0

    graph_seconds: float = 0.0
    matrix_seconds: float = 0.0
    time_matrix_seconds: float = 0.0
    cayley_seconds: float = 0.0
    defect_seconds: float = 0.0
    candidate_seconds: float = 0.0
    cost_seconds: float = 0.0
    pruning_seconds: float = 0.0
    total_seconds: float = 0.0

    _start: float = field(default=0.0,repr=False)

    def start(self):
        self._start=time.perf_counter()

    def stop(self):
        if self._start:
            self.total_seconds+=time.perf_counter()-self._start
            self._start=0.0

    @contextmanager
    def timed(self,name):
        t0=time.perf_counter()
        yield
        attr=f"{name}_seconds"
        if not hasattr(self,attr):
            raise ValueError(f"unknown timing category {name!r}")
        setattr(
            self,attr,
            getattr(self,attr)+time.perf_counter()-t0,
        )

    def observe_graph(self,update,n_basis):
        self.exact_pair_checks+=int(update.exact_pair_checks)
        self.screened_pairs+=int(update.screened_pairs)
        self.spatial_candidate_pairs+=int(
            update.spatial_candidate_pairs
        )
        self.globally_screened_pairs+=int(
            update.globally_screened_pairs
        )
        self.peak_basis_size=max(
            self.peak_basis_size,int(n_basis)
        )
        self.peak_active_edges=max(
            self.peak_active_edges,
            int(update.active_offdiagonal_edges),
        )
        self.sum_active_edges+=int(
            update.active_offdiagonal_edges
        )
        self.sum_total_edges+=int(
            update.total_offdiagonal_pairs
        )
        self.graph_samples+=1

    def record_pair_delta(self,cache,before,category="propagation"):
        delta_solves=int(
            cache.stats.canonical_solves
            -before.get("canonical_solves",0)
        )
        delta_requests=int(
            cache.stats.requests
            -before.get("requests",0)
        )
        self.pair_factorizations+=delta_solves
        self.pair_requests+=delta_requests
        if category=="propagation":
            self.propagation_pair_factorizations+=delta_solves
        elif category=="candidate":
            self.candidate_pair_factorizations+=delta_solves
        else:
            raise ValueError("category must be propagation or candidate.")

    def observe_matrices(self,mats):
        self.peak_S_nnz=max(
            self.peak_S_nnz,int(mats.S.nnz)
        )
        self.peak_H_nnz=max(
            self.peak_H_nnz,int(mats.H.nnz)
        )

    def record_electronic_cost(self,estimate):
        self.electronic_cost_units+=float(
            estimate.cost_units
        )
        if estimate.cache_hit:
            self.electronic_cache_hits+=1
        else:
            self.electronic_cache_misses+=1

    @property
    def average_edge_fraction(self):
        return float(
            self.sum_active_edges
            /max(self.sum_total_edges,1)
        )

    @property
    def average_sparsity_fraction(self):
        return float(1.0-self.average_edge_fraction)

    def as_dict(self):
        return {
            "endpoint_graph_updates":
                int(self.endpoint_graph_updates),
            "midpoint_graph_updates":
                int(self.midpoint_graph_updates),
            "exact_pair_checks":
                int(self.exact_pair_checks),
            "screened_pairs":
                int(self.screened_pairs),
            "spatial_candidate_pairs":
                int(self.spatial_candidate_pairs),
            "globally_screened_pairs":
                int(self.globally_screened_pairs),
            "pair_factorizations":
                int(self.pair_factorizations),
            "propagation_pair_factorizations":
                int(self.propagation_pair_factorizations),
            "candidate_pair_factorizations":
                int(self.candidate_pair_factorizations),
            "pair_requests":
                int(self.pair_requests),
            "sparse_matrix_builds":
                int(self.sparse_matrix_builds),
            "sparse_time_builds":
                int(self.sparse_time_builds),
            "sparse_cayley_solves":
                int(self.sparse_cayley_solves),
            "sparse_defect_solves":
                int(self.sparse_defect_solves),
            "candidate_searches":
                int(self.candidate_searches),
            "candidates_scored":
                int(self.candidates_scored),
            "cost_reranks":
                int(self.cost_reranks),
            "enrichments":int(self.enrichments),
            "pruning_events":int(self.pruning_events),
            "electronic_cache_hits":
                int(self.electronic_cache_hits),
            "electronic_cache_misses":
                int(self.electronic_cache_misses),
            "electronic_cost_units":
                float(self.electronic_cost_units),
            "peak_basis_size":
                int(self.peak_basis_size),
            "peak_active_edges":
                int(self.peak_active_edges),
            "peak_S_nnz":int(self.peak_S_nnz),
            "peak_H_nnz":int(self.peak_H_nnz),
            "average_edge_fraction":
                self.average_edge_fraction,
            "average_sparsity_fraction":
                self.average_sparsity_fraction,
            "graph_seconds":float(self.graph_seconds),
            "matrix_seconds":float(self.matrix_seconds),
            "time_matrix_seconds":
                float(self.time_matrix_seconds),
            "cayley_seconds":float(self.cayley_seconds),
            "defect_seconds":float(self.defect_seconds),
            "candidate_seconds":
                float(self.candidate_seconds),
            "cost_seconds":float(self.cost_seconds),
            "pruning_seconds":
                float(self.pruning_seconds),
            "total_seconds":float(self.total_seconds),
            "asymptotic":
                sparse_complexity_model_v16().__dict__.copy(),
        }
