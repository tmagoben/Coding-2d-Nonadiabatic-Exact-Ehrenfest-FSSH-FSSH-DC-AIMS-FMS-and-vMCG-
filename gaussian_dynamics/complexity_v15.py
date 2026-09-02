from dataclasses import dataclass, field
from contextlib import contextmanager
import time


@dataclass(frozen=True)
class AsymptoticComplexityV15:
    pair_snapshot: str
    cached_matrix_build: str
    cached_time_matrix: str
    incremental_add: str
    incremental_prune: str
    coefficient_solve: str
    defect_evaluation: str
    candidate_ranking: str
    pruning_audit: str
    memory: str


def asymptotic_complexity_v15():
    return AsymptoticComplexityV15(
        pair_snapshot=(
            "O(N^2 d^3) time and O(N^2 d^2) cached pair-moment memory; "
            "one B=A_i+A_j multi-RHS solve per canonical i<=j pair"
        ),
        cached_matrix_build=(
            "O(N^2 d^3 + N^2 s^2) overall; overlap/centroid/covariance are "
            "computed once per canonical pair and reused by S, kinetic, and V"
        ),
        cached_time_matrix=(
            "O(N^2 s^2) after a pair snapshot is available; T is ordered/non-Hermitian "
            "but reversed pair moments are obtained by conjugation rather than new dxd solves"
        ),
        incremental_add=(
            "O(N d^3 + N s^2) pair/matrix work at one fixed snapshot when one Gaussian "
            "is appended; old-old blocks are reused"
        ),
        incremental_prune=(
            "O((sN)^2) copying/slicing after the O(N^3) pruning audit; no Gaussian "
            "pair integrals are recomputed"
        ),
        coefficient_solve=(
            "O((sN)^3) time and O((sN)^2) memory for dense Cayley/Galerkin solves"
        ),
        defect_evaluation=(
            "O(N G s + s G log G + N^2 s^2 + (sN)^3) once endpoint pair moments "
            "are already cached"
        ),
        candidate_ranking=(
            "O(K G (N+s) + N^2 K + N^3) for vectorized defect contractions and "
            "orthogonalization; cost-aware reranking of a short list is O(K_short)"
        ),
        pruning_audit=(
            "O(N^3 + Ns) using one inverse/solve of the nuclear Gram matrix and exact "
            "leave-one-out projection losses"
        ),
        memory=(
            "O((sN)^2 + N^2 d^2 + K G + N G s): dense quantum matrices, pair cache, "
            "candidate grids, and diagnostic wavefunction arrays"
        ),
    )


@dataclass
class ComplexityLedgerV15:
    """Runtime, cache, and work ledger for the v0.15 adaptive runner."""

    full_matrix_builds: int = 0
    incremental_expansions: int = 0
    incremental_prunes: int = 0

    pair_snapshots: int = 0
    pair_requests: int = 0
    pair_factorizations: int = 0
    propagation_pair_factorizations: int = 0
    candidate_pair_factorizations: int = 0
    pair_direct_hits: int = 0
    pair_reverse_views: int = 0
    inherited_pairs_reused: int = 0

    v14_factorization_baseline: int = 0
    factorization_avoided: int = 0

    time_matrix_calls: int = 0
    cayley_solve_calls: int = 0
    cayley_cubic_units: int = 0
    defect_solve_calls: int = 0
    defect_cubic_units: int = 0

    defect_evaluations: int = 0
    candidate_ranking_calls: int = 0
    candidate_count_scored: int = 0
    cost_ranking_calls: int = 0
    enrichment_events: int = 0
    pruning_audits: int = 0
    pruning_events: int = 0

    peak_basis_size: int = 0
    peak_electronic_dimension: int = 0
    peak_candidate_count: int = 0

    matrix_build_seconds: float = 0.0
    time_matrix_seconds: float = 0.0
    cayley_solve_seconds: float = 0.0
    defect_seconds: float = 0.0
    candidate_ranking_seconds: float = 0.0
    cost_ranking_seconds: float = 0.0
    pruning_seconds: float = 0.0
    total_seconds: float = 0.0

    _start: float = field(default=0.0,repr=False)

    def start(self):
        self._start=time.perf_counter()

    def stop(self):
        if self._start:
            self.total_seconds+=time.perf_counter()-self._start
            self._start=0.0

    def observe_basis(self,n_basis,nstate=2):
        N=int(n_basis)
        self.peak_basis_size=max(self.peak_basis_size,N)
        self.peak_electronic_dimension=max(
            self.peak_electronic_dimension,
            int(nstate)*N,
        )

    def add_cache_stats(self,stats):
        self.pair_requests+=int(stats.requests)
        self.pair_factorizations+=int(stats.canonical_solves)
        self.pair_direct_hits+=int(stats.direct_hits)
        self.pair_reverse_views+=int(stats.reverse_views)
        self.inherited_pairs_reused+=int(stats.inherited_pairs)

    def add_cache_delta(self,stats,before,category="propagation"):
        """Accumulate only work performed since a previous cache-stat snapshot."""
        self.pair_requests+=int(stats.requests-before.get("requests",0))
        delta_factorizations=int(
            stats.canonical_solves-before.get("canonical_solves",0)
        )
        self.pair_factorizations+=delta_factorizations
        if category=="propagation":
            self.propagation_pair_factorizations+=delta_factorizations
        elif category=="candidate":
            self.candidate_pair_factorizations+=delta_factorizations
        else:
            raise ValueError("category must be propagation or candidate.")
        self.pair_direct_hits+=int(
            stats.direct_hits-before.get("direct_hits",0)
        )
        self.pair_reverse_views+=int(
            stats.reverse_views-before.get("reverse_views",0)
        )
        self.inherited_pairs_reused+=int(
            stats.inherited_pairs-before.get("inherited_pairs",0)
        )

    @contextmanager
    def timed(self,category):
        t0=time.perf_counter()
        yield
        elapsed=time.perf_counter()-t0
        attr=f"{category}_seconds"
        if not hasattr(self,attr):
            raise ValueError(f"unknown timing category {category!r}")
        setattr(self,attr,getattr(self,attr)+elapsed)

    @property
    def cache_hit_fraction(self):
        denom=max(self.pair_requests,1)
        return float(
            (self.pair_direct_hits+self.pair_reverse_views)
            /denom
        )

    def finalize_avoided(self):
        self.factorization_avoided=max(
            int(self.v14_factorization_baseline)
            -int(self.propagation_pair_factorizations),
            0,
        )

    def empirical_cost_rates(self):
        """Observed seconds per pair factorization and per cubic solve unit."""
        pair_seconds=(
            self.matrix_build_seconds
            +self.time_matrix_seconds
        )
        pair_rate=(
            pair_seconds/max(self.pair_factorizations,1)
        )

        solve_seconds=(
            self.cayley_solve_seconds
        )
        solve_rate=(
            solve_seconds/max(self.cayley_cubic_units,1)
        )

        return {
            "pair_seconds_per_factorization":
                float(pair_rate),
            "cayley_seconds_per_cubic_unit":
                float(solve_rate),
        }

    def as_dict(self):
        self.finalize_avoided()
        return {
            "full_matrix_builds":int(self.full_matrix_builds),
            "incremental_expansions":int(self.incremental_expansions),
            "incremental_prunes":int(self.incremental_prunes),
            "pair_snapshots":int(self.pair_snapshots),
            "pair_requests":int(self.pair_requests),
            "pair_factorizations":int(self.pair_factorizations),
            "propagation_pair_factorizations":
                int(self.propagation_pair_factorizations),
            "candidate_pair_factorizations":
                int(self.candidate_pair_factorizations),
            "pair_direct_hits":int(self.pair_direct_hits),
            "pair_reverse_views":int(self.pair_reverse_views),
            "inherited_pairs_reused":int(self.inherited_pairs_reused),
            "cache_hit_fraction":self.cache_hit_fraction,
            "v14_factorization_baseline":
                int(self.v14_factorization_baseline),
            "factorization_avoided":
                int(self.factorization_avoided),
            "factorization_reduction_fraction":
                float(
                    self.factorization_avoided
                    /max(self.v14_factorization_baseline,1)
                ),
            "time_matrix_calls":int(self.time_matrix_calls),
            "cayley_solve_calls":int(self.cayley_solve_calls),
            "cayley_cubic_units":int(self.cayley_cubic_units),
            "defect_solve_calls":int(self.defect_solve_calls),
            "defect_cubic_units":int(self.defect_cubic_units),
            "defect_evaluations":int(self.defect_evaluations),
            "candidate_ranking_calls":
                int(self.candidate_ranking_calls),
            "candidate_count_scored":
                int(self.candidate_count_scored),
            "cost_ranking_calls":int(self.cost_ranking_calls),
            "enrichment_events":int(self.enrichment_events),
            "pruning_audits":int(self.pruning_audits),
            "pruning_events":int(self.pruning_events),
            "peak_basis_size":int(self.peak_basis_size),
            "peak_electronic_dimension":
                int(self.peak_electronic_dimension),
            "peak_candidate_count":int(self.peak_candidate_count),
            "matrix_build_seconds":float(self.matrix_build_seconds),
            "time_matrix_seconds":float(self.time_matrix_seconds),
            "cayley_solve_seconds":float(self.cayley_solve_seconds),
            "defect_seconds":float(self.defect_seconds),
            "candidate_ranking_seconds":
                float(self.candidate_ranking_seconds),
            "cost_ranking_seconds":
                float(self.cost_ranking_seconds),
            "pruning_seconds":float(self.pruning_seconds),
            "total_seconds":float(self.total_seconds),
            "empirical_cost_rates":
                self.empirical_cost_rates(),
            "asymptotic":
                asymptotic_complexity_v15().__dict__.copy(),
        }


def dense_cubic_units(n_basis,nstate=2):
    return int((int(n_basis)*int(nstate))**3)


def canonical_pair_count(n_basis):
    N=int(n_basis)
    return N*(N+1)//2


def incremental_pair_count_for_add(n_basis):
    """New canonical pairs when appending one Gaussian at a fixed snapshot."""
    return int(n_basis)+1
