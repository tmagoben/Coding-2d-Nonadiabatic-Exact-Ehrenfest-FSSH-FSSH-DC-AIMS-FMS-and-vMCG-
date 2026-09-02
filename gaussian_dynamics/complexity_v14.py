from dataclasses import dataclass, field
from contextlib import contextmanager
import time


@dataclass(frozen=True)
class AsymptoticComplexity:
    """Human-readable leading-order complexity for the v0.14 algorithms."""

    matrix_build: str
    coefficient_solve: str
    defect_evaluation: str
    prepared_candidate_ranking: str
    pruning: str
    memory: str


def asymptotic_complexity():
    r"""Return the symbolic leading-order costs.

    Symbols
    -------
    N : number of nuclear Gaussian basis functions
    s : number of electronic states per Gaussian
    d : nuclear coordinate dimension
    G : number of diagnostic grid points
    K : number of residual candidate Gaussians

    Notes
    -----
    These expressions describe the dense reference implementation in this repository.
    They are not lower bounds on all possible implementations.
    """
    return AsymptoticComplexity(
        matrix_build=(
            "O(N^2 d^3 + N^2 s^2): unequal-width Gaussian pair algebra "
            "requires dense dxd solves/inverses; electronic blocks are size sxs"
        ),
        coefficient_solve=(
            "O((sN)^3) time and O((sN)^2) memory for the dense Cayley/Galerkin "
            "linear solve"
        ),
        defect_evaluation=(
            "O(N G s + s G log G + N^2 d^3 + (sN)^3): reconstruct Psi/Psidot, "
            "apply FFT kinetic energy, and solve the projected coefficient equation"
        ),
        prepared_candidate_ranking=(
            "O(K G (N+s) + N^2 K + N^3): candidate-grid contractions, "
            "orthogonalization against the current N-Gaussian span, and one dense "
            "nuclear-overlap factorization"
        ),
        pruning=(
            "O(N^3 + Ns) per pruning audit using one inverse/solve of the nuclear "
            "overlap matrix and exact leave-one-out projection-loss scores"
        ),
        memory=(
            "O((sN)^2 + K G + N G s): dense electronic matrices, prepared "
            "candidate grid values, and diagnostic wavefunction/basis arrays"
        ),
    )


@dataclass
class ComplexityLedger:
    """Runtime/counter ledger for one adaptive propagation."""

    matrix_build_calls: int = 0
    pair_matrix_evaluations: int = 0
    ordered_pair_equivalent: int = 0
    time_matrix_calls: int = 0
    cayley_solve_calls: int = 0
    defect_evaluations: int = 0
    candidate_ranking_calls: int = 0
    candidate_count_scored: int = 0
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
    pruning_seconds: float = 0.0
    total_seconds: float = 0.0

    _start: float = field(default=0.0, repr=False)

    def start(self):
        self._start=time.perf_counter()

    def stop(self):
        if self._start:
            self.total_seconds += time.perf_counter()-self._start
            self._start=0.0

    def observe_basis(self, n_basis, nstate=2):
        n=int(n_basis)
        self.peak_basis_size=max(self.peak_basis_size,n)
        self.peak_electronic_dimension=max(
            self.peak_electronic_dimension,
            int(nstate)*n,
        )

    @contextmanager
    def timed(self, category):
        t0=time.perf_counter()
        yield
        elapsed=time.perf_counter()-t0
        attr=f"{category}_seconds"
        if not hasattr(self,attr):
            raise ValueError(f"unknown timing category {category!r}")
        setattr(self,attr,getattr(self,attr)+elapsed)

    def as_dict(self):
        return {
            "matrix_build_calls":int(self.matrix_build_calls),
            "pair_matrix_evaluations":int(self.pair_matrix_evaluations),
            "ordered_pair_equivalent":int(self.ordered_pair_equivalent),
            "time_matrix_calls":int(self.time_matrix_calls),
            "cayley_solve_calls":int(self.cayley_solve_calls),
            "defect_evaluations":int(self.defect_evaluations),
            "candidate_ranking_calls":int(self.candidate_ranking_calls),
            "candidate_count_scored":int(self.candidate_count_scored),
            "enrichment_events":int(self.enrichment_events),
            "pruning_audits":int(self.pruning_audits),
            "pruning_events":int(self.pruning_events),
            "peak_basis_size":int(self.peak_basis_size),
            "peak_electronic_dimension":int(self.peak_electronic_dimension),
            "peak_candidate_count":int(self.peak_candidate_count),
            "matrix_build_seconds":float(self.matrix_build_seconds),
            "time_matrix_seconds":float(self.time_matrix_seconds),
            "cayley_solve_seconds":float(self.cayley_solve_seconds),
            "defect_seconds":float(self.defect_seconds),
            "candidate_ranking_seconds":float(self.candidate_ranking_seconds),
            "pruning_seconds":float(self.pruning_seconds),
            "total_seconds":float(self.total_seconds),
            "asymptotic":asymptotic_complexity().__dict__.copy(),
        }


def dense_dimension_cost_proxy(n_basis, nstate=2):
    """Simple cubic dense-solve proxy (sN)^3 used only for relative comparisons."""
    m=int(n_basis)*int(nstate)
    return int(m**3)


def pair_matrix_cost_proxy(n_basis, dimension):
    """Simple N^2 d^3 proxy for unequal-width Gaussian pair algebra."""
    return int(int(n_basis)**2 * int(dimension)**3)


def candidate_ranking_cost_proxy(n_basis, n_candidates, grid_points, nstate=2):
    """Proxy for K G (N+s) + N^2 K + N^3."""
    N=int(n_basis)
    K=int(n_candidates)
    G=int(grid_points)
    s=int(nstate)
    return int(K*G*(N+s) + N*N*K + N**3)
