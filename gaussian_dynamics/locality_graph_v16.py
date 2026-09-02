from dataclasses import dataclass
import math
import numpy as np
from scipy.spatial import cKDTree

from .gaussian_general import validate_spd
from .pair_cache_v15 import GaussianPairCache


def _position_overlap_bound_from_precomputed(qi,qj,ai,aj):
    qi=np.asarray(qi,dtype=float)
    qj=np.asarray(qj,dtype=float)
    if qi.shape!=qj.shape or qi.ndim!=1:
        raise ValueError("TBF coordinates must be equal-length vectors.")
    ai=float(ai); aj=float(aj)
    if ai<=0.0 or aj<=0.0:
        raise ValueError("minimum width eigenvalues must be positive.")
    h=1.0/(1.0/ai+1.0/aj)
    dq=qi-qj
    return float(math.exp(-0.5*h*float(dq@dq)))


def conservative_position_overlap_bound(tbf_i,tbf_j):
    r"""Cheap upper bound on the magnitude of the nuclear Gaussian overlap.

    For normalized real-width Gaussians,

        |<g_i|g_j>|

    is no larger than the zero-momentum spatial overlap.  Its displacement exponent
    contains

        H_ij = (A_i^{-1}+A_j^{-1})^{-1}.

    Since

        lambda_min(H_ij)
        >= 1 / (1/lambda_min(A_i) + 1/lambda_min(A_j)),

    and the unequal-width determinant prefactor is <= 1,

        |S_ij|
        <= exp[-0.5 h_ij ||q_i-q_j||^2]

    with

        h_ij =
        1 / (1/a_i + 1/a_j),
        a_i = lambda_min(A_i).

    Momentum mismatch can only decrease the exact overlap magnitude, so it is safe to
    omit it from this upper bound.

    The bound costs only O(d) once each TBF's smallest width eigenvalue has been
    precomputed.
    """
    qi=np.asarray(tbf_i.q,dtype=float)
    qj=np.asarray(tbf_j.q,dtype=float)
    if qi.shape!=qj.shape or qi.ndim!=1:
        raise ValueError("TBF coordinates must be equal-length vectors.")

    Ai=validate_spd(tbf_i.A)
    Aj=validate_spd(tbf_j.A)

    ai=float(np.min(np.linalg.eigvalsh(Ai)))
    aj=float(np.min(np.linalg.eigvalsh(Aj)))
    if ai<=0.0 or aj<=0.0:
        raise ValueError("Gaussian widths must be positive definite.")

    return _position_overlap_bound_from_precomputed(
        qi,qj,ai,aj
    )


@dataclass(frozen=True)
class LocalityGraphSettings:
    """Hysteretic overlap-graph policy."""

    enter_overlap: float = 0.05
    exit_overlap: float = 0.025
    exact_drop: bool = True
    use_kdtree: bool = True

    def validate(self):
        if not (0.0<self.exit_overlap<=self.enter_overlap<=1.0):
            raise ValueError(
                "Require 0 < exit_overlap <= enter_overlap <= 1."
            )
        return self


@dataclass
class LocalityGraphUpdate:
    active_edges: tuple
    cache: GaussianPairCache
    spatial_candidate_pairs: int
    globally_screened_pairs: int
    candidate_pairs: int
    exact_pair_checks: int
    screened_pairs: int
    entered_edges: int
    exited_edges: int
    retained_edges: int
    total_offdiagonal_pairs: int

    @property
    def active_offdiagonal_edges(self):
        return len(self.active_edges)

    @property
    def edge_fraction(self):
        denom=max(self.total_offdiagonal_pairs,1)
        return float(len(self.active_edges)/denom)

    @property
    def sparsity_fraction(self):
        return float(1.0-self.edge_fraction)

    def as_dict(self):
        return {
            "active_offdiagonal_edges":
                int(self.active_offdiagonal_edges),
            "spatial_candidate_pairs":
                int(self.spatial_candidate_pairs),
            "globally_screened_pairs":
                int(self.globally_screened_pairs),
            "candidate_pairs":int(self.candidate_pairs),
            "exact_pair_checks":int(self.exact_pair_checks),
            "screened_pairs":int(self.screened_pairs),
            "entered_edges":int(self.entered_edges),
            "exited_edges":int(self.exited_edges),
            "retained_edges":int(self.retained_edges),
            "total_offdiagonal_pairs":
                int(self.total_offdiagonal_pairs),
            "edge_fraction":self.edge_fraction,
            "sparsity_fraction":self.sparsity_fraction,
            "pair_cache_solves":
                int(self.cache.stats.canonical_solves),
        }


def _uid_pair(a,b):
    a=int(a); b=int(b)
    return (a,b) if a<b else (b,a)


class PersistentGaussianLocalityGraph:
    r"""Persistent overlap topology with enter/exit hysteresis.

    Topology is keyed by TBF uid rather than transient list index.

    The graph uses a two-stage screen:

    1. conservative position-only upper bound;
    2. exact overlap only when the bound says the pair could be relevant.

    Existing edges use the smaller `exit_overlap` threshold.  New edges require the
    larger `enter_overlap` threshold.  This suppresses edge flicker as moving Gaussians
    hover near the locality cutoff.
    """

    def __init__(self,settings=LocalityGraphSettings()):
        self.settings=settings.validate()
        self._active_uid_edges=set()
        self.update_count=0
        self.total_exact_pair_checks=0
        self.total_screened_pairs=0
        self.total_entered_edges=0
        self.total_exited_edges=0

    @property
    def active_uid_edges(self):
        return tuple(sorted(self._active_uid_edges))

    def clear(self):
        self._active_uid_edges.clear()

    def update(self,basis,cache=None):
        basis=list(basis)
        if not basis:
            raise ValueError("basis cannot be empty.")

        uids=[int(b.uid) for b in basis]
        if len(set(uids))!=len(uids):
            raise ValueError("TBF uids must be unique.")

        uid_to_index={uid:i for i,uid in enumerate(uids)}
        live=set(uids)

        q_arrays=[
            np.asarray(b.q,dtype=float)
            for b in basis
        ]
        min_width_eigs=[
            float(np.min(np.linalg.eigvalsh(validate_spd(b.A))))
            for b in basis
        ]

        old_edges={
            edge for edge in self._active_uid_edges
            if edge[0] in live and edge[1] in live
        }

        if cache is None:
            cache=GaussianPairCache(basis)
        else:
            if len(cache)!=len(basis):
                raise ValueError("provided cache has incompatible basis size.")
            cache_uids=[int(b.uid) for b in cache.basis]
            if cache_uids!=uids:
                raise ValueError(
                    "provided cache basis uid order must match locality-update basis."
                )

        active=set()
        candidate_pairs=0
        exact_checks=0
        screened=0
        entered=0
        retained=0

        n=len(basis)
        total=n*(n-1)//2

        if self.settings.use_kdtree and n>1:
            # Safe global radius from the smallest width eigenvalue and the smaller
            # exit threshold.  Because h_ij >= a_min/2,
            #
            #   |S_ij| <= exp[-a_min ||dq||^2 / 4].
            #
            # Any pair outside this radius cannot satisfy even the exit threshold.
            amin=min(min_width_eigs)
            global_radius=math.sqrt(
                max(
                    -4.0*math.log(self.settings.exit_overlap)/amin,
                    0.0,
                )
            )
            spatial_pairs=sorted(
                tuple(sorted(x))
                for x in cKDTree(
                    np.asarray(q_arrays,float)
                ).query_pairs(
                    global_radius,
                    output_type="set",
                )
            )
        else:
            spatial_pairs=[
                (i,j)
                for i in range(n)
                for j in range(i+1,n)
            ]

        spatial_candidate_pairs=len(spatial_pairs)
        globally_screened=total-spatial_candidate_pairs

        for i,j in spatial_pairs:
            edge=_uid_pair(uids[i],uids[j])
            was_active=edge in old_edges
            threshold=(
                self.settings.exit_overlap
                if was_active
                else self.settings.enter_overlap
            )

            upper=_position_overlap_bound_from_precomputed(
                q_arrays[i],q_arrays[j],
                min_width_eigs[i],min_width_eigs[j],
            )
            if upper<threshold:
                screened+=1
                continue

            candidate_pairs+=1
            pair=cache.pair(i,j)
            exact_checks+=1
            magnitude=abs(pair.overlap)

            if self.settings.exact_drop:
                keep=magnitude>=threshold
            else:
                keep=True

            if keep:
                active.add(edge)
                if was_active:
                    retained+=1
                else:
                    entered+=1

        screened+=globally_screened
        exited=len(old_edges-active)
        self._active_uid_edges=active
        self.update_count+=1
        self.total_exact_pair_checks+=exact_checks
        self.total_screened_pairs+=screened
        self.total_entered_edges+=entered
        self.total_exited_edges+=exited

        active_indices=tuple(sorted(
            (
                min(uid_to_index[a],uid_to_index[b]),
                max(uid_to_index[a],uid_to_index[b]),
            )
            for a,b in active
        ))

        return LocalityGraphUpdate(
            active_edges=active_indices,
            cache=cache,
            spatial_candidate_pairs=
                spatial_candidate_pairs,
            globally_screened_pairs=
                globally_screened,
            candidate_pairs=candidate_pairs,
            exact_pair_checks=exact_checks,
            screened_pairs=screened,
            entered_edges=entered,
            exited_edges=exited,
            retained_edges=retained,
            total_offdiagonal_pairs=total,
        )

    def diagnostics(self):
        return {
            "updates":int(self.update_count),
            "active_edges":int(len(self._active_uid_edges)),
            "total_exact_pair_checks":
                int(self.total_exact_pair_checks),
            "total_screened_pairs":
                int(self.total_screened_pairs),
            "total_entered_edges":
                int(self.total_entered_edges),
            "total_exited_edges":
                int(self.total_exited_edges),
            "settings":{
                "enter_overlap":
                    float(self.settings.enter_overlap),
                "exit_overlap":
                    float(self.settings.exit_overlap),
                "exact_drop":
                    bool(self.settings.exact_drop),
                "use_kdtree":
                    bool(self.settings.use_kdtree),
            },
        }
