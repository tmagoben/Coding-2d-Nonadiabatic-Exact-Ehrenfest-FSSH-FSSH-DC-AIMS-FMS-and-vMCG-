from dataclasses import dataclass
import math
import numpy as np
from scipy.spatial import cKDTree

from .gaussian_general import validate_spd
from .pair_cache_v15 import GaussianPairCache
from .dynamic_graph_aims import _kinematics


@dataclass(frozen=True)
class EdgeImportanceSettingsV17:
    """Dimensionless S/H/T edge-importance policy.

    `search_overlap_floor` controls only the safe geometric pre-screen.  Final graph
    membership is controlled by the combined S/H/T score.
    """

    enter_score: float = 0.030
    exit_score: float = 0.015
    search_overlap_floor: float = 1e-4

    overlap_weight: float = 1.0
    hamiltonian_weight: float = 0.20
    time_connection_weight: float = 1.0
    energy_floor: float = 1e-6
    local_omitted_score_l2_budget: float = 0.08

    use_kdtree: bool = True

    def validate(self):
        if not (0.0 < self.exit_score <= self.enter_score):
            raise ValueError("Require 0 < exit_score <= enter_score.")
        if not (0.0 < self.search_overlap_floor <= 1.0):
            raise ValueError("search_overlap_floor must be in (0,1].")
        for name,value in (
            ("overlap_weight",self.overlap_weight),
            ("hamiltonian_weight",self.hamiltonian_weight),
            ("time_connection_weight",self.time_connection_weight),
        ):
            if value < 0.0:
                raise ValueError(f"{name} cannot be negative.")
        if self.energy_floor <= 0.0:
            raise ValueError("energy_floor must be positive.")
        if self.local_omitted_score_l2_budget < 0.0:
            raise ValueError("local_omitted_score_l2_budget cannot be negative.")
        return self


@dataclass(frozen=True)
class EdgeImportance:
    i: int
    j: int
    overlap: float
    hamiltonian_relative: float
    time_connection_dt: float
    score: float
    hamiltonian_block_norm: float
    time_pair_norm: float

    def as_dict(self):
        return {
            "i":int(self.i),
            "j":int(self.j),
            "overlap":float(self.overlap),
            "hamiltonian_relative":
                float(self.hamiltonian_relative),
            "time_connection_dt":
                float(self.time_connection_dt),
            "score":float(self.score),
            "hamiltonian_block_norm":
                float(self.hamiltonian_block_norm),
            "time_pair_norm":
                float(self.time_pair_norm),
        }


@dataclass
class EdgeControlledGraphUpdateV17:
    active_edges: tuple
    cache: GaussianPairCache
    importance: dict

    spatial_candidate_pairs: int
    globally_screened_pairs: int
    pair_bound_screened_pairs: int
    exact_pair_checks: int
    screened_pairs: int
    entered_edges: int
    exited_edges: int
    retained_edges: int
    total_offdiagonal_pairs: int

    omitted_candidate_score_l2: float
    omitted_candidate_score_max: float
    omitted_candidate_count: int
    budget_promoted_edges: int

    @property
    def candidate_pairs(self):
        return self.spatial_candidate_pairs

    @property
    def active_offdiagonal_edges(self):
        return len(self.active_edges)

    @property
    def edge_fraction(self):
        return float(
            self.active_offdiagonal_edges
            /max(self.total_offdiagonal_pairs,1)
        )

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
            "pair_bound_screened_pairs":
                int(self.pair_bound_screened_pairs),
            "exact_pair_checks":
                int(self.exact_pair_checks),
            "screened_pairs":
                int(self.screened_pairs),
            "entered_edges":int(self.entered_edges),
            "exited_edges":int(self.exited_edges),
            "retained_edges":int(self.retained_edges),
            "total_offdiagonal_pairs":
                int(self.total_offdiagonal_pairs),
            "edge_fraction":self.edge_fraction,
            "sparsity_fraction":self.sparsity_fraction,
            "omitted_candidate_score_l2":
                float(self.omitted_candidate_score_l2),
            "omitted_candidate_score_max":
                float(self.omitted_candidate_score_max),
            "omitted_candidate_count":
                int(self.omitted_candidate_count),
            "budget_promoted_edges":
                int(self.budget_promoted_edges),
            "pair_cache_solves":
                int(self.cache.stats.canonical_solves),
        }


def _uid_pair(a,b):
    a=int(a); b=int(b)
    return (a,b) if a<b else (b,a)


def _minimum_width_eigenvalues(basis):
    return np.asarray([
        np.min(np.linalg.eigvalsh(validate_spd(b.A)))
        for b in basis
    ],dtype=float)


def safe_global_overlap_radius(min_width_eigenvalue,overlap_floor):
    """Radius outside which the conservative overlap bound is below `overlap_floor`."""
    a=float(min_width_eigenvalue)
    tau=float(overlap_floor)
    if a<=0.0:
        raise ValueError("minimum width eigenvalue must be positive.")
    if not (0.0<tau<=1.0):
        raise ValueError("overlap_floor must be in (0,1].")
    return float(math.sqrt(max(-4.0*math.log(tau)/a,0.0)))


def pair_specific_overlap_upper_bound(qi,qj,ai,aj):
    qi=np.asarray(qi,float)
    qj=np.asarray(qj,float)
    ai=float(ai); aj=float(aj)
    h=1.0/(1.0/ai+1.0/aj)
    dq=qi-qj
    return float(math.exp(-0.5*h*float(dq@dq)))


def _kinematics_arrays(basis,provider):
    qdots=[]; pdots=[]
    for b in basis:
        qdot,pdot=_kinematics(b,provider)
        qdots.append(qdot)
        pdots.append(pdot)
    return np.asarray(qdots,float),np.asarray(pdots,float)


def _diagonal_h_norms(cache,provider):
    point=provider.evaluate(
        np.asarray(cache.basis[0].q,float)
    )
    Minv=np.linalg.inv(
        np.asarray(point.mass_matrix,float)
    )
    params=provider.params
    eye=np.eye(2,dtype=complex)

    out=np.zeros(len(cache),dtype=float)
    for i in range(len(cache)):
        pair=cache.pair(i,i)
        block=(
            pair.kinetic(Minv)*eye
            +pair.lvc_potential_matrix(params)
        )
        out[i]=float(np.linalg.norm(block,ord="fro"))
    return out,Minv,params


def exact_edge_importance(
    cache,
    i,
    j,
    provider,
    qdots,
    pdots,
    dt,
    settings,
    diagonal_h_norms=None,
    Minv=None,
    params=None,
):
    """Evaluate a dimensionless S/H/T edge score from exact pair data."""
    settings=settings.validate()
    i=int(i); j=int(j)

    if diagonal_h_norms is None or Minv is None or params is None:
        diagonal_h_norms,Minv,params=_diagonal_h_norms(
            cache,provider
        )

    pair=cache.pair(i,j)
    eye=np.eye(2,dtype=complex)
    Hij=(
        pair.kinetic(Minv)*eye
        +pair.lvc_potential_matrix(params)
    )
    hnorm=float(np.linalg.norm(Hij,ord="fro"))

    hscale=max(
        math.sqrt(
            max(diagonal_h_norms[i],0.0)
            *max(diagonal_h_norms[j],0.0)
        ),
        float(settings.energy_floor),
    )
    hrel=hnorm/hscale

    tij=pair.time_element(
        qdots[j],pdots[j],None
    )
    tji=cache.pair(j,i).time_element(
        qdots[i],pdots[i],None
    )
    tpair=math.sqrt(
        abs(tij)**2+abs(tji)**2
    )
    tdt=float(abs(dt)*tpair)

    sabs=float(abs(pair.overlap))
    score=math.sqrt(
        (settings.overlap_weight*sabs)**2
        +(settings.hamiltonian_weight*hrel)**2
        +(settings.time_connection_weight*tdt)**2
    )

    return EdgeImportance(
        i=i,j=j,
        overlap=sabs,
        hamiltonian_relative=float(hrel),
        time_connection_dt=tdt,
        score=float(score),
        hamiltonian_block_norm=hnorm,
        time_pair_norm=float(tpair),
    )


class ErrorControlledGaussianLocalityGraphV17:
    """Persistent graph selected by exact local S/H/T importance, not overlap alone.

    The global KD-tree pre-screen remains overlap-based and conservative.  Pairs inside
    that safe search radius receive exact S/H/T importance scores.  Periodic dense
    audits in the v0.17 runner are responsible for checking whether the geometric
    search floor itself is sufficiently conservative for H/T.
    """

    def __init__(
        self,
        provider,
        dt,
        settings=EdgeImportanceSettingsV17(),
    ):
        self.provider=provider
        self.dt=float(dt)
        if self.dt<=0.0:
            raise ValueError("dt must be positive.")
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

    def relax_scores(self,factor):
        """One-sided error-control relaxation: smaller score thresholds retain edges."""
        factor=float(factor)
        if not (0.0<factor<1.0):
            raise ValueError("relaxation factor must lie in (0,1).")
        s=self.settings
        self.settings=EdgeImportanceSettingsV17(
            enter_score=max(s.enter_score*factor,1e-8),
            exit_score=max(s.exit_score*factor,1e-8),
            search_overlap_floor=s.search_overlap_floor,
            overlap_weight=s.overlap_weight,
            hamiltonian_weight=s.hamiltonian_weight,
            time_connection_weight=s.time_connection_weight,
            energy_floor=s.energy_floor,
            local_omitted_score_l2_budget=
                s.local_omitted_score_l2_budget,
            use_kdtree=s.use_kdtree,
        ).validate()

    def relax_search_floor(self,factor):
        """Lower the geometric search floor, increasing the safe KD-tree radius."""
        factor=float(factor)
        if not (0.0<factor<1.0):
            raise ValueError("relaxation factor must lie in (0,1).")
        s=self.settings
        self.settings=EdgeImportanceSettingsV17(
            enter_score=s.enter_score,
            exit_score=s.exit_score,
            search_overlap_floor=max(
                s.search_overlap_floor*factor,
                1e-12,
            ),
            overlap_weight=s.overlap_weight,
            hamiltonian_weight=s.hamiltonian_weight,
            time_connection_weight=s.time_connection_weight,
            energy_floor=s.energy_floor,
            local_omitted_score_l2_budget=
                s.local_omitted_score_l2_budget,
            use_kdtree=s.use_kdtree,
        ).validate()

    def update(self,basis,cache=None):
        basis=list(basis)
        if not basis:
            raise ValueError("basis cannot be empty.")

        uids=[int(b.uid) for b in basis]
        if len(set(uids))!=len(uids):
            raise ValueError("TBF uids must be unique.")

        if cache is None:
            cache=GaussianPairCache(basis)
        else:
            if len(cache)!=len(basis):
                raise ValueError("provided cache has incompatible size.")
            if [int(b.uid) for b in cache.basis]!=uids:
                raise ValueError("cache basis uid order must match.")

        n=len(basis)
        total=n*(n-1)//2
        q=np.asarray([b.q for b in basis],dtype=float)
        amin=_minimum_width_eigenvalues(basis)
        qdots,pdots=_kinematics_arrays(
            basis,self.provider
        )

        old={
            edge for edge in self._active_uid_edges
            if edge[0] in set(uids)
            and edge[1] in set(uids)
        }
        uid_to_index={uid:i for i,uid in enumerate(uids)}

        if self.settings.use_kdtree and n>1:
            radius=safe_global_overlap_radius(
                float(np.min(amin)),
                self.settings.search_overlap_floor,
            )
            spatial_pairs=sorted(
                tuple(sorted(x))
                for x in cKDTree(q).query_pairs(
                    radius,output_type="set"
                )
            )
        else:
            spatial_pairs=[
                (i,j)
                for i in range(n)
                for j in range(i+1,n)
            ]

        globally_screened=total-len(spatial_pairs)

        diagonal_h,Minv,params=_diagonal_h_norms(
            cache,self.provider
        )

        active=set()
        importance={}
        entered=0
        retained=0
        omitted_items=[]

        pair_bound_screened=0
        exact_checks=0

        for i,j in spatial_pairs:
            upper=pair_specific_overlap_upper_bound(
                q[i],q[j],amin[i],amin[j]
            )
            if upper<self.settings.search_overlap_floor:
                pair_bound_screened+=1
                continue

            info=exact_edge_importance(
                cache,i,j,self.provider,
                qdots,pdots,self.dt,self.settings,
                diagonal_h_norms=diagonal_h,
                Minv=Minv,
                params=params,
            )
            exact_checks+=1
            importance[(i,j)]=info

            edge=_uid_pair(uids[i],uids[j])
            was_active=edge in old
            threshold=(
                self.settings.exit_score
                if was_active
                else self.settings.enter_score
            )
            if info.score>=threshold:
                active.add(edge)
                if was_active:
                    retained+=1
                else:
                    entered+=1
            else:
                omitted_items.append((edge,info))

        # Global local-importance budget proxy. This is not a rigorous matrix-norm
        # theorem; it prevents many individually-small local scores from accumulating
        # without bound before the periodic dense S/H audit.
        budget=float(
            self.settings.local_omitted_score_l2_budget
        )
        promoted=0
        if omitted_items and budget>=0.0:
            remaining2=sum(
                float(info.score)**2
                for _,info in omitted_items
            )
            if math.sqrt(max(remaining2,0.0))>budget:
                for edge,info in sorted(
                    omitted_items,
                    key=lambda item:item[1].score,
                    reverse=True,
                ):
                    if math.sqrt(max(remaining2,0.0))<=budget:
                        break
                    active.add(edge)
                    if edge in old:
                        retained+=1
                    else:
                        entered+=1
                    remaining2=max(
                        remaining2-float(info.score)**2,
                        0.0,
                    )
                    promoted+=1

        # Only genuinely omitted candidates contribute to the reported residual budget.
        omitted_items=[
            (edge,info)
            for edge,info in omitted_items
            if edge not in active
        ]

        exited=len(old-active)
        self._active_uid_edges=active
        self.update_count+=1
        self.total_exact_pair_checks+=exact_checks
        self.total_screened_pairs+=(
            globally_screened+pair_bound_screened
        )
        self.total_entered_edges+=entered
        self.total_exited_edges+=exited

        active_indices=tuple(sorted(
            (
                min(uid_to_index[a],uid_to_index[b]),
                max(uid_to_index[a],uid_to_index[b]),
            )
            for a,b in active
        ))

        omitted=np.asarray(
            [info.score for _,info in omitted_items],
            dtype=float,
        )
        return EdgeControlledGraphUpdateV17(
            active_edges=active_indices,
            cache=cache,
            importance=importance,
            spatial_candidate_pairs=len(spatial_pairs),
            globally_screened_pairs=globally_screened,
            pair_bound_screened_pairs=
                pair_bound_screened,
            exact_pair_checks=exact_checks,
            screened_pairs=(
                globally_screened+pair_bound_screened
            ),
            entered_edges=entered,
            exited_edges=exited,
            retained_edges=retained,
            total_offdiagonal_pairs=total,
            omitted_candidate_score_l2=float(
                np.linalg.norm(omitted)
                if omitted.size else 0.0
            ),
            omitted_candidate_score_max=float(
                np.max(omitted)
                if omitted.size else 0.0
            ),
            omitted_candidate_count=int(omitted.size),
            budget_promoted_edges=int(promoted),
        )

    def diagnostics(self):
        s=self.settings
        return {
            "updates":int(self.update_count),
            "active_edges":
                int(len(self._active_uid_edges)),
            "total_exact_pair_checks":
                int(self.total_exact_pair_checks),
            "total_screened_pairs":
                int(self.total_screened_pairs),
            "total_entered_edges":
                int(self.total_entered_edges),
            "total_exited_edges":
                int(self.total_exited_edges),
            "settings":{
                "enter_score":float(s.enter_score),
                "exit_score":float(s.exit_score),
                "search_overlap_floor":
                    float(s.search_overlap_floor),
                "overlap_weight":
                    float(s.overlap_weight),
                "hamiltonian_weight":
                    float(s.hamiltonian_weight),
                "time_connection_weight":
                    float(s.time_connection_weight),
                "local_omitted_score_l2_budget":
                    float(s.local_omitted_score_l2_budget),
            },
        }
