from dataclasses import dataclass
import math
import numpy as np
from scipy.spatial import cKDTree

from .gaussian_general import validate_spd
from .edge_importance_v17 import (
    exact_edge_importance,
    pair_specific_overlap_upper_bound,
    safe_global_overlap_radius,
    _diagonal_h_norms,
)
from .dynamic_graph_aims import _kinematics


@dataclass(frozen=True)
class SampledSparseAuditV18:
    step: int
    sampled_pairs: tuple
    priority_pairs: tuple
    random_pairs: tuple
    omitted_pair_count_estimate: int
    maximum_score: float
    rms_score: float
    maximum_overlap: float
    maximum_hamiltonian_relative: float
    maximum_time_connection_dt: float
    violation_count: int
    passed: bool

    def as_dict(self):
        return {
            "step":int(self.step),
            "sampled_pairs":[list(x) for x in self.sampled_pairs],
            "priority_pairs":[list(x) for x in self.priority_pairs],
            "random_pairs":[list(x) for x in self.random_pairs],
            "omitted_pair_count_estimate":
                int(self.omitted_pair_count_estimate),
            "sample_count":int(len(self.sampled_pairs)),
            "maximum_score":float(self.maximum_score),
            "rms_score":float(self.rms_score),
            "maximum_overlap":float(self.maximum_overlap),
            "maximum_hamiltonian_relative":
                float(self.maximum_hamiltonian_relative),
            "maximum_time_connection_dt":
                float(self.maximum_time_connection_dt),
            "violation_count":int(self.violation_count),
            "passed":bool(self.passed),
        }


def _kinematic_arrays(basis,provider):
    qdots=[]; pdots=[]
    for b in basis:
        qdot,pdot=_kinematics(b,provider)
        qdots.append(qdot)
        pdots.append(pdot)
    return np.asarray(qdots,float),np.asarray(pdots,float)


def _active_edge_set(update):
    return {
        (min(int(i),int(j)),max(int(i),int(j)))
        for i,j in update.active_edges
    }


def _priority_omitted_pairs(
    basis,
    update,
    *,
    search_overlap_floor,
    wider_search_factor,
    priority_count,
):
    """Near-cutoff omitted pairs found with a wider but still local geometric query."""
    n=len(basis)
    if n<2 or priority_count<=0:
        return []

    q=np.asarray([b.q for b in basis],dtype=float)
    amin=np.asarray([
        np.min(np.linalg.eigvalsh(validate_spd(b.A)))
        for b in basis
    ],dtype=float)

    wider_floor=max(
        float(search_overlap_floor)*float(wider_search_factor),
        1e-12,
    )
    radius=safe_global_overlap_radius(
        float(np.min(amin)),wider_floor
    )
    candidate=sorted(
        tuple(sorted(x))
        for x in cKDTree(q).query_pairs(
            radius,output_type="set"
        )
    )

    active=_active_edge_set(update)
    scored=[]
    for i,j in candidate:
        edge=(i,j)
        if edge in active:
            continue
        upper=pair_specific_overlap_upper_bound(
            q[i],q[j],amin[i],amin[j]
        )
        scored.append((float(upper),edge))

    scored.sort(key=lambda x:(-x[0],x[1]))
    return [edge for _,edge in scored[:int(priority_count)]]


def _uniform_random_omitted_pairs(
    n,
    active,
    excluded,
    count,
    seed,
):
    """Uniform pair-index rejection sampling without materializing all N(N-1)/2 pairs."""
    count=max(int(count),0)
    if count==0 or n<2:
        return []

    rng=np.random.default_rng(int(seed))
    out=set()
    excluded=set(excluded)
    active=set(active)

    total=n*(n-1)//2
    omitted=max(total-len(active),0)
    target=min(count,max(omitted-len(excluded),0))
    if target<=0:
        return []

    max_attempts=max(100,40*target)
    attempts=0
    while len(out)<target and attempts<max_attempts:
        i=int(rng.integers(0,n))
        j=int(rng.integers(0,n-1))
        if j>=i:
            j+=1
        edge=(min(i,j),max(i,j))
        attempts+=1
        if edge in active or edge in excluded:
            continue
        out.add(edge)

    # Deterministic exhaustive fallback only if rejection sampling struggled.
    if len(out)<target:
        for i in range(n):
            for j in range(i+1,n):
                edge=(i,j)
                if edge in active or edge in excluded or edge in out:
                    continue
                out.add(edge)
                if len(out)>=target:
                    break
            if len(out)>=target:
                break

    return sorted(out)


def sampled_omitted_edge_audit_v18(
    basis,
    provider,
    dt,
    update,
    graph_settings,
    *,
    step,
    priority_count=8,
    random_count=8,
    wider_search_factor=0.1,
    seed=20260813,
    violation_factor=1.0,
):
    r"""Stratified exact audit of omitted pair importance.

    Sampling has two components:

    1. priority sample: omitted pairs closest to the geometric search boundary;
    2. random sample: omitted pairs drawn approximately uniformly from all pair indices.

    A violation occurs if an omitted sampled edge has exact S/H/T score greater than

        violation_factor * current enter score.

    Such an edge should have been active and indicates that the sparse search/selection
    is too aggressive for the current state.
    """
    basis=list(basis)
    n=len(basis)
    if n==0:
        raise ValueError("basis cannot be empty.")

    priority=_priority_omitted_pairs(
        basis,update,
        search_overlap_floor=graph_settings.search_overlap_floor,
        wider_search_factor=wider_search_factor,
        priority_count=priority_count,
    )

    active=_active_edge_set(update)
    random_pairs=_uniform_random_omitted_pairs(
        n,active,priority,random_count,
        seed=int(seed)+int(step)*1009,
    )
    sample=tuple(sorted(set(priority)|set(random_pairs)))

    total=n*(n-1)//2
    omitted_count=max(total-len(active),0)

    if not sample:
        return SampledSparseAuditV18(
            step=int(step),
            sampled_pairs=(),
            priority_pairs=tuple(priority),
            random_pairs=tuple(random_pairs),
            omitted_pair_count_estimate=omitted_count,
            maximum_score=0.0,
            rms_score=0.0,
            maximum_overlap=0.0,
            maximum_hamiltonian_relative=0.0,
            maximum_time_connection_dt=0.0,
            violation_count=0,
            passed=True,
        )

    qdots,pdots=_kinematic_arrays(basis,provider)
    diag_h,Minv,params=_diagonal_h_norms(
        update.cache,provider
    )

    infos=[]
    for i,j in sample:
        infos.append(
            exact_edge_importance(
                update.cache,
                i,j,
                provider,
                qdots,pdots,
                float(dt),
                graph_settings,
                diagonal_h_norms=diag_h,
                Minv=Minv,
                params=params,
            )
        )

    scores=np.asarray([x.score for x in infos],float)
    threshold=float(
        violation_factor*graph_settings.enter_score
    )
    violations=int(np.sum(scores>threshold))

    return SampledSparseAuditV18(
        step=int(step),
        sampled_pairs=sample,
        priority_pairs=tuple(priority),
        random_pairs=tuple(random_pairs),
        omitted_pair_count_estimate=omitted_count,
        maximum_score=float(np.max(scores)),
        rms_score=float(np.sqrt(np.mean(scores**2))),
        maximum_overlap=float(max(x.overlap for x in infos)),
        maximum_hamiltonian_relative=float(
            max(x.hamiltonian_relative for x in infos)
        ),
        maximum_time_connection_dt=float(
            max(x.time_connection_dt for x in infos)
        ),
        violation_count=violations,
        passed=bool(violations==0),
    )
