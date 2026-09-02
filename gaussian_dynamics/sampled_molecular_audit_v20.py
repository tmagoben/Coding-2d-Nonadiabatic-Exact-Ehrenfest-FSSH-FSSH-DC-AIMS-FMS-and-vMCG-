from dataclasses import dataclass
import math
import numpy as np
from scipy.spatial import cKDTree

from .gaussian_general import validate_spd
from .sparse_molecular_matrices_v20 import (
    _center_kinematics,
    _position_bound,
    _safe_radius,
    molecular_pair_data_v20,
)


@dataclass(frozen=True)
class SampledMolecularAuditV20:
    step: int
    sampled_pairs: tuple
    priority_pairs: tuple
    random_pairs: tuple
    maximum_score: float
    rms_score: float
    violation_count: int
    passed: bool
    new_centroid_pairs_scored: int

    def as_dict(self):
        return {
            "step":int(self.step),
            "sampled_pairs":[list(x) for x in self.sampled_pairs],
            "priority_pairs":[list(x) for x in self.priority_pairs],
            "random_pairs":[list(x) for x in self.random_pairs],
            "maximum_score":float(self.maximum_score),
            "rms_score":float(self.rms_score),
            "violation_count":int(self.violation_count),
            "passed":bool(self.passed),
            "new_centroid_pairs_scored":
                int(self.new_centroid_pairs_scored),
        }


def _uniform_random_omitted(
    n,active,excluded,count,seed
):
    active=set(active)
    excluded=set(excluded)
    target=min(
        max(int(count),0),
        max(n*(n-1)//2-len(active)-len(excluded),0),
    )
    if target<=0:
        return []

    rng=np.random.default_rng(int(seed))
    out=set()
    attempts=0
    max_attempts=max(100,50*target)
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

    if len(out)<target:
        for i in range(n):
            for j in range(i+1,n):
                edge=(i,j)
                if (
                    edge in active
                    or edge in excluded
                    or edge in out
                ):
                    continue
                out.add(edge)
                if len(out)>=target:
                    break
            if len(out)>=target:
                break
    return sorted(out)


def sampled_molecular_edge_audit_v20(
    basis,
    provider,
    dt,
    sparse_update,
    settings,
    *,
    step,
    priority_count=8,
    random_count=8,
    audit_search_factor=0.1,
    seed=20260813,
    violation_factor=1.0,
):
    """Sample exact omitted molecular S/H/T edge scores.

    The priority sample expands the geometric search floor by `audit_search_factor`,
    while the random sample protects against a purely geometric blind spot.
    """
    basis=list(basis)
    n=len(basis)
    active=set(
        tuple(x)
        for x in sparse_update.active_edges
    )

    q=np.asarray([b.q for b in basis],float)
    amin=np.asarray([
        float(np.min(
            np.linalg.eigvalsh(
                validate_spd(b.A)
            )
        ))
        for b in basis
    ])

    wider_floor=max(
        float(settings.search_overlap_floor)
        *float(audit_search_factor),
        1e-14,
    )
    if n>1:
        radius=_safe_radius(
            float(np.min(amin)),
            wider_floor,
        )
        near=set(
            tuple(sorted(x))
            for x in cKDTree(q).query_pairs(
                radius,output_type="set"
            )
        )
    else:
        near=set()

    scored=[]
    for i,j in near:
        if (i,j) in active:
            continue
        upper=_position_bound(
            q[i],q[j],amin[i],amin[j]
        )
        scored.append((upper,(i,j)))
    scored.sort(key=lambda x:(-x[0],x[1]))
    priority=[
        edge for _,edge in scored[:max(int(priority_count),0)]
    ]

    random_pairs=_uniform_random_omitted(
        n,active,priority,random_count,
        seed=int(seed)+1009*int(step),
    )
    sample=tuple(sorted(
        set(priority)|set(random_pairs)
    ))
    if not sample:
        return SampledMolecularAuditV20(
            step=int(step),
            sampled_pairs=(),
            priority_pairs=tuple(priority),
            random_pairs=tuple(random_pairs),
            maximum_score=0.0,
            rms_score=0.0,
            violation_count=0,
            passed=True,
            new_centroid_pairs_scored=0,
        )

    qdots,pdots,centers=_center_kinematics(
        basis,provider
    )
    diag_h_abs=np.abs(
        sparse_update.diagonal_H
    )

    scores=[]
    new_count=0
    for edge in sample:
        if edge in sparse_update.pair_data:
            pair=sparse_update.pair_data[edge]
        else:
            pair=molecular_pair_data_v20(
                basis,edge[0],edge[1],
                provider,float(dt),
                qdots,pdots,centers,
                diag_h_abs,settings,
            )
            new_count+=1
        scores.append(float(pair.score))

    arr=np.asarray(scores,float)
    threshold=float(
        violation_factor*settings.enter_score
    )
    violations=int(np.sum(arr>threshold))
    return SampledMolecularAuditV20(
        step=int(step),
        sampled_pairs=sample,
        priority_pairs=tuple(priority),
        random_pairs=tuple(random_pairs),
        maximum_score=float(np.max(arr)),
        rms_score=float(
            np.sqrt(np.mean(arr**2))
        ),
        violation_count=violations,
        passed=bool(violations==0),
        new_centroid_pairs_scored=new_count,
    )
