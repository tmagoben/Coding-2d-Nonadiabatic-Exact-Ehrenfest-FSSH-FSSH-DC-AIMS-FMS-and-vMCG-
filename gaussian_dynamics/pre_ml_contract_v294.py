"""Pre-ML oracle/data contract for Gaussian Nonadiabatic Dynamics v0.29.4.

This layer adds no new dynamics. It freezes a canonical, versioned interface between
exact residual-driven curved spawning (v0.29.3) and any future ML ranker. The exact
physics evaluator remains authoritative. ML-ready records use invariant or canonical
features, carry trajectory-level split identities, preserve rejection reasons, and
include provenance hashes so labels cannot silently drift across software revisions.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
import hashlib, json
import numpy as np

from .curved_basis_adaptation_v292 import CurvedBasisSettingsV292
from .curved_basis_adaptation_v293 import evaluate_curved_spawn_candidate_v293, candidate_coefficient_block_v293
from .curved_geometry_v293 import straight_path_v293

V294_ORACLE_SCHEMA = "gnd-spawn-oracle-v1"
V294_DATASET_SCHEMA = "gnd-spawn-dataset-v1"
V294_FEATURE_SCHEMA = "gnd-spawn-features-invariant-v1"
V294_SPLIT_SCHEMA = "trajectory-group-hash-v1"

class RejectionReasonV294(str, Enum):
    NO_NEW_DIRECTION = "no-new-curved-tangent-direction"
    ILL_CONDITIONED = "ill-conditioned-curved-candidate-space"
    LOW_GAIN = "insufficient-predicted-residual-reduction"
    LOW_RELATIVE_GAIN = "insufficient-relative-residual-reduction"
    PROJECTION_LOSS = "projection-loss"
    ENERGY_JUMP = "energy-jump"
    DUPLICATE = "duplicate-or-non-novel"
    DEGENERATE_WIDTH = "degenerate-width-eigenspace"
    PATH_FAILURE = "path-or-transport-failure"
    NUMERICAL_FAILURE = "numerical-failure"


def _canon(x):
    if isinstance(x, np.generic): return _canon(x.item())
    if isinstance(x, np.ndarray): return _canon(x.tolist())
    if isinstance(x, complex): return [float(x.real), float(x.imag)]
    if isinstance(x, dict): return {str(k): _canon(v) for k,v in sorted(x.items(), key=lambda z:str(z[0]))}
    if isinstance(x, (list,tuple)): return [_canon(v) for v in x]
    if isinstance(x, float):
        if not np.isfinite(x): raise ValueError("non-finite canonical value")
        return float(x)
    if x is None or isinstance(x,(str,int,bool)): return x
    raise TypeError(type(x).__name__)


def stable_hash_v294(payload):
    raw=json.dumps(_canon(payload),sort_keys=True,separators=(",",":"))
    return hashlib.sha256(raw.encode()).hexdigest()


def canonical_state_payload_v294(state, *, packet_ids=None, path_family="straight"):
    """Deterministic serialization contract; packet rows are ordered by stable ID."""
    state=state.validate(require_normalized=False)
    ids=tuple(packet_ids or (f"g{i:06d}" for i in range(state.ngaussian)))
    if len(ids)!=state.ngaussian or len(set(ids))!=len(ids): raise ValueError("packet IDs must be unique and complete")
    order=np.argsort(np.asarray(ids,dtype=str),kind="stable")
    packets=[]
    for i in order:
        packets.append({
            "id": ids[i], "q": state.q[i], "p": state.p[i],
            "width": state.width_matrices[i], "chirp": state.chirp_matrices[i],
            "coefficients": state.coefficients[i],
        })
    return _canon({"schema":"gnd-canonical-state-v1","time_au":state.time_au,"path_family":path_family,"packets":packets})


def trajectory_split_v294(trajectory_group_id, *, train=80, validation=10, test=10):
    """Deterministic group split. Every candidate from a trajectory receives one split."""
    if train+validation+test!=100 or min(train,validation,test)<=0: raise ValueError("split percentages must be positive and total 100")
    h=int(hashlib.sha256((V294_SPLIT_SCHEMA+"|"+str(trajectory_group_id)).encode()).hexdigest()[:16],16)%100
    return "train" if h<train else ("validation" if h<train+validation else "test")


def _safe_eigvals(a): return np.linalg.eigvalsh(0.5*(np.asarray(a,float)+np.asarray(a,float).T))


def invariant_candidate_features_v294(state,candidate,snapshot):
    """Rotation/permutation-safe scalar features. No raw electronic gauge coordinates."""
    state=state.validate(require_normalized=False); candidate=candidate.validate(state.ndim); snapshot=snapshot.validate()
    i=int(candidate.source_packet)
    dq=candidate.q-state.q[i]; dp=candidate.p-state.p[i]
    gamma=state.width_matrices[i]
    invgamma=np.linalg.inv(gamma)
    width_eigs=np.sort(_safe_eigvals(candidate.width_matrices))
    chirp_eigs=np.sort(_safe_eigvals(candidate.chirp_matrices))
    metric_s=np.linalg.svd(0.5*(snapshot.metric+snapshot.metric.T),compute_uv=False)
    pop=float(np.vdot(state.coefficients[i],state.coefficients[i]).real)
    return {
        "feature_schema":V294_FEATURE_SCHEMA,
        "ngaussian":int(state.ngaussian), "ndim":int(state.ndim), "nstate":int(state.nstate),
        "source_population":pop,
        "dq_width_norm2":float(dq@gamma@dq),
        "dp_width_inverse_norm2":float(dp@invgamma@dp),
        "dq_euclidean_norm":float(np.linalg.norm(dq)),
        "dp_euclidean_norm":float(np.linalg.norm(dp)),
        "width_eigenvalues":width_eigs,
        "chirp_eigenvalues":chirp_eigs,
        "width_logdet":float(np.linalg.slogdet(candidate.width_matrices)[1]),
        "chirp_frobenius_norm":float(np.linalg.norm(candidate.chirp_matrices)),
        "tdvp_metric_rank":int(np.sum(metric_s>max(1e-12,1e-10*metric_s[0]))) if len(metric_s) else 0,
        "tdvp_metric_condition":float(metric_s[0]/metric_s[-1]) if len(metric_s) and metric_s[-1]>0 else 1e300,
        "tdse_residual_norm2":float(np.vdot(snapshot.residual_weighted,snapshot.residual_weighted).real),
        "displacement_kind":str(candidate.displacement_kind),
    }


@dataclass(frozen=True)
class SpawnOracleRecordV294:
    oracle_schema: str
    dataset_schema: str
    trajectory_group_id: str
    split: str
    state_fingerprint: str
    candidate_fingerprint: str
    oracle_config_fingerprint: str
    features: dict
    predicted_residual_reduction: float
    relative_residual_reduction: float
    schur_rank: int
    schur_condition_number: float
    schur_singular_values: tuple
    accepted_by_residual_oracle: bool
    rejection_reasons: tuple
    exact_candidate_cost_seconds: float | None = None

    def validate(self):
        if self.oracle_schema!=V294_ORACLE_SCHEMA or self.dataset_schema!=V294_DATASET_SCHEMA: raise ValueError("schema mismatch")
        if self.split not in ("train","validation","test"): raise ValueError("invalid split")
        if self.accepted_by_residual_oracle != (len(self.rejection_reasons)==0): raise ValueError("decision/reasons mismatch")
        for x in (self.predicted_residual_reduction,self.relative_residual_reduction,self.schur_condition_number):
            if not np.isfinite(float(x)): raise ValueError("non-finite oracle target")
        return self
    def as_dict(self): return _canon(asdict(self))
    @property
    def record_fingerprint(self): return stable_hash_v294(self.as_dict())


def evaluate_candidate_oracle_v294(state,connection,quadrature,snapshot,candidate,*,trajectory_group_id,packet_ids=None,path_factory=straight_path_v293,settings=CurvedBasisSettingsV292(),exact_candidate_cost_seconds=None):
    """Frozen exact spawning label API. Future ML may rank records but not redefine them."""
    state_payload=canonical_state_payload_v294(state,packet_ids=packet_ids,path_family=getattr(path_factory,"__name__","path"))
    ev=evaluate_curved_spawn_candidate_v293(state,connection,quadrature,snapshot,candidate,path_factory=path_factory,settings=settings)
    features=invariant_candidate_features_v294(state,candidate,snapshot)
    reasons=tuple(str(x) for x in ev.rejection_reasons)
    cond=float(ev.schur_condition_number)
    # Dataset records must be JSON finite; rejected rank-zero candidates get a large finite sentinel.
    if not np.isfinite(cond): cond=1e300
    return SpawnOracleRecordV294(
        V294_ORACLE_SCHEMA,V294_DATASET_SCHEMA,str(trajectory_group_id),trajectory_split_v294(trajectory_group_id),
        stable_hash_v294(state_payload),stable_hash_v294(candidate.as_dict()),stable_hash_v294(asdict(settings)),_canon(features),
        float(ev.predicted_residual_reduction),float(ev.relative_residual_reduction),int(ev.schur_rank),cond,
        tuple(float(x) for x in ev.schur_singular_values),bool(ev.admitted),reasons,
        None if exact_candidate_cost_seconds is None else float(exact_candidate_cost_seconds),
    ).validate()


def serialize_oracle_records_v294(records):
    """Deterministic JSONL, suitable for versioned datasets."""
    rows=[r.validate().as_dict() for r in records]
    rows.sort(key=lambda r:(r["trajectory_group_id"],r["state_fingerprint"],r["candidate_fingerprint"]))
    return "".join(json.dumps(r,sort_keys=True,separators=(",",":"))+"\n" for r in rows)


def dataset_manifest_v294(records):
    rows=[r.validate() for r in records]
    split_counts={s:sum(r.split==s for r in rows) for s in ("train","validation","test")}
    groups={s:sorted({r.trajectory_group_id for r in rows if r.split==s}) for s in split_counts}
    overlap=(set(groups["train"])&set(groups["validation"])) | (set(groups["train"])&set(groups["test"])) | (set(groups["validation"])&set(groups["test"]))
    if overlap: raise ValueError("trajectory leakage across splits")
    text=serialize_oracle_records_v294(rows)
    return _canon({"dataset_schema":V294_DATASET_SCHEMA,"feature_schema":V294_FEATURE_SCHEMA,"split_schema":V294_SPLIT_SCHEMA,
                   "n_records":len(rows),"split_counts":split_counts,"trajectory_groups":groups,"sha256":hashlib.sha256(text.encode()).hexdigest()})


def audit_exact_rejections_v294(records, *, audit_fraction=0.05):
    """Deterministically mark a subset of rejected rows for future exact ML false-negative audits."""
    if not 0<audit_fraction<=1: raise ValueError("audit_fraction must be in (0,1]")
    out=[]
    threshold=int(audit_fraction*(2**64-1))
    for r in records:
        if r.accepted_by_residual_oracle: continue
        h=int(hashlib.sha256(("audit|"+r.record_fingerprint).encode()).hexdigest()[:16],16)
        if h<=threshold: out.append(r.record_fingerprint)
    return tuple(sorted(out))
