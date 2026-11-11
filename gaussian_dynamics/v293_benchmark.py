"""Deterministic scientific evidence for the v0.29.3 development checkpoint."""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import numpy as np
from scipy.linalg import expm

from .correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from .correlated_basis_adaptation_v270 import generate_multidimensional_spawn_candidates_v270
from .curved_basis_adaptation_v292 import ConstantConnectionV292, CurvedBasisSettingsV292, CurvedTDVPSnapshotV292, SVDReceiptV292, actual_linear_residual_reduction_v292, coefficient_basis_v292, tensor_legendre_quadrature_v292
from .curved_geometry_v293 import ConnectionFieldV293, CoordinateGaugeV293, PatchTransitionV293, PolylinePathV293, constant_connection_field_v293, gauge_transform_connection_v293, straight_path_v293, dogleg_path_v293, wilson_path_v293, plaquette_holonomy_v293, validate_patch_transition_v293, two_patch_transport_v293
from .curved_basis_adaptation_v293 import coefficient_basis_v293, finite_difference_full_jacobian_v293, candidate_coefficient_block_v293, evaluate_curved_spawn_candidate_v293, spawn_curved_candidate_v293, adapt_curved_basis_once_v293

SCHEMA="gnd-v0.29.3-coordinate-dependent-geometry-evidence-v1"

def _canon(x):
    if isinstance(x,np.generic): return _canon(x.item())
    if isinstance(x,np.ndarray): return _canon(x.tolist())
    if isinstance(x,complex): return [float(x.real),float(x.imag)]
    if isinstance(x,dict): return {str(k):_canon(v) for k,v in sorted(x.items())}
    if isinstance(x,(list,tuple)): return [_canon(v) for v in x]
    if isinstance(x,float):
        if not np.isfinite(x): raise ValueError("nonfinite evidence")
        return float(x)
    if isinstance(x,(str,int,bool)) or x is None: return x
    raise TypeError(type(x).__name__)

def _fingerprint(payload):
    return hashlib.sha256(json.dumps(_canon(payload),sort_keys=True,separators=(",",":")).encode()).hexdigest()

def _paulis():
    sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]],complex); sz=np.array([[1,0],[0,-1]],complex)
    return sx,sy,sz

def _state():
    return CorrelatedGaussianSpinorStateV270(np.array([[-.25,.15]]),np.array([[.12,-.08]]),np.array([[[1.05,.14],[.14,.86]]]),np.array([[[.05,-.02],[-.02,.03]]]),np.array([[.84+.08j,.43-.14j]])).normalized()

def _constant():
    _,sy,sz=_paulis(); g=np.array([.19j*sy,.13j*sz]); return g,constant_connection_field_v293(g)

def _gauge():
    _,_,sz=_paulis(); x=.35j*sz
    def th(r): return .18*r[0]-.11*r[1]+.07*r[0]*r[1]
    def gr(r): return np.array([.18+.07*r[1],-.11+.07*r[0]])
    def mat(r): return expm(th(r)*x)
    def der(r):
        gg=mat(r); return np.asarray([gr(r)[m]*x@gg for m in range(2)])
    return CoordinateGaugeV293(2,2,mat,der,"benchmark-gauge")

def _variable():
    sx,sy,sz=_paulis()
    def ev(r):
        x,y=r
        return np.array([.14j*sy+.055j*y*sz,.11j*sx-.045j*x*sz])
    return ConnectionFieldV293(2,2,ev,"benchmark-variable").validate([np.zeros(2)])

def _snapshot(state,field,quad,candidate):
    j=finite_difference_full_jacobian_v293(state,field,quad,relative_step=4e-6); k=candidate_coefficient_block_v293(candidate,state.nstate,field,quad)
    residual=.13*k[:,0]-.09*k[:,1]+.025j*k[:,0]
    jr=np.vstack([j.real,j.imag]); rr=np.concatenate([residual.real,residual.imag]); rr=rr+jr@np.linalg.lstsq(jr,-rr,rcond=1e-10)[0]
    n=j.shape[0]; residual=rr[:n]+1j*rr[n:]; m=np.real(j.conj().T@j); s=np.linalg.svd(m,compute_uv=False)
    return CurvedTDVPSnapshotV292(j,m,np.zeros(j.shape[1]),np.zeros(j.shape[1]),residual,SVDReceiptV292(s,1e-12,len(s),float(s[0]/s[-1]),0,0)).validate()

def run_v293_evidence():
    state=_state(); gens,const=_constant(); old=ConstantConnectionV292(gens).validate(); q=tensor_legendre_quadrature_v292(((-2,2),(-2,2)),(5,5))
    reduction=float(np.max(np.abs(coefficient_basis_v292(state,old,q)-coefficient_basis_v293(state,const,q))))
    gauge=_gauge(); target=gauge_transform_connection_v293(const,gauge); a=np.array([-.4,.2]); b=np.array([.55,-.35])
    w=wilson_path_v293(const,straight_path_v293(a,b)); wp=wilson_path_v293(target,straight_path_v293(a,b)); gauge_err=float(np.linalg.norm(wp-gauge.at(b).conj().T@w@gauge.at(a)))
    variable=_variable(); a2=np.array([-.35,-.25]); b2=np.array([.45,.4]); ws=wilson_path_v293(variable,straight_path_v293(a2,b2)); wd=wilson_path_v293(variable,dogleg_path_v293(a2,b2,axis=0)); path_sep=float(np.linalg.norm(ws-wd))
    p=dogleg_path_v293(np.array([-.2,.1]),np.array([.5,-.3]),axis=1); reverse_err=float(np.linalg.norm(wilson_path_v293(variable,p.reversed())@wilson_path_v293(variable,p)-np.eye(2)))
    r=np.array([.13,-.21]); h=1.5e-3; plaq=plaquette_holonomy_v293(variable,r,0,1,h); f=variable.curvature(r,0,1,step=2e-6); curv_err=float(np.linalg.norm((plaq-np.eye(2))/(h*h)+f)/np.linalg.norm(f))
    gt=_gauge(); patch_target=gauge_transform_connection_v293(variable,gt); trans=PatchTransitionV293("alpha","beta",gt); probes=[np.array([-.2,.1]),np.array([.1,-.15]),np.array([.4,.2])]; patch_local=float(validate_patch_transition_v293(variable,patch_target,trans,probes))
    pa=np.array([-.45,-.1]); px=np.array([.02,.18]); pb=np.array([.42,.35]); mixed=two_patch_transport_v293(variable,patch_target,trans,pa,px,pb); same=PolylinePathV293(np.vstack([pa,px,pb]),"same").validate(); patch_transport=float(np.linalg.norm(mixed-gt.at(pb).conj().T@wilson_path_v293(variable,same)))
    q2=tensor_legendre_quadrature_v292(((-1.8,1.8),(-1.8,1.8)),(4,4)); candidate=generate_multidimensional_spawn_candidates_v270(state)[0]; snap=_snapshot(state,variable,q2,candidate); settings=CurvedBasisSettingsV292(minimum_predicted_residual_reduction=1e-14,minimum_relative_residual_reduction=1e-12); ev=evaluate_curved_spawn_candidate_v293(state,variable,q2,snap,candidate,settings=settings); k=candidate_coefficient_block_v293(candidate,state.nstate,variable,q2); actual=actual_linear_residual_reduction_v292(snap,k); gain_rel=float(abs(ev.predicted_residual_reduction-actual)/actual)
    spawned,receipt=spawn_curved_candidate_v293(state,candidate,variable,q2,settings=settings)
    def ih(_s,psi,_q): return np.asarray(psi,complex)
    event=adapt_curved_basis_once_v293(state,variable,q2,snap,ih,settings=settings)
    metrics={"constant_reduction_max_abs":reduction,"wilson_gauge_covariance_error":gauge_err,"straight_dogleg_separation":path_sep,"reverse_path_inverse_error":reverse_err,"plaquette_curvature_relative_error":curv_err,"patch_connection_relation_error":patch_local,"patch_same_path_transport_error":patch_transport,"spawn_gain_prediction_relative_error":gain_rel,"spawn_projection_relative_loss":float(receipt.relative_loss),"spawn_fidelity":float(receipt.fidelity),"newborn_coefficient_norm":float(np.linalg.norm(spawned.coefficients[-1])),"spawn_event_kind":event.event_kind,"newborn_shape_dormant":bool(not event.active_shape_mask_after[-1])}
    gates={"constant_reduction":reduction<1e-12,"gauge_covariance":gauge_err<1e-10,"path_dependence_detected":path_sep>1e-4,"reverse_transport":reverse_err<1e-9,"curvature_holonomy":curv_err<2e-2,"patch_relation":patch_local<1e-6,"patch_transport":patch_transport<1e-8,"spawn_gain":gain_rel<5e-6,"projection":receipt.relative_loss<1e-18 and receipt.fidelity>1-1e-11,"newborn":np.linalg.norm(spawned.coefficients[-1])<1e-8 and event.event_kind=="spawn" and not event.active_shape_mask_after[-1]}
    payload={"schema":SCHEMA,"metrics":metrics,"gates":gates,"passed":bool(all(gates.values()))}; payload["fingerprint"]=_fingerprint(payload)
    return payload

def save_v293_evidence(path):
    payload=run_v293_evidence(); Path(path).write_text(json.dumps(_canon(payload),indent=2,sort_keys=True)+"\n"); return payload
