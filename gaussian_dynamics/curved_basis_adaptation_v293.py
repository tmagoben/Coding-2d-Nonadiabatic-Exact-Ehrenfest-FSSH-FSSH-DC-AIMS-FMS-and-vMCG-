"""v0.29.3 curved TDVP/adaptation layer for coordinate-dependent connections.

This module preserves the v0.29.2 residual-driven Schur-complement spawning rule while
replacing the constant straight-section geometry with a coordinate-dependent connection
field and explicit path family.  It deliberately reuses the certified v0.27 correlated
Gaussian state and the v0.29.2 SVD/admission contracts.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from .correlated_gaussian_tdvp_v270 import (
    CorrelatedGaussianSpinorStateV270,
    pack_correlated_parameters_v270,
    state_from_correlated_parameters_v270,
)
from .correlated_basis_adaptation_v270 import (
    MultidimensionalSpawnCandidateV270,
    generate_multidimensional_spawn_candidates_v270,
)
from .curved_basis_adaptation_v292 import (
    CurvedQuadratureV292,
    CurvedBasisSettingsV292,
    CurvedTDVPSnapshotV292,
    CurvedCandidateEvaluationV292,
    CurvedProjectionReceiptV292,
    CurvedBasisEventV292,
    SVDReceiptV292,
    _gaussian_value,
    _weighted_flat,
    _svd_real_solve,
    _parameter_active_indices_v292,
    _evaluate_complete_coefficient_tangent_v292,
)
from .curved_geometry_v293 import ConnectionFieldV293, straight_path_v293, wilson_path_v293

V293_SCHEMA="gnd-coordinate-dependent-curved-basis-v0.29.3-dev"


def coefficient_basis_v293(state, connection: ConnectionFieldV293, quadrature: CurvedQuadratureV292, *, path_factory=straight_path_v293):
    state=state.validate(require_normalized=False); quadrature=quadrature.validate(); connection.validate()
    if (state.ndim,state.nstate)!=(connection.ndim,connection.nstate) or quadrature.ndim!=state.ndim:
        raise ValueError("state/connection/quadrature dimensions disagree")
    phi=np.empty((quadrature.npoint,state.nstate,state.ngaussian,state.nstate),complex)
    for i in range(state.ngaussian):
        g=_gaussian_value(state,i,quadrature.points)
        for n,r in enumerate(quadrature.points):
            w=wilson_path_v293(connection,path_factory(state.q[i],r))
            phi[n,:,i,:]=g[n]*w
    return phi


def evaluate_curved_state_v293(state,connection,quadrature,*,path_factory=straight_path_v293):
    phi=coefficient_basis_v293(state,connection,quadrature,path_factory=path_factory)
    return np.einsum('nsia,ia->ns',phi,state.coefficients,optimize=True)


def active_parameter_indices_v293(state, active_shape_mask=None):
    """Keep every coefficient and only enabled shape blocks, in packing order."""
    if active_shape_mask is None:
        return np.arange(state.parameter_count, dtype=int)
    mask = np.asarray(active_shape_mask)
    if mask.shape != (state.ngaussian,) or mask.dtype.kind != "b":
        raise ValueError("active shape mask must contain one Boolean per packet")
    return np.sort(_parameter_active_indices_v292(state, mask))


def finite_difference_full_jacobian_v293(state,connection,quadrature,*,path_factory=straight_path_v293,relative_step=2e-6,active_shape_mask=None):
    state=state.validate(require_normalized=False)
    theta=pack_correlated_parameters_v270(state); cols=[]
    for k in active_parameter_indices_v293(state, active_shape_mask):
        value = theta[k]
        h=relative_step*max(1.0,abs(float(value)))
        plus=theta.copy(); minus=theta.copy(); plus[k]+=h; minus[k]-=h
        sp=state_from_correlated_parameters_v270(plus,ngaussian=state.ngaussian,ndim=state.ndim,nstate=state.nstate,time_au=state.time_au)
        sm=state_from_correlated_parameters_v270(minus,ngaussian=state.ngaussian,ndim=state.ndim,nstate=state.nstate,time_au=state.time_au)
        fp=evaluate_curved_state_v293(sp,connection,quadrature,path_factory=path_factory)
        fm=evaluate_curved_state_v293(sm,connection,quadrature,path_factory=path_factory)
        cols.append(_weighted_flat((fp-fm)/(2*h),quadrature))
    return np.stack(cols,axis=1)


def build_tdvp_snapshot_v293(state,connection,quadrature,hpsi,*,path_factory=straight_path_v293,relative_step=2e-6,active_shape_mask=None):
    hpsi=np.asarray(hpsi,complex); quadrature=quadrature.validate()
    if hpsi.shape!=(quadrature.npoint,state.nstate): raise ValueError("hpsi has wrong shape")
    j=finite_difference_full_jacobian_v293(state,connection,quadrature,path_factory=path_factory,relative_step=relative_step,active_shape_mask=active_shape_mask)
    hw=_weighted_flat(hpsi,quadrature); m=np.real(j.conj().T@j); b=np.imag(j.conj().T@hw)
    v,receipt=_svd_real_solve(m,b); residual=j@v+1j*hw
    return CurvedTDVPSnapshotV292(j,m,b,v,residual,receipt).validate()


def candidate_coefficient_block_v293(candidate,nstate,connection,quadrature,*,path_factory=straight_path_v293):
    candidate=candidate.validate(quadrature.ndim)
    tmp=CorrelatedGaussianSpinorStateV270(candidate.q[None,:],candidate.p[None,:],candidate.width_matrices[None,:,:],candidate.chirp_matrices[None,:,:],np.ones((1,nstate),complex),0.0).validate(require_normalized=False)
    g=_gaussian_value(tmp,0,quadrature.points); cols=[]
    for a in range(nstate):
        e=np.zeros(nstate,complex); e[a]=1
        field=np.empty((quadrature.npoint,nstate),complex)
        for n,r in enumerate(quadrature.points):
            field[n]=g[n]*(wilson_path_v293(connection,path_factory(candidate.q,r))@e)
        cols.append(_weighted_flat(field,quadrature))
    return np.stack(cols,axis=1)


def evaluate_curved_spawn_candidate_v293(state,connection,quadrature,snapshot,candidate,*,path_factory=straight_path_v293,settings=CurvedBasisSettingsV292()):
    state = state.validate(require_normalized=False)
    connection = connection.validate()
    quadrature = quadrature.validate()
    k = candidate_coefficient_block_v293(
        candidate, state.nstate, connection, quadrature, path_factory=path_factory,
    )
    return _evaluate_complete_coefficient_tangent_v292(
        snapshot, candidate, k, settings=settings,
    )


def select_curved_spawn_candidate_v293(state,connection,quadrature,snapshot,*,path_factory=straight_path_v293,settings=CurvedBasisSettingsV292()):
    candidates=generate_multidimensional_spawn_candidates_v270(state,settings=settings.inherited)
    evals=[evaluate_curved_spawn_candidate_v293(state,connection,quadrature,snapshot,c,path_factory=path_factory,settings=settings) for c in candidates]
    admitted=[e for e in evals if e.admitted]
    if not admitted: return None,tuple(evals)
    admitted.sort(key=lambda e:(-e.predicted_residual_reduction,e.candidate.canonical_key()))
    return admitted[0],tuple(evals)


def project_curved_state_v293(source_state,target_state,connection,quadrature,*,path_factory=straight_path_v293,settings=CurvedBasisSettingsV292()):
    source_state=source_state.validate(require_normalized=False); target_state=target_state.validate(require_normalized=False); settings=settings.validate()
    psi=_weighted_flat(evaluate_curved_state_v293(source_state,connection,quadrature,path_factory=path_factory),quadrature)
    phi=coefficient_basis_v293(target_state,connection,quadrature,path_factory=path_factory)
    a=np.stack([_weighted_flat(phi[:,:,i,s],quadrature) for i in range(target_state.ngaussian) for s in range(target_state.nstate)],axis=1)
    u,sv,vh=np.linalg.svd(a,full_matrices=False); cut=max(settings.inherited.basis_absolute_cutoff,settings.inherited.basis_relative_cutoff*float(sv[0])); keep=sv>cut
    rank=int(np.count_nonzero(keep)); cond=float(sv[0]/sv[rank-1]) if rank else float('inf')
    if rank!=a.shape[1] or cond>settings.inherited.maximum_basis_condition_number: raise ValueError("target curved coefficient basis is rank deficient or ill conditioned")
    coeff=vh[keep].conj().T@((u[:,keep].conj().T@psi)/sv[keep]); projected=a@coeff
    sn=float(np.vdot(psi,psi).real); pn=float(np.vdot(projected,projected).real); ov=np.vdot(psi,projected)
    fidelity=float(abs(ov)**2/max(sn*pn,1e-30)); loss=float(np.vdot(psi-projected,psi-projected).real/max(sn,1e-30))
    out=CorrelatedGaussianSpinorStateV270(target_state.q,target_state.p,target_state.width_matrices,target_state.chirp_matrices,coeff.reshape(target_state.ngaussian,target_state.nstate),target_state.time_au).validate(require_normalized=False)
    return out,CurvedProjectionReceiptV292(sn,pn,fidelity,loss,sv,rank,cond)


def spawn_curved_candidate_v293(state,candidate,connection,quadrature,*,path_factory=straight_path_v293,settings=CurvedBasisSettingsV292()):
    state=state.validate(require_normalized=False); c=candidate.validate(state.ndim)
    target=CorrelatedGaussianSpinorStateV270(
        np.vstack([state.q,c.q]),np.vstack([state.p,c.p]),
        np.concatenate([state.width_matrices,c.width_matrices[None,:,:]]),
        np.concatenate([state.chirp_matrices,c.chirp_matrices[None,:,:]]),
        np.vstack([state.coefficients,np.zeros((1,state.nstate),complex)]),state.time_au).validate(require_normalized=False)
    out,receipt=project_curved_state_v293(state,target,connection,quadrature,path_factory=path_factory,settings=settings)
    if receipt.relative_loss>settings.inherited.maximum_projection_loss: raise ValueError("spawn projection loss exceeds gate")
    return out,receipt


def packet_subspace_overlap_v293(state,i,j,connection,quadrature,*,path_factory=straight_path_v293):
    phi=coefficient_basis_v293(state,connection,quadrature,path_factory=path_factory)
    ai=np.stack([_weighted_flat(phi[:,:,i,s],quadrature) for s in range(state.nstate)],axis=1)
    aj=np.stack([_weighted_flat(phi[:,:,j,s],quadrature) for s in range(state.nstate)],axis=1)
    block=ai.conj().T@aj
    return float(np.min(np.linalg.svd(block,compute_uv=False)))


def curved_energy_v293(state,connection,quadrature,hamiltonian_action,*,path_factory=straight_path_v293):
    psi=evaluate_curved_state_v293(state,connection,quadrature,path_factory=path_factory); hpsi=np.asarray(hamiltonian_action(state,psi,quadrature),complex)
    p=_weighted_flat(psi,quadrature); hp=_weighted_flat(hpsi,quadrature)
    return float((np.vdot(p,hp)/np.vdot(p,p)).real)


def _remove_target(state,remove):
    keep=np.array([i for i in range(state.ngaussian) if i!=remove],int)
    if not len(keep): raise ValueError("cannot remove final packet")
    target=CorrelatedGaussianSpinorStateV270(state.q[keep],state.p[keep],state.width_matrices[keep],state.chirp_matrices[keep],state.coefficients[keep],state.time_au).validate(require_normalized=False)
    return keep,target


def _project_removal(state,remove,connection,quadrature,hamiltonian_action,path_factory,settings):
    keep,target=_remove_target(state,remove)
    projected,receipt=project_curved_state_v293(state,target,connection,quadrature,path_factory=path_factory,settings=settings)
    de=curved_energy_v293(projected,connection,quadrature,hamiltonian_action,path_factory=path_factory)-curved_energy_v293(state,connection,quadrature,hamiltonian_action,path_factory=path_factory)
    return keep,projected,receipt,de


def curved_velocity_with_shape_mask_v293(state,connection,quadrature,hpsi,active_shape_mask,*,path_factory=straight_path_v293,relative_step=2e-6):
    jfull=finite_difference_full_jacobian_v293(state,connection,quadrature,path_factory=path_factory,relative_step=relative_step); idx=_parameter_active_indices_v292(state,active_shape_mask); j=jfull[:,idx]
    hw=_weighted_flat(np.asarray(hpsi,complex),quadrature); m=np.real(j.conj().T@j); b=np.imag(j.conj().T@hw); v,receipt=_svd_real_solve(m,b)
    full=np.zeros(jfull.shape[1]); full[idx]=v
    return full,receipt


def adapt_curved_basis_once_v293(state,connection,quadrature,snapshot,hamiltonian_action,*,path_factory=straight_path_v293,packet_ids=None,packet_ages=None,next_packet_serial=None,active_shape_mask=None,settings=CurvedBasisSettingsV292()):
    state=state.validate(require_normalized=False); settings=settings.validate(); snapshot=snapshot.validate(); quadrature=quadrature.validate(); connection.validate()
    if packet_ids is None: packet_ids=tuple(f"g{i:06d}" for i in range(state.ngaussian))
    if packet_ages is None: packet_ages=tuple(0 for _ in range(state.ngaussian))
    packet_ids=tuple(packet_ids); packet_ages=tuple(int(x) for x in packet_ages)
    if next_packet_serial is None: next_packet_serial=state.ngaussian
    if active_shape_mask is None: active_shape_mask=np.ones(state.ngaussian,bool)
    active_shape_mask=np.asarray(active_shape_mask,bool)
    def none(reason,evals=()):
        return CurvedBasisEventV292("none",reason,state,state,packet_ids,packet_ids,packet_ages,packet_ages,active_shape_mask,active_shape_mask,candidate_evaluations=tuple(evals)).validate()
    for i in range(state.ngaussian):
        for j in range(i+1,state.ngaussian):
            if packet_subspace_overlap_v293(state,i,j,connection,quadrature,path_factory=path_factory)<settings.inherited.minimum_merge_overlap: continue
            pi=float(np.sum(abs(state.coefficients[i])**2)); pj=float(np.sum(abs(state.coefficients[j])**2)); remove=j if pj<=pi else i
            try: keep,out,receipt,de=_project_removal(state,remove,connection,quadrature,hamiltonian_action,path_factory,settings)
            except ValueError: continue
            if receipt.relative_loss<=settings.inherited.maximum_projection_loss and abs(de)<=settings.inherited.maximum_event_energy_jump_hartree:
                return CurvedBasisEventV292("merge","coordinate-dependent curved overlap and projection gates passed",state,out,packet_ids,tuple(packet_ids[k] for k in keep),packet_ages,tuple(packet_ages[k] for k in keep),active_shape_mask,active_shape_mask[keep],projection=receipt,energy_jump_hartree=de,removed_packet_id=packet_ids[remove]).validate()
    order=sorted(range(state.ngaussian),key=lambda i:(float(np.sum(abs(state.coefficients[i])**2)),packet_ids[i]))
    for remove in order:
        if state.ngaussian<=1: break
        if packet_ages[remove]<settings.inherited.minimum_prune_age or float(np.sum(abs(state.coefficients[remove])**2))>settings.inherited.prune_population_threshold: continue
        try: keep,out,receipt,de=_project_removal(state,remove,connection,quadrature,hamiltonian_action,path_factory,settings)
        except ValueError: continue
        if receipt.relative_loss<=settings.inherited.maximum_projection_loss and abs(de)<=settings.inherited.maximum_event_energy_jump_hartree:
            return CurvedBasisEventV292("prune","coordinate-dependent curved projection gates passed",state,out,packet_ids,tuple(packet_ids[k] for k in keep),packet_ages,tuple(packet_ages[k] for k in keep),active_shape_mask,active_shape_mask[keep],projection=receipt,energy_jump_hartree=de,removed_packet_id=packet_ids[remove]).validate()
    if np.any(~active_shape_mask): return none("dormant packet activation pending")
    if state.ngaussian>=settings.inherited.maximum_packet_count: return none("maximum packet count reached")
    selected,evals=select_curved_spawn_candidate_v293(state,connection,quadrature,snapshot,path_factory=path_factory,settings=settings)
    if selected is None: return none("no coordinate-dependent curved spawn candidate passed every gate",evals)
    try: out,receipt=spawn_curved_candidate_v293(state,selected.candidate,connection,quadrature,path_factory=path_factory,settings=settings)
    except ValueError: return none("selected coordinate-dependent curved candidate failed projection",evals)
    added=f"g{int(next_packet_serial):06d}"
    return CurvedBasisEventV292("spawn","largest coordinate-dependent full-tangent residual reduction passed",state,out,packet_ids,packet_ids+(added,),packet_ages,packet_ages+(0,),active_shape_mask,np.concatenate([active_shape_mask,[False]]),projection=receipt,selected_candidate=selected.candidate,candidate_evaluations=evals,added_packet_id=added).validate()
