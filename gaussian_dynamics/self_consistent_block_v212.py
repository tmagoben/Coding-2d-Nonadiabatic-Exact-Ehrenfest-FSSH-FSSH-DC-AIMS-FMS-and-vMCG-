from dataclasses import dataclass, asdict
import numpy as np
from scipy import sparse

from .block_sparse_molecular_v21 import (
    BlockMolecularTBFV21, BlockSparseSettingsV21, BlockSparseMolecularGraphV21,
    build_block_sparse_matrices_v21, build_dense_block_reference_v21,
)
from .block_basis_lifecycle_v212 import insert_zero_block_v212, prune_block_projected_v212
from .sparse_pair_matrices_v16 import (
    sparse_metric_compatible_connection, sparse_moving_basis_midpoint_cayley_step,
    sparse_generalized_norm,
)


@dataclass(frozen=True)
class MeanFieldGuidanceSettingsV212:
    minimum_local_amplitude: float=1e-12
    low_amplitude_policy: str="zero_force"

    def validate(self):
        if not np.isfinite(self.minimum_local_amplitude) or self.minimum_local_amplitude<0.0:
            raise ValueError("minimum_local_amplitude cannot be negative.")
        if self.low_amplitude_policy != "zero_force":
            raise ValueError(
                "v0.21.3 removes the gauge-dependent lowest-eigenvector fallback; "
                "use zero_force or the density-matrix v0.21.3 guidance engine."
            )
        return self


class BlockMeanFieldGuidanceV212:
    """Representation-neutral local electronic mean-field guidance.

    For Gaussian i, normalize its local electronic coefficient block c_i and use

        F_i,a = - c_i^dagger K_a(q_i) c_i.

    The rule is gauge invariant under c_i' = G_i^dagger c_i and K_a'=G_i^dagger K_aG_i.
    It is a nuclear *guidance policy*, not a claim of full variational AIMS nuclear
    equations.
    """
    def __init__(self,settings=MeanFieldGuidanceSettingsV212()):
        self.settings=settings.validate()

    def forces_and_masses(self,basis,coefficients,provider,nstate):
        basis=list(basis)
        C=np.asarray(coefficients,dtype=complex)
        s=int(nstate)
        if C.shape!=(len(basis)*s,):
            raise ValueError("coefficient dimension is inconsistent with basis*nstate.")
        forces=[]; masses=[]; local_norms=[]
        for i,b in enumerate(basis):
            snap=provider.evaluate_snapshot(b.q)
            c=C[s*i:s*(i+1)]
            amp=float(np.real(np.vdot(c,c)))
            local_norms.append(amp)
            if amp<=self.settings.minimum_local_amplitude:
                f=np.zeros(snap.point.nq,dtype=float)
            else:
                f=snap.point.force_expectation(c)
            forces.append(f)
            masses.append(np.asarray(snap.point.mass_matrix_q_au,float))
        return np.asarray(forces,float),tuple(masses),np.asarray(local_norms,float)


@dataclass(frozen=True)
class SelfConsistentBlockSettingsV212:
    graph: BlockSparseSettingsV21=BlockSparseSettingsV21()
    guidance: MeanFieldGuidanceSettingsV212=MeanFieldGuidanceSettingsV212()
    use_dense_reference: bool=True
    corrector_iterations: int=2
    momentum_tolerance: float=1e-10

    def validate(self):
        self.graph.validate(); self.guidance.validate()
        if self.corrector_iterations<1 or int(self.corrector_iterations)!=self.corrector_iterations:
            raise ValueError("corrector_iterations must be >=1.")
        if not np.isfinite(self.momentum_tolerance) or self.momentum_tolerance<=0.0:
            raise ValueError("momentum_tolerance must be positive.")
        return self


def _copy_basis(basis):
    return [b.copy() for b in basis]


def _matrices(basis,provider,dt,qdots,pdots,settings,graph=None):
    if settings.use_dense_reference:
        d=build_dense_block_reference_v21(basis,provider,dt,qdots,pdots,settings.graph)
        return {
            "S":sparse.csr_matrix(d["S"]),"H":sparse.csr_matrix(d["H"]),
            "T":sparse.csr_matrix(d["T_seed"]),"update":None,
        }
    u=graph.update(basis,qdots,pdots)
    m=build_block_sparse_matrices_v21(basis,u)
    return {"S":m.S,"H":m.H,"T":m.T_seed,"update":u}


def _velocities(basis,masses):
    return np.asarray([np.linalg.solve(M,b.p) for b,M in zip(basis,masses)],float)


def _metric_phase_error(a,b,S):
    z=np.vdot(a,S@b)
    phase=1.0+0j if abs(z)<1e-30 else np.exp(-1j*np.angle(z))
    d=phase*b-a
    return float(np.sqrt(max(np.real(np.vdot(d,S@d)),0.0))/np.sqrt(max(np.real(np.vdot(a,S@a)),1e-30)))


def run_self_consistent_block_dynamics_v212(
    initial_basis,C0,provider,*,dt=0.002,steps=20,
    settings=SelfConsistentBlockSettingsV212(),store_every=5,adaptation_policy=None,
    guidance_engine=None,graph_active_uid_edges=None,
):
    settings=settings.validate(); dt=float(dt)
    if int(steps)!=steps:
        raise ValueError("steps must be an integer.")
    steps=int(steps)
    if int(store_every)!=store_every or int(store_every)<1:
        raise ValueError("store_every must be a positive integer.")
    store_every=int(store_every)
    if not np.isfinite(dt) or dt<=0.0 or steps<0:
        raise ValueError("invalid propagation interval.")
    basis=_copy_basis(initial_basis)
    if not basis:
        raise ValueError("basis cannot be empty.")
    s=provider.evaluate_snapshot(basis[0].q).point.nstate
    C=np.asarray(C0,dtype=complex).copy()
    if C.shape!=(len(basis)*s,):
        raise ValueError("C0 has incompatible dimension.")
    if not np.all(np.isfinite(C)):
        raise ValueError("C0 contains non-finite data.")

    guidance=(
        BlockMeanFieldGuidanceV212(settings.guidance)
        if guidance_engine is None
        else guidance_engine
    )
    graph=None if settings.use_dense_reference else BlockSparseMolecularGraphV21(provider,dt,settings.graph)
    if graph is None:
        if graph_active_uid_edges not in (None, (), []):
            raise ValueError("dense propagation cannot restore sparse graph edges.")
    elif graph_active_uid_edges is not None:
        graph.restore_active_uid_edges_v214(
            graph_active_uid_edges,[item.uid for item in basis]
        )

    F,masses,local_norms=guidance.forces_and_masses(basis,C,provider,s)
    qdot=_velocities(basis,masses)
    current=_matrices(basis,provider,dt,qdot,F,settings,graph)
    C=C/np.sqrt(sparse_generalized_norm(C,current["S"]))
    # recompute force after metric normalization (direction unchanged locally)
    F,masses,local_norms=guidance.forces_and_masses(basis,C,provider,s)

    records=[]; corrector_history=[]; adaptation_events=[]
    guidance_trial_state_rollbacks=0
    def record(step):
        records.append({
            "step":int(step),"time":float(step*dt),
            "norm":float(sparse_generalized_norm(C,current["S"])),
            "condition_number":float(np.linalg.cond(current["S"].toarray())),
            "maximum_force":float(np.max(np.linalg.norm(F,axis=1))),
            "minimum_local_amplitude":float(np.min(local_norms)),
            "active_edges":None if current["update"] is None else int(len(current["update"].active_edges)),
        })
    record(0)

    for step in range(1,steps+1):
        old_basis=_copy_basis(basis); old_C=C.copy(); old=current
        accepted_guidance_state=(
            guidance.checkpoint_state()
            if hasattr(guidance,"checkpoint_state") and hasattr(guidance,"restore_state")
            else None
        )
        F0=F.copy(); masses0=masses
        p_half=[]; q_new=[]
        for b,f,M in zip(old_basis,F0,masses0):
            ph=b.p+0.5*dt*f
            p_half.append(ph)
            q_new.append(b.q+dt*np.linalg.solve(M,ph))
        p_half=np.asarray(p_half,float); q_new=np.asarray(q_new,float)

        # Explicit predictor for the endpoint momentum.
        p_trial=p_half+0.5*dt*F0
        F_trial=F0.copy(); C_trial=old_C.copy(); endpoint=None
        converged=False
        for it in range(settings.corrector_iterations):
            if accepted_guidance_state is not None:
                guidance.restore_state(accepted_guidance_state)
                guidance_trial_state_rollbacks+=1
            trial_basis=[BlockMolecularTBFV21(b.uid,q_new[i],p_trial[i],b.A.copy()) for i,b in enumerate(old_basis)]
            # Require the generalized mass matrix to remain constant over one step.
            snaps=[provider.evaluate_snapshot(b.q) for b in trial_basis]
            masses1=tuple(np.asarray(x.point.mass_matrix_q_au,float) for x in snaps)
            for M0,M1 in zip(masses0,masses1):
                if not np.allclose(M0,M1,rtol=1e-10,atol=1e-12):
                    raise ValueError("v0.21.2 self-consistent Verlet currently requires a constant generalized mass matrix over each step.")
            qdot1=_velocities(trial_basis,masses1)
            endpoint=_matrices(trial_basis,provider,dt,qdot1,F_trial,settings,graph)
            seed=0.5*(old["T"]+endpoint["T"])
            Tmid=sparse_metric_compatible_connection(old["S"],endpoint["S"],dt,seed)
            C_trial=sparse_moving_basis_midpoint_cayley_step(
                old_C,old["S"],old["H"],endpoint["S"],endpoint["H"],Tmid,dt
            )
            F1,_,local1=guidance.forces_and_masses(trial_basis,C_trial,provider,s)
            p_new=p_half+0.5*dt*F1
            delta=float(np.max(np.linalg.norm(p_new-p_trial,axis=1)))
            corrector_history.append({"step":int(step),"iteration":int(it+1),"momentum_change":delta})
            p_trial=p_new; F_trial=F1
            if delta<=settings.momentum_tolerance:
                converged=True
                break

        # Final consistent endpoint rebuild at the corrected momentum/force.  Trial
        # guide densities are discarded; only this accepted endpoint is committed.
        if accepted_guidance_state is not None:
            guidance.restore_state(accepted_guidance_state)
            guidance_trial_state_rollbacks+=1
        basis=[BlockMolecularTBFV21(b.uid,q_new[i],p_trial[i],b.A.copy()) for i,b in enumerate(old_basis)]
        masses=tuple(np.asarray(provider.evaluate_snapshot(b.q).point.mass_matrix_q_au,float) for b in basis)
        qdot1=_velocities(basis,masses)
        current=_matrices(basis,provider,dt,qdot1,F_trial,settings,graph)
        seed=0.5*(old["T"]+current["T"])
        Tmid=sparse_metric_compatible_connection(old["S"],current["S"],dt,seed)
        C=sparse_moving_basis_midpoint_cayley_step(
            old_C,old["S"],old["H"],current["S"],current["H"],Tmid,dt
        )
        F,masses,local_norms=guidance.forces_and_masses(basis,C,provider,s)

        # Optional representation-neutral adaptive basis action at a converged step boundary.
        if adaptation_policy is not None:
            action=adaptation_policy(step,tuple(basis),C.copy(),current["S"])
            actions=[] if action is None else (list(action) if isinstance(action,(list,tuple)) else [action])
            for item in actions:
                if "insert" in item:
                    basis_tuple,C=insert_zero_block_v212(basis,C,item["insert"],s,item.get("index"))
                    basis=list(basis_tuple)
                    if hasattr(guidance,"on_insert"):
                        guidance.on_insert(
                            item["insert"],
                            provider,
                            parent_uid=item.get("guide_parent_uid"),
                            guide_density=item.get("guide_density"),
                        )
                    adaptation_events.append({"step":int(step),"kind":"insert","uid":int(item["insert"].uid),"basis_size":len(basis)})
                elif "prune_index" in item:
                    result=prune_block_projected_v212(basis,C,current["S"],s,item["prune_index"])
                    basis=list(result.basis); C=result.coefficients
                    if hasattr(guidance,"on_prune"):
                        guidance.on_prune(result.removed_uid)
                    adaptation_events.append({"step":int(step),"kind":"prune","uid":int(result.removed_uid),"projection_loss":float(result.projection_loss),"basis_size":len(basis)})
                else:
                    raise ValueError("adaptation action must contain 'insert' or 'prune_index'.")

            if actions:
                F,masses,local_norms=guidance.forces_and_masses(basis,C,provider,s)
                qdot=_velocities(basis,masses)
                current=_matrices(basis,provider,dt,qdot,F,settings,graph)

        if step%store_every==0:
            record(step)

    output={
        "records":records,"corrector_history":corrector_history,"adaptation_events":adaptation_events,
        "final_basis":basis,"final_coefficients":C,"final_S":current["S"],
        "final_H":current["H"],"final_T_seed":current["T"],
        "settings":{"dt":dt,"steps":steps,"control":asdict(settings)},
        "maximum_norm_drift":float(max(abs(r["norm"]-1.0) for r in records)),
        "guidance_trial_state_rollbacks":int(guidance_trial_state_rollbacks),
        "final_active_uid_edges":(
            () if graph is None else graph.active_uid_edges_v214()
        ),
    }
    if hasattr(guidance,"diagnostics_dict"):
        output["guidance_diagnostics"]=guidance.diagnostics_dict()
    return output
