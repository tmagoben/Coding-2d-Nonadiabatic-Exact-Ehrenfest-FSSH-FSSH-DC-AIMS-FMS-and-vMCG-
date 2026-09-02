from dataclasses import dataclass, asdict, replace
import numpy as np
from scipy import sparse

from .direct_dynamics_nd import (
    maybe_spawn_once,
)
from .local_gaussian_nd import (
    LocalAdiabaticTBF,
)
from .sparse_pair_matrices_v16 import (
    sparse_metric_compatible_connection,
    sparse_moving_basis_midpoint_cayley_step,
    sparse_generalized_norm,
)
from .sparse_molecular_matrices_v20 import (
    SparseMolecularTBFV20,
    MolecularSparseSettingsV20,
    SparseMolecularEdgeGraphV20,
    build_sparse_molecular_matrices_v20,
    build_dense_molecular_reference_v20,
    relative_frobenius_error,
)
from .sampled_molecular_audit_v20 import (
    sampled_molecular_edge_audit_v20,
)


@dataclass(frozen=True)
class SparseMolecularDynamicsSettingsV20:
    graph: MolecularSparseSettingsV20=MolecularSparseSettingsV20()
    sampled_audit_interval: int=10
    sampled_audit_priority_pairs: int=6
    sampled_audit_random_pairs: int=6
    sampled_audit_search_factor: float=0.1
    sampled_audit_seed: int=20260813
    sampled_audit_violation_factor: float=1.0
    sampled_audit_relaxation_factor: float=0.5
    max_sampled_audit_relaxations: int=3

    dense_sentinel_S_limit: float=0.01
    dense_sentinel_H_limit: float=0.01
    dense_sentinel_T_limit: float=0.02

    allow_spawning: bool=False
    spawn_threshold: float=1e-3
    overlap_block: float=0.999999
    max_basis: int=32

    def validate(self):
        self.graph.validate()
        if self.sampled_audit_interval<=0:
            raise ValueError("sampled_audit_interval must be positive.")
        if self.sampled_audit_priority_pairs<0 or self.sampled_audit_random_pairs<0:
            raise ValueError("audit counts cannot be negative.")
        if not (0.0<self.sampled_audit_search_factor<=1.0):
            raise ValueError("sampled_audit_search_factor must lie in (0,1].")
        if not (0.0<self.sampled_audit_relaxation_factor<1.0):
            raise ValueError("sampled_audit_relaxation_factor must lie in (0,1).")
        if self.max_sampled_audit_relaxations<0:
            raise ValueError("max_sampled_audit_relaxations cannot be negative.")
        if self.max_basis<=0:
            raise ValueError("max_basis must be positive.")
        return self


def _copy_basis(basis):
    out=[]
    for i,b in enumerate(basis):
        if isinstance(b,SparseMolecularTBFV20):
            out.append(b.copy())
        else:
            out.append(
                SparseMolecularTBFV20(
                    uid=int(getattr(b,"uid",i)),
                    state=int(b.state),
                    q=np.asarray(b.q,float).copy(),
                    p=np.asarray(b.p,float).copy(),
                    A=np.asarray(b.A,float).copy(),
                )
            )
    return out


def _to_local(basis):
    return [
        LocalAdiabaticTBF(
            b.state,b.q.copy(),b.p.copy(),b.A.copy()
        )
        for b in basis
    ]


def _nuclear_velocity_verlet_step(
    basis,
    provider,
    dt,
):
    new=_copy_basis(basis)
    for old,b in zip(basis,new):
        point0=provider.evaluate(old.q)
        M0=point0.mass_matrix_q_au
        force0=-point0.gradients_q[old.state]

        p_half=old.p+0.5*dt*force0
        q_new=old.q+dt*np.linalg.solve(
            M0,p_half
        )

        point1=provider.evaluate(q_new)
        if not np.allclose(
            point1.mass_matrix_q_au,
            M0,atol=1e-10,rtol=1e-10
        ):
            raise ValueError(
                "v0.20 velocity-Verlet runner currently requires a constant generalized mass matrix."
            )
        force1=-point1.gradients_q[old.state]
        p_new=p_half+0.5*dt*force1

        b.q=q_new
        b.p=p_new
    return new


def _dense_sentinel(
    basis,
    provider,
    dt,
    sparse_mats,
    settings,
):
    dense=build_dense_molecular_reference_v20(
        basis,provider,dt,settings.graph
    )
    row={
        "relative_S_frobenius_error":
            relative_frobenius_error(
                sparse_mats.S.toarray(),
                dense["S"],
            ),
        "relative_H_frobenius_error":
            relative_frobenius_error(
                sparse_mats.H.toarray(),
                dense["H"],
            ),
        "relative_Tseed_frobenius_error":
            relative_frobenius_error(
                sparse_mats.T_seed.toarray(),
                dense["T_seed"],
            ),
        "dense_pair_count":
            int(dense["pair_count"]),
        "active_edge_count":
            int(len(sparse_mats.active_edges)),
    }
    row["passed"]=bool(
        row["relative_S_frobenius_error"]
        <=settings.dense_sentinel_S_limit
        and row["relative_H_frobenius_error"]
        <=settings.dense_sentinel_H_limit
        and row["relative_Tseed_frobenius_error"]
        <=settings.dense_sentinel_T_limit
    )
    return row


def run_sparse_molecular_dynamics_v20(
    initial_basis,
    C0,
    provider,
    *,
    dt=0.01,
    steps=50,
    settings=SparseMolecularDynamicsSettingsV20(),
    store_every=5,
    sentinel_provider_factory=None,
):
    """End-to-end sparse molecular Gaussian propagation.

    Nuclear centers follow independent active-surface velocity-Verlet guidance.
    Coefficients use sparse moving-basis midpoint/Cayley propagation.  Molecular
    pair-centroid electronic calculations are created only for geometrically local
    S/H/T candidates, active edges, and fixed-size sampled audits.
    """
    settings=settings.validate()
    dt=float(dt)
    if dt<=0.0:
        raise ValueError("dt must be positive.")
    steps=int(steps)
    basis=_copy_basis(initial_basis)
    C=np.asarray(C0,dtype=complex).copy()
    if C.shape!=(len(basis),):
        raise ValueError("C0 must contain one coefficient per TBF.")

    graph=SparseMolecularEdgeGraphV20(
        provider,dt,settings.graph
    )
    update=graph.update(basis)
    mats=build_sparse_molecular_matrices_v20(
        basis,update
    )

    norm0=sparse_generalized_norm(
        C,mats.S
    )
    if norm0<=0.0:
        raise ValueError("initial generalized norm must be positive.")
    C=C/np.sqrt(norm0)

    sentinel_provider_initial=(
        sentinel_provider_factory()
        if sentinel_provider_factory is not None
        else provider
    )
    initial_sentinel=_dense_sentinel(
        basis,
        sentinel_provider_initial,
        dt,mats,settings
    )
    if not initial_sentinel["passed"]:
        raise RuntimeError(
            "Initial v0.20 dense sentinel failed."
        )

    records=[]
    sampled_audits=[]
    events=[]

    def record(step,current_mats,current_update):
        records.append({
            "step":int(step),
            "time":float(step*dt),
            "basis_size":int(len(basis)),
            "norm":sparse_generalized_norm(
                C,current_mats.S
            ),
            "condition_number":float(
                np.linalg.cond(
                    current_mats.S.toarray()
                )
            ),
            "active_edges":
                int(len(current_mats.active_edges)),
            "total_pairs":
                int(current_update.total_offdiagonal_pairs),
            "sparsity_fraction":
                float(current_update.sparsity_fraction),
            "exact_pair_checks":
                int(current_update.exact_pair_checks),
            "omitted_score_l2":
                float(current_update.omitted_score_l2),
        })

    record(0,mats,update)

    for step in range(1,steps+1):
        new_basis=_nuclear_velocity_verlet_step(
            basis,provider,dt
        )
        new_update=graph.update(new_basis)
        new_mats=build_sparse_molecular_matrices_v20(
            new_basis,new_update
        )

        seed=0.5*(
            mats.T_seed+new_mats.T_seed
        )
        Tmid=sparse_metric_compatible_connection(
            mats.S,new_mats.S,dt,seed
        )
        C=sparse_moving_basis_midpoint_cayley_step(
            C,
            mats.S,mats.H,
            new_mats.S,new_mats.H,
            Tmid,dt,
        )
        basis=new_basis
        update=new_update
        mats=new_mats

        if (
            settings.allow_spawning
            and len(basis)<settings.max_basis
        ):
            local=_to_local(basis)
            parent,child=maybe_spawn_once(
                local,provider,
                threshold=settings.spawn_threshold,
                overlap_block=settings.overlap_block,
            )
            if child is not None:
                next_uid=max(
                    b.uid for b in basis
                )+1
                basis.append(
                    SparseMolecularTBFV20(
                        next_uid,child.state,
                        child.q.copy(),
                        child.p.copy(),
                        child.A.copy(),
                    )
                )
                C=np.concatenate([
                    C,np.array([0.0+0.0j])
                ])
                update=graph.update(basis)
                mats=build_sparse_molecular_matrices_v20(
                    basis,update
                )
                events.append({
                    "step":int(step),
                    "time":float(step*dt),
                    "parent_index":int(parent),
                    "new_index":int(len(basis)-1),
                    "new_uid":int(next_uid),
                    "target_state":int(child.state),
                })

        if step%settings.sampled_audit_interval==0:
            attempt=0
            while True:
                audit=sampled_molecular_edge_audit_v20(
                    basis,provider,dt,
                    update,graph.settings,
                    step=step,
                    priority_count=
                        settings.sampled_audit_priority_pairs,
                    random_count=
                        settings.sampled_audit_random_pairs,
                    audit_search_factor=
                        settings.sampled_audit_search_factor,
                    seed=settings.sampled_audit_seed,
                    violation_factor=
                        settings.sampled_audit_violation_factor,
                )
                row={
                    **audit.as_dict(),
                    "attempt":int(attempt),
                    "search_overlap_floor":
                        float(graph.settings.search_overlap_floor),
                }
                sampled_audits.append(row)
                if audit.passed:
                    break
                if attempt>=settings.max_sampled_audit_relaxations:
                    events.append({
                        "kind":"sampled_molecular_audit_unresolved",
                        **row,
                    })
                    break

                previous=float(
                    graph.settings.search_overlap_floor
                )
                graph.relax_search_floor(
                    settings.sampled_audit_relaxation_factor
                )
                update=graph.update(basis)
                mats=build_sparse_molecular_matrices_v20(
                    basis,update
                )
                events.append({
                    "kind":"sampled_molecular_search_relaxation",
                    "step":int(step),
                    "time":float(step*dt),
                    "attempt":int(attempt),
                    "previous_search_overlap_floor":
                        previous,
                    "new_search_overlap_floor":
                        float(graph.settings.search_overlap_floor),
                    "maximum_sampled_score":
                        float(audit.maximum_score),
                    "violation_count":
                        int(audit.violation_count),
                })
                attempt+=1

        if step%int(store_every)==0:
            record(step,mats,update)

    final_settings=replace(
        settings,
        graph=graph.settings,
    )
    sentinel_provider_final=(
        sentinel_provider_factory()
        if sentinel_provider_factory is not None
        else provider
    )
    final_sentinel=_dense_sentinel(
        basis,
        sentinel_provider_final,
        dt,mats,final_settings
    )

    sentinel_diagnostics={
        "initial":(
            sentinel_provider_initial.diagnostics_dict()
            if (
                sentinel_provider_factory is not None
                and hasattr(
                    sentinel_provider_initial,
                    "diagnostics_dict",
                )
            )
            else None
        ),
        "final":(
            sentinel_provider_final.diagnostics_dict()
            if (
                sentinel_provider_factory is not None
                and hasattr(
                    sentinel_provider_final,
                    "diagnostics_dict",
                )
            )
            else None
        ),
    }

    return {
        "records":records,
        "sampled_audits":sampled_audits,
        "sentinels":{
            "initial":initial_sentinel,
            "final":final_sentinel,
        },
        "sentinel_provider_diagnostics":
            sentinel_diagnostics,
        "events":events,
        "final_basis":basis,
        "final_coefficients":C,
        "final_matrices":mats,
        "final_update":update,
        "graph_total_exact_pair_checks":
            int(graph.total_exact_pair_checks),
        "settings":{
            "dt":float(dt),
            "steps":int(steps),
            "control":asdict(settings),
        },
        "provider_diagnostics":
            provider.diagnostics_dict()
            if hasattr(provider,"diagnostics_dict")
            else None,
    }


def run_dense_molecular_reference_dynamics_v20(
    initial_basis,
    C0,
    provider,
    *,
    dt=0.01,
    steps=50,
    store_every=5,
):
    """Dense molecular reference using the same pair-centroid approximation."""
    basis=_copy_basis(initial_basis)
    C=np.asarray(C0,dtype=complex).copy()
    graph_settings=MolecularSparseSettingsV20(
        enter_score=1e-30,
        exit_score=1e-30,
        search_overlap_floor=1e-14,
        local_omitted_score_l2_budget=0.0,
        use_kdtree=False,
    )

    dense=build_dense_molecular_reference_v20(
        basis,provider,dt,graph_settings
    )
    S=sparse.csr_matrix(dense["S"])
    H=sparse.csr_matrix(dense["H"])
    T=sparse.csr_matrix(dense["T_seed"])
    norm=float(np.real(np.vdot(C,S@C)))
    C=C/np.sqrt(norm)

    records=[]
    def record(step):
        records.append({
            "step":int(step),
            "time":float(step*dt),
            "norm":float(
                np.real(np.vdot(C,S@C))
            ),
            "condition_number":
                float(np.linalg.cond(S.toarray())),
        })
    record(0)

    for step in range(1,int(steps)+1):
        new_basis=_nuclear_velocity_verlet_step(
            basis,provider,float(dt)
        )
        dense_new=build_dense_molecular_reference_v20(
            new_basis,provider,dt,graph_settings
        )
        Snew=sparse.csr_matrix(
            dense_new["S"]
        )
        Hnew=sparse.csr_matrix(
            dense_new["H"]
        )
        Tnew=sparse.csr_matrix(
            dense_new["T_seed"]
        )
        seed=0.5*(T+Tnew)
        Tmid=sparse_metric_compatible_connection(
            S,Snew,dt,seed
        )
        C=sparse_moving_basis_midpoint_cayley_step(
            C,S,H,Snew,Hnew,Tmid,dt
        )
        basis=new_basis
        S,H,T=Snew,Hnew,Tnew
        if step%int(store_every)==0:
            record(step)

    return {
        "records":records,
        "final_basis":basis,
        "final_coefficients":C,
        "final_S":S,
        "provider_diagnostics":
            provider.diagnostics_dict()
            if hasattr(provider,"diagnostics_dict")
            else None,
    }
