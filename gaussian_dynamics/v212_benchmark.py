from dataclasses import dataclass,asdict
from pathlib import Path
import importlib.util
import numpy as np

from .synthetic_operator_provider_v21 import SyntheticLinearOperatorProviderV21,SyntheticLinearOperatorConfigV21
from .complex_gauge_v21 import PhaseMixingGaugeV21,GaugeTransformedOperatorProviderV21,random_unitary_v21
from .block_sparse_molecular_v21 import (
    BlockMolecularTBFV21,BlockSparseSettingsV21,BlockSparseMolecularGraphV21,
    build_block_sparse_matrices_v21,build_dense_block_reference_v21,
)
from .block_dynamics_v21 import gauge_block_matrices_v21,gauge_covariance_errors_v21,gauge_mapped_coefficient_error_v21
from .self_consistent_block_v212 import SelfConsistentBlockSettingsV212,run_self_consistent_block_dynamics_v212
from .electronic_observables_v212 import ElectronicObservableV212,build_electronic_observable_matrix_v212,observable_expectation_v212
from .subspace_provider_v212 import SubspaceAwareOperatorProviderV212,SubspaceTrackingSettingsV212
from .block_basis_lifecycle_v212 import insert_zero_block_v212,prune_block_projected_v212
from .complex_dtype_audit_v212 import audit_pre_soc_complex_core_v212
from .v21_benchmark import run_v021_release_benchmark


@dataclass(frozen=True)
class V212AcceptanceThresholds:
    max_unequal_width_S_covariance:float=2e-11
    max_unequal_width_H_covariance:float=2e-11
    max_unequal_width_T_covariance:float=2e-11
    max_all_edge_sparse_dense_error:float=2e-12
    max_self_consistent_position_error:float=2e-11
    max_self_consistent_momentum_error:float=2e-11
    max_self_consistent_gauge_error:float=1e-11
    min_self_consistent_gauge_order:float=1.8
    max_self_consistent_norm_drift:float=2e-9
    max_observable_gauge_error:float=2e-10
    min_subspace_singular_value:float=0.999999
    max_zero_birth_change:float=1e-14
    max_zero_birth_prune_loss:float=1e-18
    require_dtype_audit:bool=True
    require_imaginary_signal:float=1e-6
    require_inherited_v021:bool=True


def _provider(nstate=3,seed=21201,mass=25.0):
    return SyntheticLinearOperatorProviderV21(
        SyntheticLinearOperatorConfigV21(nstate=nstate,nq=2,seed=seed,mass=mass,base_scale=.03,derivative_scale=.012)
    )


def _gauge(nstate,seed=21202):
    return PhaseMixingGaugeV21(
        random_unitary_v21(nstate,seed),
        np.column_stack([np.linspace(.04,.12,nstate),np.linspace(-.08,.07,nstate)]),
        np.linspace(-.25,.3,nstate),
    )


def _unequal_width_campaign():
    basis=[
        BlockMolecularTBFV21(0,np.array([-.6,.2]),np.array([.2,.05]),np.array([[1.2,.1],[.1,1.7]])),
        BlockMolecularTBFV21(1,np.array([.1,.35]),np.array([-.1,.08]),np.array([[1.8,-.12],[-.12,1.1]])),
        BlockMolecularTBFV21(2,np.array([.75,.15]),np.array([.15,-.04]),np.diag([.9,2.0])),
    ]
    qdot=np.array([[.08,.02],[-.03,.05],[.06,-.02]])
    pdot=np.zeros_like(qdot)
    gauge=_gauge(3)
    settings=BlockSparseSettingsV21(enter_score=1e-14,exit_score=1e-14,search_overlap_floor=1e-14,local_omitted_score_l2_budget=0,use_kdtree=False)
    pb=_provider(); pg=GaugeTransformedOperatorProviderV21(_provider(),gauge)
    A=build_dense_block_reference_v21(basis,pb,.01,qdot,pdot,settings)
    B=build_dense_block_reference_v21(basis,pg,.01,qdot,pdot,settings)
    G,dG=gauge_block_matrices_v21(gauge,basis,qdot)
    cov=gauge_covariance_errors_v21(A,B,G,dG)

    ps=_provider()
    graph=BlockSparseMolecularGraphV21(ps,.01,settings)
    update=graph.update(basis,qdot,pdot)
    sparse=build_block_sparse_matrices_v21(basis,update)
    dense=build_dense_block_reference_v21(basis,ps,.01,qdot,pdot,settings)
    def rel(X,Y): return float(np.linalg.norm(X-Y,'fro')/max(np.linalg.norm(Y,'fro'),1e-30))
    exact={
        'S':rel(sparse.S.toarray(),dense['S']),
        'H':rel(sparse.H.toarray(),dense['H']),
        'T':rel(sparse.T_seed.toarray(),dense['T_seed']),
    }
    return {'covariance':cov,'all_edge_sparse_dense':exact,'width_eigenvalues':[np.linalg.eigvalsh(b.A).tolist() for b in basis]}


def _self_consistent_campaign():
    def base(): return SyntheticLinearOperatorProviderV21(SyntheticLinearOperatorConfigV21(nstate=2,nq=2,seed=21212,mass=30.0,base_scale=.025,derivative_scale=.01))
    gauge=PhaseMixingGaugeV21(random_unitary_v21(2,21213),np.array([[.14,-.06],[-.08,.11]]),np.array([.2,-.3]))
    def basis(): return [
        BlockMolecularTBFV21(0,np.array([-.7,.2]),np.array([.3,.05]),np.diag([1.2,1.5])),
        BlockMolecularTBFV21(1,np.array([.05,.3]),np.array([-.15,.08]),np.array([[1.7,.15],[.15,1.3]])),
        BlockMolecularTBFV21(2,np.array([.8,.15]),np.array([.2,-.04]),np.diag([1.0,1.8])),
    ]
    C0=np.array([.7+.1j,.25-.15j,.4+.2j,-.1+.3j,.3-.2j,.2+.1j])
    settings=SelfConsistentBlockSettingsV212(
        graph=BlockSparseSettingsV21(enter_score=1e-14,exit_score=1e-14,search_overlap_floor=1e-14,local_omitted_score_l2_budget=0,use_kdtree=False),
        use_dense_reference=True,corrector_iterations=3,momentum_tolerance=1e-12,
    )
    rows=[]
    for dt in (.01,.005,.0025):
        steps=int(round(.05/dt))
        ob=run_self_consistent_block_dynamics_v212(basis(),C0,base(),dt=dt,steps=steps,settings=settings,store_every=steps)
        gp=GaugeTransformedOperatorProviderV21(base(),gauge)
        G0,_=gauge_block_matrices_v21(gauge,basis(),np.zeros((3,2)))
        og=run_self_consistent_block_dynamics_v212(basis(),G0.conj().T@C0,gp,dt=dt,steps=steps,settings=settings,store_every=steps)
        qerr=max(np.linalg.norm(a.q-b.q) for a,b in zip(ob['final_basis'],og['final_basis']))
        perr=max(np.linalg.norm(a.p-b.p) for a,b in zip(ob['final_basis'],og['final_basis']))
        Gf,_=gauge_block_matrices_v21(gauge,ob['final_basis'],np.zeros((3,2)))
        cerr=gauge_mapped_coefficient_error_v21(ob['final_coefficients'],og['final_coefficients'],ob['final_S'],Gf)
        rows.append({'dt':dt,'steps':steps,'position_error':float(qerr),'momentum_error':float(perr),'coefficient_error':float(cerr),'base_norm_drift':ob['maximum_norm_drift'],'gauge_norm_drift':og['maximum_norm_drift']})
    orders=[float(np.log(rows[i]['coefficient_error']/rows[i+1]['coefficient_error'])/np.log(2)) for i in range(len(rows)-1)]
    return {'rows':rows,'observed_orders':orders,'minimum_order':min(orders)}


def _observable_campaign():
    gauge=_gauge(3,21221)
    def base(): return _provider(3,21220,20.0)
    basis=[
        BlockMolecularTBFV21(0,np.array([-.45,.2]),np.array([.1,0]),np.eye(2)),
        BlockMolecularTBFV21(1,np.array([.55,.25]),np.array([-.08,.03]),1.3*np.eye(2)),
    ]
    qdot=np.zeros((2,2)); pdot=np.zeros_like(qdot)
    settings=BlockSparseSettingsV21(enter_score=1e-14,exit_score=1e-14,search_overlap_floor=1e-14,local_omitted_score_l2_budget=0,use_kdtree=False)
    pb=base(); pg=GaugeTransformedOperatorProviderV21(base(),gauge)
    db=build_dense_block_reference_v21(basis,pb,.01,qdot,pdot,settings)
    dg=build_dense_block_reference_v21(basis,pg,.01,qdot,pdot,settings)
    obs=ElectronicObservableV212('physical dH/dq0',lambda snap:snap.point.dH_dq[0])
    Ob=build_electronic_observable_matrix_v212(basis,pb,obs)
    Og=build_electronic_observable_matrix_v212(basis,pg,obs)
    C=np.array([.6+.1j,.2-.2j,-.1+.3j,.3+.2j,.15-.1j,.25+.05j])
    G,_=gauge_block_matrices_v21(gauge,basis,qdot)
    eb=observable_expectation_v212(C,db['S'],Ob); eg=observable_expectation_v212(G.conj().T@C,dg['S'],Og)
    p=pg.evaluate(basis[0].q)
    return {'base_expectation':eb,'gauge_expectation':eg,'absolute_error':abs(eb-eg),'imaginary_H_signal':float(np.max(np.abs(np.imag(p.H)))),'imaginary_dH_signal':float(np.max(np.abs(np.imag(p.dH_dq)))),'imaginary_connection_signal':float(np.max(np.abs(np.imag(p.connection_q))))}


def _subspace_campaign():
    base=_provider(5,21230)
    gauge=_gauge(5,21231)
    raw=GaugeTransformedOperatorProviderV21(base,gauge)
    p=SubspaceAwareOperatorProviderV212(raw,dimension=2,settings=SubspaceTrackingSettingsV212(minimum_singular_value=.999999,ambiguity_policy='raise',rebuild_batch=3))
    points=[np.array([x,.2+.01*np.sin(3*x)]) for x in np.linspace(-.5,.5,9)]
    for q in points: p.evaluate_snapshot(q)
    return p.diagnostics_dict()['subspace']


def _lifecycle_campaign():
    s=2; A=np.eye(2)
    basis=[BlockMolecularTBFV21(0,np.array([-.5,0]),np.zeros(2),A),BlockMolecularTBFV21(1,np.array([.5,0]),np.zeros(2),A)]
    C=np.array([.8+.1j,.2-.3j,-.1+.2j,.4+.05j])
    new=BlockMolecularTBFV21(2,np.array([1.2,0]),np.zeros(2),1.4*A)
    born,Cb=insert_zero_block_v212(basis,C,new,s)
    birth_change=float(np.linalg.norm(Cb[:4]-C)+np.linalg.norm(Cb[4:]))
    X=np.eye(6,dtype=complex); X[0,4]=.04; X[4,2]=.12; X[1,5]=.03; X[5,3]=.10
    S=X.conj().T@X
    pr=prune_block_projected_v212(born,Cb,S,s,2)
    return {'birth_state_change':birth_change,'prune_projection_loss':pr.projection_loss,'pruned_coefficient_error':float(np.linalg.norm(pr.coefficients-C)),'retained_condition_number':pr.retained_condition_number}


def run_v0212_release_benchmark():
    unequal=_unequal_width_campaign()
    selfc=_self_consistent_campaign()
    obs=_observable_campaign()
    sub=_subspace_campaign()
    lifecycle=_lifecycle_campaign()
    dtype=audit_pre_soc_complex_core_v212(Path(__file__).resolve().parent).as_dict()
    inherited=run_v021_release_benchmark()['acceptance']['passed']
    thresholds=V212AcceptanceThresholds()
    checks={
        'unequal_width_S_covariance':unequal['covariance']['S_relative_error']<=thresholds.max_unequal_width_S_covariance,
        'unequal_width_H_covariance':unequal['covariance']['H_relative_error']<=thresholds.max_unequal_width_H_covariance,
        'unequal_width_T_covariance':unequal['covariance']['T_relative_error']<=thresholds.max_unequal_width_T_covariance,
        'all_edge_sparse_dense':max(unequal['all_edge_sparse_dense'].values())<=thresholds.max_all_edge_sparse_dense_error,
        'self_consistent_positions':max(r['position_error'] for r in selfc['rows'])<=thresholds.max_self_consistent_position_error,
        'self_consistent_momenta':max(r['momentum_error'] for r in selfc['rows'])<=thresholds.max_self_consistent_momentum_error,
        'self_consistent_coefficients':selfc['rows'][-1]['coefficient_error']<=thresholds.max_self_consistent_gauge_error,
        'self_consistent_order':selfc['minimum_order']>=thresholds.min_self_consistent_gauge_order,
        'self_consistent_norm':max(max(r['base_norm_drift'],r['gauge_norm_drift']) for r in selfc['rows'])<=thresholds.max_self_consistent_norm_drift,
        'observable_gauge_invariance':obs['absolute_error']<=thresholds.max_observable_gauge_error,
        'complex_signal_preserved':min(obs['imaginary_H_signal'],obs['imaginary_dH_signal'],obs['imaginary_connection_signal'])>=thresholds.require_imaginary_signal,
        'subspace_continuity':sub['minimum_seen_singular_value']>=thresholds.min_subspace_singular_value and sub['subspace_ambiguities']==0,
        'zero_block_birth':lifecycle['birth_state_change']<=thresholds.max_zero_birth_change,
        'zero_block_prune':lifecycle['prune_projection_loss']<=thresholds.max_zero_birth_prune_loss,
        'dtype_audit':dtype['passed'] if thresholds.require_dtype_audit else True,
        'inherited_v021':bool(inherited) if thresholds.require_inherited_v021 else True,
    }
    return {
        'release':'v0.21.2',
        'theme':'pre-SOC integration hardening; generalized spin-neutral framework',
        'unequal_width_block':unequal,
        'self_consistent_block_dynamics':selfc,
        'electronic_observables':obs,
        'subspace_provider':sub,
        'adaptive_block_lifecycle':lifecycle,
        'complex_dtype_audit':dtype,
        'inherited_v021_acceptance':bool(inherited),
        'pyscf':{
            'installed_in_build_environment':bool(importlib.util.find_spec('pyscf') is not None),
            'runtime_validated':False,
            'note':'Real PySCF runtime validation remains a separate empirical milestone; it is not required for the analytic first SOC model.'
        },
        'acceptance':{'passed':bool(all(checks.values())),'checks':checks,'thresholds':asdict(thresholds)},
    }
