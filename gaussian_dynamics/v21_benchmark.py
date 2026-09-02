from dataclasses import dataclass, asdict
import importlib.util, time, numpy as np
from .analytic_molecular_backend_v19 import AnalyticMolecularLVCBackendV19, default_diatomic_two_mode_map_v19
from .indexed_molecular_provider_v20 import IndexedTrackedMolecularDirectProviderV20
from .electronic_operator_v21 import ElectronicOperatorProviderAdapterV21
from .complex_gauge_v21 import PhaseMixingGaugeV21, GaugeTransformedOperatorProviderV21, random_unitary_v21
from .subspace_tracking_v21 import procrustes_subspace_alignment_v21
from .wilson_loop_v21 import gauge_transform_cycle_links_v21, sorted_wilson_eigenphases_v21
from .synthetic_operator_provider_v21 import SyntheticLinearOperatorConfigV21, SyntheticLinearOperatorProviderV21
from .block_sparse_molecular_v21 import BlockMolecularTBFV21, BlockSparseSettingsV21, BlockSparseMolecularGraphV21, build_block_sparse_matrices_v21, build_dense_block_reference_v21
from .block_dynamics_v21 import PrescribedBlockDynamicsSettingsV21, prescribed_linear_basis_v21, run_prescribed_block_dynamics_v21, gauge_block_matrices_v21, gauge_covariance_errors_v21, gauge_mapped_coefficient_error_v21

@dataclass(frozen=True)
class V21AcceptanceThresholds:
    max_point_covariance_error:float=1e-11; max_force_covariance_error:float=1e-10
    max_block_covariance_error:float=1e-11; max_block_score_error:float=1e-11
    max_finest_gauge_propagation_error:float=2e-9; min_gauge_propagation_order:float=1.8; max_gauge_norm_drift:float=1e-7
    max_subspace_residual:float=1e-11; min_subspace_singular_value:float=0.999999; max_wilson_phase_error:float=1e-10
    max_zero_budget_error:float=1e-5; min_topology_entered_edges:int=10; min_topology_exited_edges:int=5; required_max_nstate:int=8; max_nstate_nnz_scaling_error:float=1e-12

def _molecular_provider():
    g=default_diatomic_two_mode_map_v19()
    return ElectronicOperatorProviderAdapterV21(IndexedTrackedMolecularDirectProviderV20(AnalyticMolecularLVCBackendV19(g),g,rebuild_batch=8))

def _gauge(s=2,seed=219):
    if s==2:
        B=np.array([[.35,-.18],[-.21,.27]]); off=np.array([.1,-.3])
    else:
        B=np.column_stack([np.linspace(.05,.15,s),np.linspace(-.07,.09,s)]); off=np.linspace(-.2,.3,s)
    return PhaseMixingGaugeV21(random_unitary_v21(s,seed),B,off)

def _all_edges(): return BlockSparseSettingsV21(1e-14,1e-14,1e-14,local_omitted_score_l2_budget=0,use_kdtree=False)

def _point_covariance():
    base=_molecular_provider(); gauge=_gauge(); gp=GaugeTransformedOperatorProviderV21(_molecular_provider(),gauge); q=np.array([-.31,.42]); a=base.evaluate(q); b=gp.evaluate(q); G=gauge.matrix(q)
    rel=lambda A,B:float(np.linalg.norm(A-B,'fro')/max(np.linalg.norm(B,'fro'),1e-30))
    c=np.array([.8+.1j,-.2+.5j]); c/=np.linalg.norm(c); cg=G.conj().T@c
    return {'H_relative_error':rel(b.H,G.conj().T@a.H@G),'maximum_dH_relative_error':max(rel(b.dH_dq[k],G.conj().T@a.dH_dq[k]@G) for k in range(a.nq)),'force_error':float(np.linalg.norm(a.force_expectation(c)-b.force_expectation(cg)))}

def _block_covariance():
    A=1.3*np.eye(2); basis=[BlockMolecularTBFV21(0,np.array([-.7,.35]),np.array([1.2,.1]),A),BlockMolecularTBFV21(1,np.array([0,.42]),np.array([-.4,.2]),A),BlockMolecularTBFV21(2,np.array([.75,.33]),np.array([.7,-.1]),A)]
    qd=np.array([[.12,.03],[-.05,.08],[.09,-.04]]); pd=np.array([[.01,0],[-.02,.01],[0,-.01]]); gauge=_gauge()
    base=build_dense_block_reference_v21(basis,_molecular_provider(),.01,qd,pd,_all_edges()); gp=GaugeTransformedOperatorProviderV21(_molecular_provider(),gauge); trans=build_dense_block_reference_v21(basis,gp,.01,qd,pd,_all_edges()); G,dG=gauge_block_matrices_v21(gauge,basis,qd); out=gauge_covariance_errors_v21(base,trans,G,dG)
    out['maximum_edge_score_error']=max(abs(base['pair_data'][e].score-trans['pair_data'][e].score) for e in base['pair_data']); return out

def _gauge_propagation():
    A=1.2*np.eye(2); basis=[BlockMolecularTBFV21(0,np.array([-.75,.31]),np.array([1,.1]),A),BlockMolecularTBFV21(1,np.array([-.1,.43]),np.array([-.2,.15]),A),BlockMolecularTBFV21(2,np.array([.62,.36]),np.array([.5,-.12]),A)]; vel=np.array([[.12,.025],[-.04,.05],[.08,-.03]]); gauge=_gauge(); settings=PrescribedBlockDynamicsSettingsV21(_all_edges(),True); rng=np.random.default_rng(123); C0=rng.normal(size=6)+1j*rng.normal(size=6); rows=[]
    for dt in (.02,.01,.005):
        steps=round(.10/dt); G0,_=gauge_block_matrices_v21(gauge,basis,vel); Cg0=G0.conj().T@C0
        b=run_prescribed_block_dynamics_v21(basis,C0,_molecular_provider(),vel,dt=dt,steps=steps,settings=settings,store_every=steps)
        t=run_prescribed_block_dynamics_v21(basis,Cg0,GaugeTransformedOperatorProviderV21(_molecular_provider(),gauge),vel,dt=dt,steps=steps,settings=settings,store_every=steps)
        Gf,_=gauge_block_matrices_v21(gauge,b['final_basis'],vel); err=gauge_mapped_coefficient_error_v21(b['final_coefficients'],t['final_coefficients'],b['final_S'],Gf); drift=max([abs(x['norm']-1) for x in b['records']+t['records']]); rows.append({'dt':dt,'steps':steps,'gauge_mapped_coefficient_error':err,'maximum_norm_drift':drift})
    orders=[float(np.log(rows[i]['gauge_mapped_coefficient_error']/rows[i+1]['gauge_mapped_coefficient_error'])/np.log(2)) for i in range(2)]
    return {'rows':rows,'observed_orders':orders,'minimum_observed_order':min(orders)}

def _geometry():
    sub=procrustes_subspace_alignment_v21(random_unitary_v21(8,2111)); links=[random_unitary_v21(4,2200+k) for k in range(6)]; gauges=[random_unitary_v21(4,2300+k) for k in range(6)]; a=sorted_wilson_eigenphases_v21(links); b=sorted_wilson_eigenphases_v21(gauge_transform_cycle_links_v21(links,gauges)); return {'subspace':sub.as_dict(),'wilson_max_phase_error':float(np.max(abs(a-b)))}

def _synthetic(s,seed=404): return GaugeTransformedOperatorProviderV21(SyntheticLinearOperatorProviderV21(SyntheticLinearOperatorConfigV21(nstate=s,nq=2,seed=seed)),_gauge(s,seed+1))

def _block_convergence():
    N=16; spacing=.85; A=2.4*np.eye(2); x0=-.5*spacing*(N-1); basis=[BlockMolecularTBFV21(i,np.array([x0+i*spacing,.22+.03*np.sin(.7*i)]),np.array([.1*(-1)**i,.01*np.cos(i)]),A) for i in range(N)]; qd=np.array([[.04*(-1)**i,.01*np.sin(i)] for i in range(N)]); pd=np.zeros_like(qd); dense=build_dense_block_reference_v21(basis,_synthetic(4),.01,qd,pd,_all_edges()); rel=lambda A,B:float(np.linalg.norm(A-B,'fro')/max(np.linalg.norm(B,'fro'),1e-30))
    def row(th,budget):
        st=BlockSparseSettingsV21(th,.7*th,1e-6,local_omitted_score_l2_budget=budget); u=BlockSparseMolecularGraphV21(_synthetic(4),.01,st).update(basis,qd,pd); m=build_block_sparse_matrices_v21(basis,u); return {'active_edges':len(u.active_edges),'S_error':rel(m.S.toarray(),dense['S']),'H_error':rel(m.H.toarray(),dense['H']),'T_error':rel(m.T_seed.toarray(),dense['T_seed']),'promoted_edges':u.budget_promoted_edges,'remaining_score_l2':u.omitted_score_l2}
    trs=[]
    for x in (.5,.35,.25,.18,.12,.08,.04): r=row(x,1e9); r['enter_score']=x; trs.append(r)
    brs=[]
    for b in (1e9,.2,.1,.05,.02,.01,0): r=row(.25,b); r['budget']=b; brs.append(r)
    mono=lambda rows,key:all(rows[i+1][key]<=rows[i][key]+1e-14 for i in range(len(rows)-1))
    return {'threshold_rows':trs,'budget_rows':brs,'threshold_monotone':{k:mono(trs,k) for k in ('S_error','H_error','T_error')},'budget_monotone':{k:mono(brs,k) for k in ('S_error','H_error','T_error')}}

def _topology():
    p=_synthetic(4,435); A=2*np.eye(2); qs=[[-2.8,.3],[-2.1,.34],[-1.4,.38],[1.4,.36],[2.1,.32],[2.8,.28]]; vs=np.array([[.9,0],[.7,0],[.5,0],[-.5,0],[-.7,0],[-.9,0]],float); basis=[BlockMolecularTBFV21(i,np.array(q,float),np.array([.1,0]),A) for i,q in enumerate(qs)]; g=BlockSparseMolecularGraphV21(p,.02,BlockSparseSettingsV21(.5,.35,1e-4,local_omitted_score_l2_budget=1e9)); events=[]; ent=ext=maxa=0
    for step in range(51):
        u=g.update(prescribed_linear_basis_v21(basis,vs,.1*step),vs,np.zeros_like(vs)); ent+=u.entered_edges; ext+=u.exited_edges; maxa=max(maxa,len(u.active_edges));
        if u.entered_edges or u.exited_edges: events.append({'step':step,'time':.1*step,'active_edges':len(u.active_edges),'entered_edges':u.entered_edges,'exited_edges':u.exited_edges})
    return {'events':events,'total_entered_edges':ent,'total_exited_edges':ext,'maximum_active_edges':maxa,'final_active_edges':events[-1]['active_edges'],'total_exact_pair_checks':g.total_exact_pair_checks}

def _nstate():
    N=24; spacing=1.2; A=2.2*np.eye(2); x0=-.5*spacing*(N-1); basis=[BlockMolecularTBFV21(i,np.array([x0+i*spacing,.25+.02*np.sin(i)]),np.array([.08*(-1)**i,0]),A) for i in range(N)]; qd=np.array([[.03*(-1)**i,.005*np.cos(i)] for i in range(N)]); rows=[]
    for s in (2,4,8):
        st=BlockSparseSettingsV21(.25,.15,1e-4,local_omitted_score_l2_budget=.02); u=BlockSparseMolecularGraphV21(_synthetic(s,500+s),.01,st).update(basis,qd,np.zeros_like(qd)); m=build_block_sparse_matrices_v21(basis,u); rows.append({'nstate':s,'dimension':m.dimension,'active_edges':len(u.active_edges),'H_nnz':m.H.nnz,'S_nnz':m.S.nnz,'T_nnz':m.T_seed.nnz,'H_density':m.H.nnz/(m.dimension**2)})
    norm=[r['H_nnz']/r['nstate']**2 for r in rows]; err=(max(norm)-min(norm))/max(max(norm),1); return {'rows':rows,'H_nnz_over_nstate_squared':norm,'relative_nnz_scaling_error':err}

def run_v021_release_benchmark():
    t=time.perf_counter(); point=_point_covariance(); block=_block_covariance(); prop=_gauge_propagation(); geom=_geometry(); conv=_block_convergence(); topo=_topology(); ns=_nstate(); th=V21AcceptanceThresholds(); zb=conv['budget_rows'][-1]
    checks={'point':max(point['H_relative_error'],point['maximum_dH_relative_error'])<=th.max_point_covariance_error and point['force_error']<=th.max_force_covariance_error,'block':max(block['S_relative_error'],block['H_relative_error'],block['T_relative_error'])<=th.max_block_covariance_error and block['maximum_edge_score_error']<=th.max_block_score_error,'propagation':prop['rows'][-1]['gauge_mapped_coefficient_error']<=th.max_finest_gauge_propagation_error and prop['minimum_observed_order']>=th.min_gauge_propagation_order and max(r['maximum_norm_drift'] for r in prop['rows'])<=th.max_gauge_norm_drift,'geometry':geom['subspace']['antihermitian_residual']<=th.max_subspace_residual and geom['subspace']['minimum_singular_value']>=th.min_subspace_singular_value and geom['wilson_max_phase_error']<=th.max_wilson_phase_error,'sparse_convergence':all(conv['threshold_monotone'].values()) and all(conv['budget_monotone'].values()) and max(zb['S_error'],zb['H_error'],zb['T_error'])<=th.max_zero_budget_error,'topology':topo['total_entered_edges']>=th.min_topology_entered_edges and topo['total_exited_edges']>=th.min_topology_exited_edges,'nstate':max(r['nstate'] for r in ns['rows'])>=th.required_max_nstate and ns['relative_nnz_scaling_error']<=th.max_nstate_nnz_scaling_error}
    return {'release':'v0.21','theme':'complex/block/gauge-ready generalized framework; no spin physics introduced','point_covariance':point,'block_covariance':block,'gauge_propagation':prop,'subspace_and_wilson':geom,'block_sparse_convergence':conv,'dynamic_topology':topo,'nstate_scaling':ns,'pyscf':{'installed_in_build_environment':bool(importlib.util.find_spec('pyscf')),'runtime_validated':False},'benchmark_wall_seconds':time.perf_counter()-t,'acceptance':{'passed':all(checks.values()),'checks':checks,'thresholds':asdict(th)}}
