from dataclasses import dataclass, asdict
import numpy as np
from scipy import sparse
from .block_sparse_molecular_v21 import BlockMolecularTBFV21, BlockSparseSettingsV21, BlockSparseMolecularGraphV21, build_block_sparse_matrices_v21, build_dense_block_reference_v21, block_diagonal_gauge_v21
from .sparse_pair_matrices_v16 import sparse_metric_compatible_connection, sparse_moving_basis_midpoint_cayley_step, sparse_generalized_norm

@dataclass(frozen=True)
class PrescribedBlockDynamicsSettingsV21:
    graph: BlockSparseSettingsV21=BlockSparseSettingsV21()
    use_dense_reference: bool=False
    def validate(self): self.graph.validate(); return self

def prescribed_linear_basis_v21(initial_basis, velocities, time):
    initial_basis=list(initial_basis); velocities=np.asarray(velocities,float); t=float(time)
    return [BlockMolecularTBFV21(b.uid,b.q+t*velocities[i],b.p.copy(),b.A.copy()) for i,b in enumerate(initial_basis)]

def _dense_as_sparse_v21(basis,provider,dt,qdots,pdots,graph_settings):
    d=build_dense_block_reference_v21(basis,provider,dt,qdots,pdots,graph_settings)
    return {'S':sparse.csr_matrix(d['S']),'H':sparse.csr_matrix(d['H']),'T_seed':sparse.csr_matrix(d['T_seed']),'active_edges':tuple((i,j) for i in range(len(basis)) for j in range(i+1,len(basis))),'update':None,'dense':d}

def run_prescribed_block_dynamics_v21(initial_basis,C0,provider,velocities,*,dt=0.005,steps=20,settings=PrescribedBlockDynamicsSettingsV21(),store_every=5):
    settings=settings.validate(); initial_basis=list(initial_basis); velocities=np.asarray(velocities,float); pdots=np.zeros_like(velocities); dt=float(dt); steps=int(steps)
    graph=None if settings.use_dense_reference else BlockSparseMolecularGraphV21(provider,dt,settings.graph)
    def matrices_at(step):
        basis=prescribed_linear_basis_v21(initial_basis,velocities,step*dt)
        if settings.use_dense_reference: return basis,_dense_as_sparse_v21(basis,provider,dt,velocities,pdots,settings.graph)
        u=graph.update(basis,velocities,pdots); m=build_block_sparse_matrices_v21(basis,u)
        return basis,{'S':m.S,'H':m.H,'T_seed':m.T_seed,'active_edges':m.active_edges,'update':u,'object':m}
    basis,current=matrices_at(0); C=np.asarray(C0,complex).copy(); C=C/np.sqrt(sparse_generalized_norm(C,current['S']))
    records=[]; events=[]
    def record(step):
        u=current['update']; records.append({'step':int(step),'time':float(step*dt),'norm':sparse_generalized_norm(C,current['S']),'condition_number':float(np.linalg.cond(current['S'].toarray())),'active_edges':len(current['active_edges']),'entered_edges':0 if u is None else int(u.entered_edges),'exited_edges':0 if u is None else int(u.exited_edges),'sparsity_fraction':0.0 if u is None else float(u.sparsity_fraction)})
    record(0)
    for step in range(1,steps+1):
        new_basis,new=matrices_at(step); seed=.5*(current['T_seed']+new['T_seed']); Tmid=sparse_metric_compatible_connection(current['S'],new['S'],dt,seed)
        C=sparse_moving_basis_midpoint_cayley_step(C,current['S'],current['H'],new['S'],new['H'],Tmid,dt)
        if new['update'] is not None and (new['update'].entered_edges or new['update'].exited_edges): events.append({'step':step,'time':step*dt,'entered_edges':new['update'].entered_edges,'exited_edges':new['update'].exited_edges,'active_edges':len(new['active_edges'])})
        basis,current=new_basis,new
        if step%int(store_every)==0: record(step)
    return {'records':records,'topology_events':events,'final_basis':basis,'final_coefficients':C,'final_S':current['S'],'final_H':current['H'],'final_T_seed':current['T_seed'],'final_active_edges':current['active_edges'],'graph_total_exact_pair_checks':None if graph is None else graph.total_exact_pair_checks,'settings':{'dt':dt,'steps':steps,'control':asdict(settings)}}

def gauge_block_matrices_v21(gauge,basis,velocities):
    G=[]; dG=[]
    for b,v in zip(basis,np.asarray(velocities,float)):
        G.append(gauge.matrix(b.q)); dG.append(gauge.velocity_derivative(b.q,v))
    return block_diagonal_gauge_v21(G), sparse.block_diag([sparse.csr_matrix(x) for x in dG],format='csr')

def gauge_covariance_errors_v21(base_dense,gauge_dense,gauge_block,gauge_dot_block):
    G=gauge_block.toarray(); dG=gauge_dot_block.toarray(); Sb=np.asarray(base_dense['S'],complex); Hb=np.asarray(base_dense['H'],complex); Tb=np.asarray(base_dense['T_seed'],complex)
    expected=(G.conj().T@Sb@G, G.conj().T@Hb@G, G.conj().T@Tb@G+G.conj().T@Sb@dG)
    actual=(np.asarray(gauge_dense['S'],complex),np.asarray(gauge_dense['H'],complex),np.asarray(gauge_dense['T_seed'],complex))
    def rel(A,B): return float(np.linalg.norm(A-B,'fro')/max(np.linalg.norm(B,'fro'),1e-30))
    return {'S_relative_error':rel(actual[0],expected[0]),'H_relative_error':rel(actual[1],expected[1]),'T_relative_error':rel(actual[2],expected[2])}

def gauge_mapped_coefficient_error_v21(base_coefficients,gauge_coefficients,base_metric,final_gauge_block):
    cb=np.asarray(base_coefficients,complex); mapped=final_gauge_block.toarray()@np.asarray(gauge_coefficients,complex); inner=np.vdot(cb,base_metric@mapped); phase=1 if abs(inner)<1e-30 else np.exp(-1j*np.angle(inner)); diff=phase*mapped-cb
    return float(np.sqrt(max(float(np.real(np.vdot(diff,base_metric@diff))),0))/np.sqrt(max(float(np.real(np.vdot(cb,base_metric@cb))),1e-30)))
