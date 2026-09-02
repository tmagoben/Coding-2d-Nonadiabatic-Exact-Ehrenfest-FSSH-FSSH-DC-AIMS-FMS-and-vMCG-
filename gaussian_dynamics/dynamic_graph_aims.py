from dataclasses import dataclass, field
import numpy as np
from .dynamic_gauge_graph import IncrementalElectronicGraph, AnalyticCI2DFrameProvider
from .graph_gaussian import GraphGaussianTBF, build_static_graph_gaussian_matrices
from .gaussian_nd import analytic_overlap_equal_width
from .moving_graph_gaussian import nuclear_seed_basis_time_matrix,metric_compatible_basis_connection,moving_basis_coefficient_step,generalized_norm

@dataclass
class DynamicGraphTBF:
    uid:int; state:int; q:np.ndarray; p:np.ndarray; A:np.ndarray; node:object; spawned_targets:set=field(default_factory=set)
    def __post_init__(self):
        self.uid=int(self.uid); self.state=int(self.state); self.q=np.asarray(self.q,float); self.p=np.asarray(self.p,float); self.A=np.asarray(self.A,float)
        if self.q.shape!=self.p.shape: raise ValueError('q and p must have equal shape')
        if self.A.shape!=(len(self.q),len(self.q)): raise ValueError('A incompatible')
    def graph_tbf(self,dimension):
        c=np.zeros(dimension,complex); c[self.state]=1
        return GraphGaussianTBF(self.node,self.q.copy(),self.p.copy(),self.A.copy(),c)

def _center_node(step,uid): return ('tbf',int(uid),int(step))
def _centroid_node(step,ui,uj):
    a,b=sorted((int(ui),int(uj))); return ('centroid',int(step),a,b)

def _verlet(tbf,provider,dt):
    a=provider.evaluate(tbf.q); force=-a.gradients[tbf.state]; ph=tbf.p+0.5*dt*force; q=tbf.q+dt*np.linalg.solve(a.mass_matrix,ph)
    b=provider.evaluate(q); return q, ph+0.5*dt*(-b.gradients[tbf.state])

def _kinematics(tbf,provider):
    a=provider.evaluate(tbf.q); return np.linalg.solve(a.mass_matrix,tbf.p),-a.gradients[tbf.state]

def _indicator(parent,target,provider):
    a=provider.evaluate(parent.q); return float(abs(np.linalg.solve(a.mass_matrix,parent.p)@a.nac[parent.state,target]))

def _child_momentum(parent,target,provider):
    a=provider.evaluate(parent.q); d=a.nac[parent.state,target]; dn=np.linalg.norm(d)
    if dn<1e-14:return None
    n=d/dn; B=np.linalg.inv(a.mass_matrix); de=a.energies[target]-a.energies[parent.state]
    A=float(n@B@n); Bc=float(parent.p@B@n); disc=Bc*Bc-2*A*de
    if disc<0:return None
    r=np.sqrt(disc); roots=[(-Bc+r)/A,(-Bc-r)/A]; lam=min(roots,key=abs); return parent.p+lam*n

def maybe_spawn_dynamic(basis,provider,next_uid,threshold,overlap_block=.85):
    ns=provider.evaluate(basis[0].q).frame.shape[1]
    for parent in basis:
        for target in range(ns):
            if target==parent.state or target in parent.spawned_targets: continue
            if _indicator(parent,target,provider)<=threshold: continue
            pc=_child_momentum(parent,target,provider)
            if pc is None: continue
            redundant=False
            for b in basis:
                if b.state==target and np.allclose(b.A,parent.A,atol=1e-12):
                    ov=abs(analytic_overlap_equal_width(b.q,b.p,parent.q,pc,parent.A))
                    if ov>=overlap_block: redundant=True; break
            if redundant: continue
            parent.spawned_targets.add(target)
            return DynamicGraphTBF(next_uid,target,parent.q.copy(),pc,parent.A.copy(),parent.node),parent.uid
    return None,None

def _ensure_graph(step,basis,manager,provider):
    for b in basis:
        node=_center_node(step,b.uid)
        if node not in manager.frames:
            connect=[b.node] if b.node in manager.frames and b.node!=node else []
            manager.add_from_provider(node,b.q,provider,connect_to=connect)
        b.node=node
    refs={}
    for i in range(len(basis)):
        refs[(i,i)]=basis[i].node
        for j in range(i+1,len(basis)):
            n=_centroid_node(step,basis[i].uid,basis[j].uid)
            if n not in manager.frames:
                manager.add_from_provider(n,0.5*(basis[i].q+basis[j].q),provider,connect_to=[basis[i].node,basis[j].node])
            refs[(i,j)]=n; refs[(j,i)]=n
    return refs

def _gbasis(basis,dimension): return [b.graph_tbf(dimension) for b in basis]

def run_dynamic_graph_aims(initial_basis,C0,provider=None,dt=2e-4,steps=40,spawn_threshold=1e-6,overlap_block=.9,max_basis=4,store_every=5):
    provider=provider or AnalyticCI2DFrameProvider(); dim=provider.evaluate(initial_basis[0].q).frame.shape[1]; manager=IncrementalElectronicGraph(dim)
    basis=[DynamicGraphTBF(b.uid,b.state,b.q.copy(),b.p.copy(),b.A.copy(),b.node,set(b.spawned_targets)) for b in initial_basis]
    C=np.asarray(C0,complex).copy()
    if len(C)!=len(basis): raise ValueError('C0 length mismatch')
    refs=_ensure_graph(0,basis,manager,provider); gb=_gbasis(basis,dim); M=provider.evaluate(basis[0].q).mass_matrix
    S,H=build_static_graph_gaussian_matrices(gb,manager.registry,M,reference_selector=lambda i,j:refs[(i,j)])
    C/=np.sqrt(generalized_norm(C,S)); next_uid=max(b.uid for b in basis)+1; events=[]; records=[]
    def record(step):
        norm=generalized_norm(C,S); pops=np.zeros(dim)
        for st in range(dim):
            idx=[i for i,b in enumerate(basis) if b.state==st]
            if idx:
                block=S[np.ix_(idx,idx)]; cc=C[idx]; pops[st]=np.real(np.vdot(cc,block@cc))
        if norm>0:pops/=norm
        records.append({'step':step,'time':step*dt,'norm':norm,'basis_size':len(basis),'state_populations':pops,'graph':manager.summary(),'condition_number':float(np.linalg.cond(S))})
    record(0)
    for step in range(1,steps+1):
        old=[DynamicGraphTBF(b.uid,b.state,b.q.copy(),b.p.copy(),b.A.copy(),b.node,set(b.spawned_targets)) for b in basis]; oldrefs=refs; oldgb=_gbasis(old,dim); S0,H0=S,H
        for b in basis:b.q,b.p=_verlet(b,provider,dt)
        refs=_ensure_graph(step,basis,manager,provider); gb=_gbasis(basis,dim)
        S1,H1=build_static_graph_gaussian_matrices(gb,manager.registry,M,reference_selector=lambda i,j:refs[(i,j)])
        qd=[];pd=[]
        for b in old:
            a,c=_kinematics(b,provider); qd.append(a);pd.append(c)
        seed=nuclear_seed_basis_time_matrix(oldgb,manager.registry,lambda i,j:oldrefs[(i,j)],np.asarray(qd),np.asarray(pd))
        T=metric_compatible_basis_connection(S0,S1,dt,seed=seed); C=moving_basis_coefficient_step(C,S0,H0,S1,H1,T,dt); S,H=S1,H1
        if len(basis)<max_basis:
            child,parent=maybe_spawn_dynamic(basis,provider,next_uid,spawn_threshold,overlap_block)
            if child is not None:
                basis.append(child); C=np.concatenate([C,[0j]]); events.append({'step':step,'time':step*dt,'parent_uid':parent,'child_uid':next_uid,'target_state':child.state}); next_uid+=1
                refs=_ensure_graph(step,basis,manager,provider); gb=_gbasis(basis,dim); S,H=build_static_graph_gaussian_matrices(gb,manager.registry,M,reference_selector=lambda i,j:refs[(i,j)])
        if step%store_every==0:record(step)
    return {'records':records,'events':events,'final_basis':basis,'final_coefficients':C,'graph':manager.graph,'registry':manager.registry}
