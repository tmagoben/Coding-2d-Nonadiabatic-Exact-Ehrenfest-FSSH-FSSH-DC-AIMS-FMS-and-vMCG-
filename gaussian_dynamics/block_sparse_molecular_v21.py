from dataclasses import dataclass
import math
import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from .gaussian_general import validate_spd
from .gaussian_general import (
    gaussian_overlap_general, kinetic_matrix_element_general,
    basis_time_matrix_element_general, real_overlap_saddle_point,
)
from .gauge_graph import nearest_unitary
from .sparse_molecular_matrices_v20 import _position_bound, _safe_radius

@dataclass
class BlockMolecularTBFV21:
    uid:int; q:np.ndarray; p:np.ndarray; A:np.ndarray
    def __post_init__(self):
        self.uid=int(self.uid); self.q=np.asarray(self.q,float); self.p=np.asarray(self.p,float); self.A=validate_spd(self.A)
        if self.q.shape!=self.p.shape or self.A.shape!=(len(self.q),len(self.q)): raise ValueError("incompatible TBF shapes")
    def copy(self): return BlockMolecularTBFV21(self.uid,self.q.copy(),self.p.copy(),self.A.copy())

@dataclass(frozen=True)
class BlockSparseSettingsV21:
    enter_score:float=0.030; exit_score:float=0.015; search_overlap_floor:float=1e-5
    overlap_weight:float=1.0; hamiltonian_weight:float=0.20; nuclear_time_weight:float=1.0; subspace_mismatch_weight:float=0.20
    local_omitted_score_l2_budget:float=0.010; hamiltonian_floor:float=1e-10; use_kdtree:bool=True
    def validate(self):
        if not (0<self.exit_score<=self.enter_score): raise ValueError("Require 0 < exit_score <= enter_score")
        if not (0<self.search_overlap_floor<1): raise ValueError("bad search floor")
        if self.local_omitted_score_l2_budget<0: raise ValueError("bad omission budget")
        return self

@dataclass(frozen=True)
class BlockPairDataV21:
    i:int; j:int; nuclear_overlap:complex; electronic_overlap:np.ndarray; S_block:np.ndarray; H_block:np.ndarray
    T_ij_block:np.ndarray; T_ji_block:np.ndarray; centroid_q:np.ndarray
    score_overlap:float; score_hamiltonian:float; score_nuclear_time:float; score_subspace_mismatch:float; score:float

@dataclass
class BlockSparseUpdateV21:
    nstate:int; active_edges:tuple; pair_data:dict; diagonal_S:tuple; diagonal_H:tuple; diagonal_T:tuple
    qdots:np.ndarray; pdots:np.ndarray; center_snapshots:tuple; exact_pair_checks:int; spatial_candidate_pairs:int
    budget_promoted_edges:int; omitted_score_l2:float; total_offdiagonal_pairs:int; entered_edges:int; exited_edges:int; retained_edges:int
    @property
    def sparsity_fraction(self): return float(1-len(self.active_edges)/max(self.total_offdiagonal_pairs,1))

@dataclass
class BlockSparseMatricesV21:
    S:sparse.csr_matrix; H:sparse.csr_matrix; T_seed:sparse.csr_matrix; active_edges:tuple; n_basis:int; nstate:int; update:BlockSparseUpdateV21
    @property
    def dimension(self): return self.n_basis*self.nstate


def _uid_edge(a,b): return (int(a),int(b)) if int(a)<int(b) else (int(b),int(a))

def _append_block(rows,cols,data,i,j,block,s):
    B=np.asarray(block,complex)
    for a in range(s):
        for b in range(s):
            if B[a,b]!=0:
                rows.append(s*i+a); cols.append(s*j+b); data.append(B[a,b])

def _block_scale(A):
    A=np.asarray(A,complex); return float(np.linalg.norm(A,'fro')/math.sqrt(max(A.shape[0],1)))

def _electronic_connection_velocity(snapshot,qdot): return np.tensordot(np.asarray(qdot,float),np.asarray(snapshot.point.connection_q,complex),axes=(0,0))

def _center_data(basis,provider,qdots,pdots):
    basis=list(basis); qdots=np.asarray(qdots,float); pdots=np.asarray(pdots,float)
    snaps=[]; diagS=[]; diagH=[]; diagT=[]; ns=None
    for i,b in enumerate(basis):
        snap=provider.evaluate_snapshot(b.q); ns=snap.point.nstate if ns is None else ns
        if snap.point.nstate!=ns: raise ValueError("electronic dimension changed")
        M=snap.point.mass_matrix_q_au
        Tn=kinetic_matrix_element_general(b.q,b.p,b.A,b.q,b.p,b.A,M)
        tau=basis_time_matrix_element_general(b.q,b.p,b.A,b.q,b.p,b.A,qdots[i],pdots[i])
        Gamma=_electronic_connection_velocity(snap,qdots[i]); I=np.eye(ns,dtype=complex)
        diagS.append(I); diagH.append(Tn*I+snap.point.H); diagT.append(tau*I+Gamma); snaps.append(snap)
    return ns,tuple(snaps),tuple(diagS),tuple(diagH),tuple(diagT)

def block_pair_data_v21(basis,i,j,provider,dt,qdots,pdots,center_snapshots,diagonal_H,settings):
    bi,bj=basis[i],basis[j]; si,sj=center_snapshots[i],center_snapshots[j]; ns=si.point.nstate
    Oij=np.asarray(provider.snapshot_overlap(si,sj),complex)
    qbar=real_overlap_saddle_point(bi.q,bi.A,bj.q,bj.A); sc=provider.evaluate_snapshot(qbar)
    Uci=nearest_unitary(np.asarray(provider.snapshot_overlap(sc,si),complex)); Ucj=nearest_unitary(np.asarray(provider.snapshot_overlap(sc,sj),complex))
    He=Uci.conj().T@sc.point.H@Ucj
    Sn=gaussian_overlap_general(bi.q,bi.p,bi.A,bj.q,bj.p,bj.A)
    Tn=kinetic_matrix_element_general(bi.q,bi.p,bi.A,bj.q,bj.p,bj.A,sc.point.mass_matrix_q_au)
    S=Sn*Oij; H=Tn*Oij+Sn*He
    tauij=basis_time_matrix_element_general(bi.q,bi.p,bi.A,bj.q,bj.p,bj.A,qdots[j],pdots[j])
    tauji=basis_time_matrix_element_general(bj.q,bj.p,bj.A,bi.q,bi.p,bi.A,qdots[i],pdots[i])
    Gi=_electronic_connection_velocity(si,qdots[i]); Gj=_electronic_connection_velocity(sj,qdots[j])
    Oji=np.asarray(provider.snapshot_overlap(sj,si),complex)
    Tij=tauij*Oij+Sn*(Oij@Gj); Tji=tauji*Oji+np.conj(Sn)*(Oji@Gi)
    snorm=_block_scale(S); hnorm=_block_scale(H)
    hi=max(_block_scale(diagonal_H[i]),settings.hamiltonian_floor); hj=max(_block_scale(diagonal_H[j]),settings.hamiltonian_floor)
    hrel=hnorm/max(math.sqrt(hi*hj),settings.hamiltonian_floor)
    overlap_scale=_block_scale(Oij)
    nuclear_time=dt*math.sqrt(abs(tauij)**2+abs(tauji)**2)*overlap_scale/math.sqrt(2.0)
    sv=np.clip(np.real(np.linalg.svd(Oij,compute_uv=False)),0,1)
    mismatch=float(math.sqrt(np.mean(np.maximum(1-sv**2,0))))
    score=float(math.sqrt((settings.overlap_weight*snorm)**2+(settings.hamiltonian_weight*hrel)**2+(settings.nuclear_time_weight*nuclear_time)**2+(settings.subspace_mismatch_weight*mismatch)**2))
    return BlockPairDataV21(i,j,Sn,Oij,S,H,Tij,Tji,qbar.copy(),snorm,hrel,nuclear_time,mismatch,score)

class BlockSparseMolecularGraphV21:
    def __init__(self,provider,dt,settings=BlockSparseSettingsV21()):
        self.provider=provider; self.dt=float(dt); self.settings=settings.validate(); self._active_uid_edges=set(); self.total_exact_pair_checks=0
    def active_uid_edges_v214(self):
        """Return the hysteretic graph state needed for a deterministic restart."""
        return tuple(sorted(self._active_uid_edges))
    def restore_active_uid_edges_v214(self,edges,live_uids=None):
        """Restore validated uid edges without exposing mutable graph internals."""
        live=None if live_uids is None else {int(uid) for uid in live_uids}
        restored=set()
        for raw in edges:
            if len(raw)!=2:
                raise ValueError("every restored graph edge must contain two uids.")
            a,b=raw
            ia,ib=int(a),int(b)
            if ia!=a or ib!=b or ia==ib:
                raise ValueError("restored graph edges require distinct integer uids.")
            edge=_uid_edge(ia,ib)
            if live is not None and not set(edge)<=live:
                raise ValueError("a restored graph edge references a non-live Gaussian uid.")
            if edge in restored:
                raise ValueError("restored graph edges must be unique.")
            restored.add(edge)
        self._active_uid_edges=restored
        return self.active_uid_edges_v214()
    def update(self,basis,qdots,pdots):
        basis=list(basis); uids=[b.uid for b in basis]; uid_to_index={u:i for i,u in enumerate(uids)}; live=set(uids)
        old={e for e in self._active_uid_edges if e[0] in live and e[1] in live}
        q=np.asarray([b.q for b in basis],float); amin=np.asarray([np.min(np.linalg.eigvalsh(validate_spd(b.A))) for b in basis])
        ns,centers,diagS,diagH,diagT=_center_data(basis,self.provider,qdots,pdots)
        n=len(basis); total=n*(n-1)//2
        if self.settings.use_kdtree and n>1:
            cand={tuple(sorted(x)) for x in cKDTree(q).query_pairs(_safe_radius(float(np.min(amin)),self.settings.search_overlap_floor),output_type='set')}
        else: cand={(i,j) for i in range(n) for j in range(i+1,n)}
        for a,b in old: cand.add(tuple(sorted((uid_to_index[a],uid_to_index[b]))))
        spatial=len(cand); filtered=[]
        for i,j in sorted(cand):
            e=_uid_edge(uids[i],uids[j])
            if e in old or _position_bound(q[i],q[j],amin[i],amin[j])>=self.settings.search_overlap_floor: filtered.append((i,j))
        pdata={}; active=set(); omitted=[]; entered=retained=0
        for i,j in filtered:
            e=_uid_edge(uids[i],uids[j]); p=block_pair_data_v21(basis,i,j,self.provider,self.dt,np.asarray(qdots,float),np.asarray(pdots,float),centers,diagH,self.settings); pdata[(i,j)]=p
            threshold=self.settings.exit_score if e in old else self.settings.enter_score
            if p.score>=threshold:
                active.add(e); retained += int(e in old); entered += int(e not in old)
            else: omitted.append((p.score,e,(i,j)))
        sq=sum(float(x[0])**2 for x in omitted); promoted=0; budget=self.settings.local_omitted_score_l2_budget
        if math.sqrt(sq)>budget:
            for score,e,ij in sorted(omitted,key=lambda x:(-x[0],x[1])):
                if math.sqrt(max(sq,0))<=budget: break
                active.add(e); sq-=float(score)**2; promoted+=1; retained+=int(e in old); entered+=int(e not in old)
        self._active_uid_edges=active; self.total_exact_pair_checks+=len(filtered)
        active_indices=tuple(sorted((min(uid_to_index[a],uid_to_index[b]),max(uid_to_index[a],uid_to_index[b])) for a,b in active))
        return BlockSparseUpdateV21(ns,active_indices,pdata,diagS,diagH,diagT,np.asarray(qdots,float).copy(),np.asarray(pdots,float).copy(),centers,len(filtered),spatial,promoted,float(math.sqrt(max(sq,0))),total,entered,len(old-active),retained)

def build_block_sparse_matrices_v21(basis,update):
    n=len(basis); s=update.nstate; rs=[];cs=[];ds=[]; rh=[];ch=[];dh=[]; rt=[];ct=[];dt=[]
    for i in range(n):
        _append_block(rs,cs,ds,i,i,update.diagonal_S[i],s); _append_block(rh,ch,dh,i,i,update.diagonal_H[i],s); _append_block(rt,ct,dt,i,i,update.diagonal_T[i],s)
    for i,j in update.active_edges:
        p=update.pair_data[(i,j)]
        _append_block(rs,cs,ds,i,j,p.S_block,s); _append_block(rs,cs,ds,j,i,p.S_block.conj().T,s)
        _append_block(rh,ch,dh,i,j,p.H_block,s); _append_block(rh,ch,dh,j,i,p.H_block.conj().T,s)
        _append_block(rt,ct,dt,i,j,p.T_ij_block,s); _append_block(rt,ct,dt,j,i,p.T_ji_block,s)
    shape=(n*s,n*s)
    return BlockSparseMatricesV21(sparse.coo_matrix((ds,(rs,cs)),shape=shape).tocsr(),sparse.coo_matrix((dh,(rh,ch)),shape=shape).tocsr(),sparse.coo_matrix((dt,(rt,ct)),shape=shape).tocsr(),update.active_edges,n,s,update)

def build_dense_block_reference_v21(basis,provider,dt,qdots,pdots,settings=BlockSparseSettingsV21()):
    basis=list(basis); ns,centers,diagS,diagH,diagT=_center_data(basis,provider,qdots,pdots); n=len(basis); dim=n*ns
    S=np.zeros((dim,dim),complex); H=np.zeros((dim,dim),complex); T=np.zeros((dim,dim),complex)
    def put(A,i,j,B): A[ns*i:ns*(i+1),ns*j:ns*(j+1)]=B
    for i in range(n): put(S,i,i,diagS[i]); put(H,i,i,diagH[i]); put(T,i,i,diagT[i])
    pdata={}
    for i in range(n):
        for j in range(i+1,n):
            p=block_pair_data_v21(basis,i,j,provider,float(dt),np.asarray(qdots,float),np.asarray(pdots,float),centers,diagH,settings); pdata[(i,j)]=p
            put(S,i,j,p.S_block); put(S,j,i,p.S_block.conj().T); put(H,i,j,p.H_block); put(H,j,i,p.H_block.conj().T); put(T,i,j,p.T_ij_block); put(T,j,i,p.T_ji_block)
    return {'S':S,'H':H,'T_seed':T,'pair_data':pdata,'nstate':ns,'pair_count':n*(n-1)//2,'center_snapshots':centers}

def block_diagonal_gauge_v21(gauges): return sparse.block_diag([sparse.csr_matrix(np.asarray(G,complex)) for G in gauges],format='csr')
