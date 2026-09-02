from dataclasses import dataclass
import numpy as np
from .electronic_operator_v21 import ElectronicOperatorPointV21, ElectronicOperatorSnapshotV21

def _hermitian_from_rng(rng, n, scale):
    X = rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n)); H = 0.5*(X+X.conj().T)
    return float(scale)*H/max(np.linalg.norm(H,ord=2),1e-30)

@dataclass(frozen=True)
class SyntheticLinearOperatorConfigV21:
    nstate:int=4; nq:int=2; mass:float=20.0; seed:int=21021; base_scale:float=0.04; derivative_scale:float=0.02

class SyntheticLinearOperatorProviderV21:
    def __init__(self, config=SyntheticLinearOperatorConfigV21()):
        self.config=config; n=int(config.nstate); nq=int(config.nq); rng=np.random.default_rng(int(config.seed))
        self.H0=_hermitian_from_rng(rng,n,config.base_scale)
        self.derivatives=np.asarray([_hermitian_from_rng(rng,n,config.derivative_scale) for _ in range(nq)])
        self.mass_matrix=float(config.mass)*np.eye(nq); self.frame=np.eye(n,dtype=complex); self.calls=0
    def evaluate_snapshot(self,q):
        self.calls+=1; q=np.asarray(q,dtype=float)
        H=self.H0+np.tensordot(q,self.derivatives,axes=(0,0)); D=np.zeros_like(self.derivatives,dtype=complex)
        point=ElectronicOperatorPointV21(q.copy(),H,self.derivatives.copy(),D,self.mass_matrix.copy(),{"backend":"SyntheticLinearOperatorProviderV21","spin_physics":False}).validate()
        return ElectronicOperatorSnapshotV21(point,state_vectors=self.frame.copy(),metadata={"fixed_electronic_frame":True}).validate()
    def evaluate(self,q): return self.evaluate_snapshot(q).point
    def snapshot_overlap(self,left,right): return left.state_vectors.conj().T@right.state_vectors
    def diagnostics_dict(self): return {"calls":int(self.calls),"nstate":int(self.config.nstate),"nq":int(self.config.nq)}
