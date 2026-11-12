import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import json, hashlib
import numpy as np
from gaussian_dynamics.correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from gaussian_dynamics.correlated_basis_adaptation_v270 import generate_multidimensional_spawn_candidates_v270
from gaussian_dynamics.curved_basis_adaptation_v292 import *

state=CorrelatedGaussianSpinorStateV270(
 q=np.array([[-.35,.20]]),p=np.array([[.18,-.11]]),
 width_matrices=np.array([[[1.10,.18],[.18,.82]]]),chirp_matrices=np.array([[[.07,-.03],[-.03,.04]]]),
 coefficients=np.array([[.82+.10j,.42-.16j]])).normalized()
sy=np.array([[0,-1j],[1j,0]],complex); sz=np.array([[1,0],[0,-1]],complex)
conn=ConstantConnectionV292(np.array([.22j*sy,.17j*sz])).validate()
quad=tensor_legendre_quadrature_v292(((-4,4),(-4,4)),(18,18))
cand=generate_multidimensional_spawn_candidates_v270(state)[0]
j=finite_difference_full_jacobian_v292(state,conn,quad,relative_step=3e-6)
k=candidate_coefficient_block_v292(cand,state.nstate,conn,quad)
res=.16*k[:,0]-.11*k[:,1]+.03j*k[:,0]
jr=np.vstack([j.real,j.imag]); rr=np.concatenate([res.real,res.imag])
correction=np.linalg.lstsq(jr,-rr,rcond=1e-10)[0]; rr=rr+jr@correction
n=j.shape[0]; res=rr[:n]+1j*rr[n:]
m=np.real(j.conj().T@j); s=np.linalg.svd(m,compute_uv=False)
snap=CurvedTDVPSnapshotV292(j,m,np.zeros(j.shape[1]),np.zeros(j.shape[1]),res,SVDReceiptV292(s,1e-12,len(s),float(s[0]/s[-1]),0,0)).validate()
settings=CurvedBasisSettingsV292(minimum_predicted_residual_reduction=1e-14,minimum_relative_residual_reduction=1e-12)
ev=evaluate_curved_spawn_candidate_v292(state,conn,quad,snap,cand,settings=settings)
actual=actual_linear_residual_reduction_v292(snap,k)
out,proj=spawn_curved_candidate_v292(state,cand,conn,quad)
payload={"schema":V292_SCHEMA,"predicted_reduction":ev.predicted_residual_reduction,"actual_linear_reduction":actual,"relative_difference":abs(ev.predicted_residual_reduction-actual)/actual,"schur_rank":ev.schur_rank,"schur_condition":ev.schur_condition_number,"spawn_projection_loss":proj.relative_loss,"spawn_fidelity":proj.fidelity,"newborn_coefficient_norm":float(np.linalg.norm(out.coefficients[-1])),"no_event_identity":no_event_reduction_v292(state) is state}
payload["fingerprint"]=checkpoint_fingerprint_v292(payload)
print(json.dumps(payload,indent=2,sort_keys=True))
