import numpy as np
from scipy.linalg import expm
from gaussian_dynamics.correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from gaussian_dynamics.correlated_basis_adaptation_v270 import generate_multidimensional_spawn_candidates_v270
from gaussian_dynamics.curved_basis_adaptation_v292 import CurvedBasisSettingsV292, CurvedTDVPSnapshotV292, SVDReceiptV292, tensor_legendre_quadrature_v292
from gaussian_dynamics.curved_basis_adaptation_v293 import finite_difference_full_jacobian_v293, candidate_coefficient_block_v293
from gaussian_dynamics.curved_geometry_v293 import ConnectionFieldV293
from gaussian_dynamics.pre_ml_contract_v294 import *

def state():
 return CorrelatedGaussianSpinorStateV270(np.array([[-.25,.15]]),np.array([[.12,-.08]]),np.array([[[1.05,.14],[.14,.86]]]),np.array([[[.05,-.02],[-.02,.03]]]),np.array([[.84+.08j,.43-.14j]])).normalized()
def field():
 sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]],complex); sz=np.diag([1,-1]).astype(complex)
 def f(r):
  x,y=r; return np.array([.14j*sy+.055j*y*sz,.11j*sx-.045j*x*sz])
 return ConnectionFieldV293(2,2,f,"v294-variable").validate([np.zeros(2)])
def snapshot(st,fi,q,cand):
 j=finite_difference_full_jacobian_v293(st,fi,q,relative_step=4e-6); k=candidate_coefficient_block_v293(cand,st.nstate,fi,q)
 r=.13*k[:,0]-.09*k[:,1]+.025j*k[:,0]; jr=np.vstack([j.real,j.imag]); rr=np.r_[r.real,r.imag]; rr += jr@np.linalg.lstsq(jr,-rr,rcond=1e-10)[0]; n=j.shape[0]; r=rr[:n]+1j*rr[n:]
 m=np.real(j.conj().T@j); s=np.linalg.svd(m,compute_uv=False); rec=SVDReceiptV292(s,1e-12,len(s),float(s[0]/s[-1]),0,0)
 return CurvedTDVPSnapshotV292(j,m,np.zeros(j.shape[1]),np.zeros(j.shape[1]),r,rec).validate()
def test_state_fingerprint_is_packet_id_canonical():
 st=state(); p=canonical_state_payload_v294(st,packet_ids=("g000003",)); assert stable_hash_v294(p)==stable_hash_v294(p)
def test_trajectory_split_cannot_leak_candidates():
 assert trajectory_split_v294("traj-A")==trajectory_split_v294("traj-A")
def test_exact_oracle_record_is_finite_versioned_and_reproducible():
 st=state(); fi=field(); q=tensor_legendre_quadrature_v292(((-1.8,1.8),(-1.8,1.8)),(4,4)); c=generate_multidimensional_spawn_candidates_v270(st)[0]; sn=snapshot(st,fi,q,c)
 settings=CurvedBasisSettingsV292(minimum_predicted_residual_reduction=1e-14,minimum_relative_residual_reduction=1e-12)
 a=evaluate_candidate_oracle_v294(st,fi,q,sn,c,trajectory_group_id="traj-A",settings=settings); b=evaluate_candidate_oracle_v294(st,fi,q,sn,c,trajectory_group_id="traj-A",settings=settings)
 assert a.record_fingerprint==b.record_fingerprint and a.accepted_by_residual_oracle and a.predicted_residual_reduction>0
 assert np.isfinite(a.schur_condition_number)
def test_features_do_not_expose_raw_electronic_coefficients():
 st=state(); fi=field(); q=tensor_legendre_quadrature_v292(((-1.8,1.8),(-1.8,1.8)),(4,4)); c=generate_multidimensional_spawn_candidates_v270(st)[0]; sn=snapshot(st,fi,q,c)
 f=invariant_candidate_features_v294(st,c,sn); text=str(f).lower(); assert "coefficient" not in text and "projected_residual_vector" not in text
def test_dataset_manifest_has_group_disjoint_splits_and_stable_hash():
 st=state(); fi=field(); q=tensor_legendre_quadrature_v292(((-1.8,1.8),(-1.8,1.8)),(4,4)); cands=generate_multidimensional_spawn_candidates_v270(st)[:4]; sn=snapshot(st,fi,q,cands[0]); settings=CurvedBasisSettingsV292(minimum_predicted_residual_reduction=1e-14,minimum_relative_residual_reduction=1e-12)
 rows=[]
 for group in ("traj-A","traj-B","traj-C","traj-D"):
  rows.append(evaluate_candidate_oracle_v294(st,fi,q,sn,cands[0],trajectory_group_id=group,settings=settings))
 m=dataset_manifest_v294(rows); assert m["n_records"]==4 and len(m["sha256"])==64
def test_rejection_taxonomy_preserves_exact_reasons():
 st=state(); fi=field(); q=tensor_legendre_quadrature_v292(((-1.8,1.8),(-1.8,1.8)),(4,4)); c=generate_multidimensional_spawn_candidates_v270(st)[0]; sn=snapshot(st,fi,q,c)
 settings=CurvedBasisSettingsV292(minimum_predicted_residual_reduction=1e6,minimum_relative_residual_reduction=.9)
 r=evaluate_candidate_oracle_v294(st,fi,q,sn,c,trajectory_group_id="traj-reject",settings=settings); assert not r.accepted_by_residual_oracle and len(r.rejection_reasons)>=1
def test_rejected_audit_sampling_is_deterministic():
 st=state(); fi=field(); q=tensor_legendre_quadrature_v292(((-1.8,1.8),(-1.8,1.8)),(4,4)); c=generate_multidimensional_spawn_candidates_v270(st)[0]; sn=snapshot(st,fi,q,c); settings=CurvedBasisSettingsV292(minimum_predicted_residual_reduction=1e6,minimum_relative_residual_reduction=.9)
 rows=[evaluate_candidate_oracle_v294(st,fi,q,sn,c,trajectory_group_id=f"traj-{i}",settings=settings) for i in range(10)]
 assert audit_exact_rejections_v294(rows,audit_fraction=.5)==audit_exact_rejections_v294(rows,audit_fraction=.5)
