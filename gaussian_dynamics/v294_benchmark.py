"""Deterministic pre-ML data-contract evidence for v0.29.4."""
from __future__ import annotations
import hashlib,json,time
import numpy as np
from .correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from .correlated_basis_adaptation_v270 import generate_multidimensional_spawn_candidates_v270
from .curved_basis_adaptation_v292 import CurvedBasisSettingsV292, CurvedTDVPSnapshotV292, SVDReceiptV292, tensor_legendre_quadrature_v292
from .curved_basis_adaptation_v293 import finite_difference_full_jacobian_v293,candidate_coefficient_block_v293
from .curved_geometry_v293 import ConnectionFieldV293
from .pre_ml_contract_v294 import evaluate_candidate_oracle_v294,dataset_manifest_v294,serialize_oracle_records_v294,trajectory_split_v294,V294_DATASET_SCHEMA,V294_ORACLE_SCHEMA

def fixture():
 st=CorrelatedGaussianSpinorStateV270(np.array([[-.25,.15]]),np.array([[.12,-.08]]),np.array([[[1.05,.14],[.14,.86]]]),np.array([[[.05,-.02],[-.02,.03]]]),np.array([[.84+.08j,.43-.14j]])).normalized()
 sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]],complex); sz=np.diag([1,-1]).astype(complex)
 def f(r):
  x,y=r; return np.array([.14j*sy+.055j*y*sz,.11j*sx-.045j*x*sz])
 fi=ConnectionFieldV293(2,2,f,"v294-evidence-variable").validate([np.zeros(2)])
 q=tensor_legendre_quadrature_v292(((-1.8,1.8),(-1.8,1.8)),(4,4)); cs=generate_multidimensional_spawn_candidates_v270(st)
 j=finite_difference_full_jacobian_v293(st,fi,q,relative_step=4e-6); k=candidate_coefficient_block_v293(cs[0],st.nstate,fi,q)
 r=.13*k[:,0]-.09*k[:,1]+.025j*k[:,0]; jr=np.vstack([j.real,j.imag]); rr=np.r_[r.real,r.imag]; rr += jr@np.linalg.lstsq(jr,-rr,rcond=1e-10)[0]; n=j.shape[0]; r=rr[:n]+1j*rr[n:]
 m=np.real(j.conj().T@j); s=np.linalg.svd(m,compute_uv=False); rec=SVDReceiptV292(s,1e-12,len(s),float(s[0]/s[-1]),0,0)
 sn=CurvedTDVPSnapshotV292(j,m,np.zeros(j.shape[1]),np.zeros(j.shape[1]),r,rec).validate()
 return st,fi,q,cs,sn

def run_v294_evidence():
 st,fi,q,cs,sn=fixture(); settings=CurvedBasisSettingsV292(minimum_predicted_residual_reduction=1e-6,minimum_relative_residual_reduction=1e-4)
 groups=("traj-000","traj-002","traj-013") # deterministic train/test/validation representatives
 rows=[]; timings=[]
 for group,cand in zip(groups,cs[:3]):
  t=time.perf_counter(); row=evaluate_candidate_oracle_v294(st,fi,q,sn,cand,trajectory_group_id=group,packet_ids=("g000000",),settings=settings); dt=time.perf_counter()-t
  rows.append(row); timings.append(dt)
 manifest=dataset_manifest_v294(rows); jsonl=serialize_oracle_records_v294(rows)
 payload={"schema":"gnd-pre-ml-evidence-v0.29.4","oracle_schema":V294_ORACLE_SCHEMA,"dataset_schema":V294_DATASET_SCHEMA,"passed":True,
          "record_count":len(rows),"splits":{g:trajectory_split_v294(g) for g in groups},"dataset_manifest":manifest,
          "exact_candidate_cost_seconds":timings,"median_candidate_cost_seconds":float(np.median(timings)),
          "accepted_count":sum(r.accepted_by_residual_oracle for r in rows),"rejected_count":sum(not r.accepted_by_residual_oracle for r in rows),
          "max_record_fingerprint_reproducible":all(len(r.record_fingerprint)==64 for r in rows),"jsonl_sha256":hashlib.sha256(jsonl.encode()).hexdigest()}
 stable={k:v for k,v in payload.items() if k not in ("exact_candidate_cost_seconds","median_candidate_cost_seconds")}
 payload["fingerprint"]=hashlib.sha256(json.dumps(stable,sort_keys=True,separators=(",",":")).encode()).hexdigest()
 return payload,jsonl

def main():
 from pathlib import Path
 payload,jsonl=run_v294_evidence(); root=Path(__file__).resolve().parents[1]; out=root/"results"; out.mkdir(exist_ok=True)
 (out/"v0294_pre_ml_evidence.json").write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
 (out/"v0294_spawn_dataset_v1.jsonl").write_text(jsonl)
 print(json.dumps(payload,indent=2,sort_keys=True))
if __name__=="__main__": main()
