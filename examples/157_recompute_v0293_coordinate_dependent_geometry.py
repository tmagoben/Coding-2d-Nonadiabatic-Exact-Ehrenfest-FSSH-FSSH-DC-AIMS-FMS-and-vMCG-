from pathlib import Path
from gaussian_dynamics.v293_benchmark import save_v293_evidence
p=Path(__file__).resolve().parents[1]/"results"/"v0293_coordinate_dependent_geometry_evidence.json"
e=save_v293_evidence(p)
print("passed:",e["passed"])
print("fingerprint:",e["fingerprint"])
for k,v in e["metrics"].items(): print(f"{k}: {v}")
