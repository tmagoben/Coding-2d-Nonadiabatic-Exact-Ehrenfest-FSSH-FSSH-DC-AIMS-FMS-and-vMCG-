from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v016_sparse_locality_campaign.json").read_text()
)
a=data["final_sparse_matrix_audit"]

print("v0.16 dense endpoint sparse-matrix audit")
print("---------------------------------------")
for key,value in a.items():
    print(f"{key:38s}: {value}")

print(
    "\nThe overlap cutoff is rigorous for |S_ij| screening, but the dense H audit "
    "is retained because small overlap alone is not a universal Hamiltonian-error bound."
)
