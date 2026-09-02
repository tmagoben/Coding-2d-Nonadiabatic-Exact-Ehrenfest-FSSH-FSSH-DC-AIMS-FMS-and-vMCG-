from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v020_sparse_molecular_campaign.json").read_text()
)

print("v0.20 PySCF sparse-molecular status")
print("----------------------------------")
for key,value in data["pyscf"].items():
    print(f"{key}: {value}")

print(
    "\nSee V20_PYSCF_SPARSE_MOLECULAR_PROTOCOL.md "
    "for the real-backend integration path."
)
