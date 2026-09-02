from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v019_molecular_direct_dynamics_campaign.json").read_text()
)

print("v0.19 PySCF molecular bridge status")
print("-----------------------------------")
for key,value in data["pyscf"].items():
    print(f"{key}: {value}")

print(
    "\nUse PySCFRawSnapshotBackendV19 with "
    "pyscf_snapshot_overlap_engine_v19 in an environment where PySCF is installed."
)
