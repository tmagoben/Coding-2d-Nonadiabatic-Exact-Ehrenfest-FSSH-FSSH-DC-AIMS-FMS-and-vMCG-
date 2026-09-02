from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v020_sparse_molecular_campaign.json").read_text()
)
ctrl=data["controller_demo"]

print("v0.20 sampled molecular search controller")
print("-----------------------------------------")
print("relaxations:",ctrl["relaxation_events"])
print("audit history:")
for row in ctrl["audit_history"]:
    print(row)
print("final audit passed:",ctrl["final_audit_passed"])
