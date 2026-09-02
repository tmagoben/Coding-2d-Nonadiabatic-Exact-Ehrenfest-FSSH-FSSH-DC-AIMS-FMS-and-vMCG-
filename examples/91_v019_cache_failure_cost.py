from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v019_molecular_direct_dynamics_campaign.json").read_text()
)

print("v0.19 cache / failure / cost policy")
print("-----------------------------------")
print("provider cost:")
for key,value in data["provider_cost"].items():
    print(key,value)

print("\nfailure fallback diagnostics:")
diag=data["failure_fallback"]["diagnostics"]
print("backend failures:",diag["backend_failures"])
print("fallback uses:",diag["fallback_uses"])
print("history:",diag["history"])
