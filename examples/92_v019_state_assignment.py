from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v019_molecular_direct_dynamics_campaign.json").read_text()
)
s=data["state_tracking_scaling"]

print("v0.19 scalable maximum-overlap assignment")
print("-----------------------------------------")
print("nstate:",s["nstate"])
print("valid permutation:",s["permutation_is_valid"])
print("best score:",s["best_score"])
print("second-best score:",s["second_best_score"])
print("diagnostic seconds:",s["assignment_seconds"])
print("legacy complexity:",s["legacy_complexity"])
print("v0.19 complexity:",s["complexity"])
