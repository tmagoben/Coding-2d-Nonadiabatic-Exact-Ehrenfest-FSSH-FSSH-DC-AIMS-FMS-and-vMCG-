from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v019_molecular_direct_dynamics_campaign.json").read_text()
)
g=data["gauge_graph"]

print("v0.19 molecular center-centroid gauge graph")
print("-------------------------------------------")
for key,value in g.items():
    print(f"{key}: {value}")
