from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v019_molecular_direct_dynamics_campaign.json").read_text()
)
d=data["direct_dynamics"]

print("v0.19 spawned molecular direct dynamics")
print("---------------------------------------")
print("spawn events:",d["scrambled_spawn_events"])
print("coefficient difference:",
      d["coefficient_difference"])
print("center difference:",
      d["center_difference"])
print("momentum difference:",
      d["momentum_difference"])
print("norm drift:",
      d["maximum_norm_drift"])
print("final graph audit:",
      d["final_graph_audit"])
