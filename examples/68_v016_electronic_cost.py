from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v016_sparse_locality_campaign.json").read_text()
)
d=data["electronic_cost_demo"]

print("v0.16 electronic-structure cache cost demonstration")
print("---------------------------------------------------")
print("cached geometry:")
print(d["cached_geometry"])
print("\nnew geometry:")
print(d["new_geometry"])

print(
    "\nThe analytic LVC release assigns zero provider cost in production; this "
    "demonstration validates the interface intended for future PySCF timing data."
)
