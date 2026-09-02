from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v018_convergence_complete_campaign.json").read_text()
)

print("v0.18 coarse -> medium -> fine refinement")
print("-----------------------------------------")
for label,row in zip(
    ["coarse","medium","fine"],
    data["refinement_ladder"]["rows"],
):
    print(
        f"{label:6s}",
        row["coordinates"],
        "fidelity=",
        row["projected_wavefunction_fidelity"],
        "L2=",
        row["projected_wavefunction_l2_error"],
    )

print("\nsummary:")
print(data["refinement_ladder"]["summary"])
