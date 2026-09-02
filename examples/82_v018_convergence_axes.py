from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v018_convergence_complete_campaign.json").read_text()
)

for name in ["basis","edge_budget","enrich_threshold"]:
    print(f"\n{name}")
    print("-"*len(name))
    for row in data["axes"][name]:
        print(
            row["coordinates"],
            "fidelity=",
            row["projected_wavefunction_fidelity"],
            "L2=",
            row["projected_wavefunction_l2_error"],
        )
