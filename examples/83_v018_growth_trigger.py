from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v018_convergence_complete_campaign.json").read_text()
)

print("v0.18 adaptive-growth threshold sensitivity")
print("-------------------------------------------")
for row in data["growth_trigger_axis"]:
    print(
        f"threshold={row['enrich_relative_threshold']:.3f}  "
        f"N={row['final_basis_size']:2d}  "
        f"events={row['enrichment_steps']}  "
        f"L2={row['wavefunction_projected']['phase_aligned_l2_error']:.9f}"
    )
