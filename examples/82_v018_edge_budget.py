from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v018_convergence_complete_campaign.json").read_text()
)

print("v0.18 sparse-edge-budget convergence")
print("------------------------------------")
for row in data["edge_budget_axis"]:
    print(
        f"B_local={row['local_score_budget']:.3f}  "
        f"L2={row['wavefunction_projected']['phase_aligned_l2_error']:.9f}  "
        f"sparsity={row['average_graph_sparsity']:.6f}  "
        f"H_audit={row['final_sentinel_H_error']:.6g}"
    )
