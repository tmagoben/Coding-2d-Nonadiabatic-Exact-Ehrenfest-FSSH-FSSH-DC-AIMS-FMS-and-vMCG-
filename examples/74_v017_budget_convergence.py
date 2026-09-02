from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v017_sparse_error_control_campaign.json").read_text()
)

print("v0.17 local importance-budget convergence")
print("-----------------------------------------")
for row in data["budget_sweep"]:
    print(
        f"budget={row['budget']:.4g} "
        f"edges={row['active_edges']:2d} "
        f"promoted={row['budget_promoted_edges']:2d} "
        f"omittedL2={row['omitted_score_l2']:.6g} "
        f"Serr={row['relative_S_frobenius_error']:.6g} "
        f"Herr={row['relative_H_frobenius_error']:.6g}"
    )
