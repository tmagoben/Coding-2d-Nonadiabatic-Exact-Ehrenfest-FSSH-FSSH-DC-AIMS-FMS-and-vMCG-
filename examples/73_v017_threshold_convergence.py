from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v017_sparse_error_control_campaign.json").read_text()
)

print("v0.17 edge-score threshold convergence")
print("--------------------------------------")
for row in data["threshold_sweep"]:
    print(
        f"enter={row['enter_score']:.3f} "
        f"edges={row['active_edges']:2d} "
        f"omittedL2={row['omitted_score_l2']:.6g} "
        f"Serr={row['relative_S_frobenius_error']:.6g} "
        f"Herr={row['relative_H_frobenius_error']:.6g}"
    )

print("\nSummary")
print(data["sweep_summary"])
