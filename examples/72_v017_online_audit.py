from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v017_sparse_error_control_campaign.json").read_text()
)

print("v0.17 online sparse-matrix audit history")
print("----------------------------------------")
for row in data["adaptive"]["audit_history"]:
    print(
        f"step={row['step']:3d} "
        f"attempt={row['attempt']} "
        f"pass={row['passed']} "
        f"enter={row['enter_score']:.5f} "
        f"Serr={row['relative_S_frobenius_error']:.6g} "
        f"Herr={row['relative_H_frobenius_error']:.6g}"
    )

print("\nRelaxation events")
for event in data["adaptive"]["events"]:
    if event["kind"]=="sparse_audit_relaxation":
        print(event)
