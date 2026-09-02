from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v018_convergence_complete_campaign.json").read_text()
)
c=data["canonical"]

print("v0.18 sampled omitted-edge audits")
print("---------------------------------")
for row in c["sampled_audits"]:
    print(
        f"step={row['step']:3d} "
        f"sample={row['sample_count']:2d} "
        f"max_score={row['maximum_score']:.8g} "
        f"pass={row['passed']}"
    )

print("\nDense sentinels")
for row in c["sentinel_audits"]:
    print(
        f"{row['label']:7s} "
        f"Serr={row['relative_S_frobenius_error']:.8g} "
        f"Herr={row['relative_H_frobenius_error']:.8g}"
    )
