from pathlib import Path

from gaussian_dynamics import (
    run_v020_release_benchmark,
)
from gaussian_dynamics.campaign_io import (
    save_campaign_json,
)

root=Path(__file__).resolve().parents[1]
out=run_v020_release_benchmark()
path=save_campaign_json(
    root/"results"/"v020_sparse_molecular_campaign_recomputed.json",
    out,
)

print("Recomputed v0.20 campaign.")
print("Saved:",path)
print("Acceptance:",out["acceptance"])
