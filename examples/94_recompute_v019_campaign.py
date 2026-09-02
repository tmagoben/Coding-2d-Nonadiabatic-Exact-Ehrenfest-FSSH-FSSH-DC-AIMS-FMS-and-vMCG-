from pathlib import Path

from gaussian_dynamics import (
    run_v019_release_benchmark,
)
from gaussian_dynamics.campaign_io import (
    save_campaign_json,
)

root=Path(__file__).resolve().parents[1]

result=run_v019_release_benchmark()
path=save_campaign_json(
    root/"results"/"v019_molecular_direct_dynamics_campaign_recomputed.json",
    result,
)

print("Recomputed v0.19 campaign.")
print("Saved:",path)
print("Acceptance:",result["acceptance"])
