from pathlib import Path

from gaussian_dynamics import (
    run_v015_release_benchmark,
)
from gaussian_dynamics.campaign_io import (
    save_campaign_json,
)

root=Path(__file__).resolve().parents[1]

result=run_v015_release_benchmark(
    repository_root=root,
)
path=save_campaign_json(
    root/"results"/"v015_cost_aware_cache_campaign_recomputed.json",
    result,
)

print("Recomputed v0.15 campaign.")
print("Saved:",path)
print("Acceptance:",result["acceptance"])
print("Target density error:",result["reference"]["target_density_error"])
print(
    "Pair-factorization reduction:",
    result["adaptive"]["complexity"]["factorization_reduction_fraction"],
)
