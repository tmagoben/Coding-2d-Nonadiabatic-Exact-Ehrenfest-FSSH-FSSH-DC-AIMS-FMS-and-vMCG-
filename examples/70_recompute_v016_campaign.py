from pathlib import Path

from gaussian_dynamics import (
    run_v016_release_benchmark,
)
from gaussian_dynamics.campaign_io import (
    save_campaign_json,
)

root=Path(__file__).resolve().parents[1]

result=run_v016_release_benchmark(
    repository_root=root,
)
path=save_campaign_json(
    root/"results"/"v016_sparse_locality_campaign_recomputed.json",
    result,
)

print("Recomputed v0.16 sparse-locality campaign.")
print("Saved:",path)
print("Acceptance:",result["acceptance"])
print(
    "Target density error:",
    result["reference"]["target_density_error"],
)
print(
    "N=80 pair reduction:",
    result["sparse_scaling"][-1]["pair_reduction_fraction"],
)
