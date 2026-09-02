from pathlib import Path

from gaussian_dynamics.v14_benchmark import (
    run_v014_release_benchmark,
)
from gaussian_dynamics.campaign_io import (
    save_campaign_json,
)

root=Path(__file__).resolve().parents[1]

result=run_v014_release_benchmark(
    repository_root=root,
)
path=save_campaign_json(
    root/"results"/"v014_time_adaptive_defect_campaign_recomputed.json",
    result,
)

print("Recomputed v0.14 campaign.")
print("Saved:",path)
print("Acceptance:",result["acceptance"])
print(
    "Target density error:",
    result["reference"]["target_density_error"],
)
print(
    "Projected-state dynamics error:",
    result["reference"]["projected_dynamics_density_error"],
)
