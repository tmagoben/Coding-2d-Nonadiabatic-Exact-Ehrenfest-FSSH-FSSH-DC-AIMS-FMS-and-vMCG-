from pathlib import Path

from gaussian_dynamics import (
    run_v017_release_benchmark,
)
from gaussian_dynamics.campaign_io import (
    save_campaign_json,
)

root=Path(__file__).resolve().parents[1]

result=run_v017_release_benchmark(
    repository_root=root,
)
path=save_campaign_json(
    root/"results"/"v017_sparse_error_control_campaign_recomputed.json",
    result,
)

print("Recomputed v0.17 campaign.")
print("Saved:",path)
print("Acceptance:",result["acceptance"])
print(
    "Target density error:",
    result["reference"]["target_density_error"],
)
print(
    "Final H audit error:",
    result["final_audit"]["relative_H_frobenius_error"],
)
