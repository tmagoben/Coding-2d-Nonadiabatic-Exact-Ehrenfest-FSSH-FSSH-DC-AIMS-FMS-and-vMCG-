from pathlib import Path

from gaussian_dynamics import (
    run_v018_release_benchmark,
)
from gaussian_dynamics.campaign_io import (
    save_campaign_json,
)

root=Path(__file__).resolve().parents[1]

print(
    "This recomputes the complete v0.18 multi-axis campaign in one process. "
    "For release-grade timing isolation, prefer running individual "
    "ConvergenceCoordinatesV18 points with run_coordinate_worker_v18."
)

result=run_v018_release_benchmark(
    repository_root=root
)
path=save_campaign_json(
    root/"results"/"v018_convergence_complete_campaign_recomputed.json",
    result,
)

print("Saved:",path)
print("Acceptance:",result["acceptance"])
