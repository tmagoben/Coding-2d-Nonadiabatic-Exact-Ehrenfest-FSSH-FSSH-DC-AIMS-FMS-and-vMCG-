from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v015_cost_aware_cache_campaign.json").read_text()
)

old=data["v14_context"]["reference"]
new=data["reference"]
timing=data["timing_comparison"]

print("v0.14 -> v0.15 comparison")
print("-------------------------")
for label,key in [
    ("projection fidelity","projection_fidelity"),
    ("initial density error","initial_density_error"),
    ("projected dynamics error","projected_dynamics_density_error"),
    ("target density error","target_density_error"),
    ("target population error","target_population_error"),
    ("purity","purity"),
    ("coherence phase error","coherence_phase_error"),
]:
    print(
        f"{label:28s}"
        f" v0.14={old[key]:.12g}"
        f" v0.15={new[key]:.12g}"
    )

print("\nmaximum physical metric difference:")
print(data["acceptance"]["v14_reference_difference"]["maximum"])

print("\nsaved timing diagnostic:")
print("v0.14 seconds:",timing["v14_saved_adaptive_seconds"])
print("v0.15 seconds:",timing["v15_adaptive_seconds"])
print("speedup:",timing["saved_benchmark_speedup"])
