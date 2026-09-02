from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v016_sparse_locality_campaign.json").read_text()
)

old=data["v15_context"]["reference"]
new=data["reference"]
acc=data["acceptance"]

print("v0.15 dense-cache -> v0.16 sparse-locality")
print("-------------------------------------------")
for label,key in [
    ("projected dynamics error","projected_dynamics_density_error"),
    ("target density error","target_density_error"),
    ("population error","target_population_error"),
    ("coherence phase error","coherence_phase_error"),
    ("purity","purity"),
]:
    print(
        f"{label:28s}"
        f" v0.15={old[key]:.12g}"
        f" v0.16={new[key]:.12g}"
    )

print("\nfinal rho difference:",acc["final_rho_difference_vs_v15"])
print(
    "pair factorization reduction:",
    acc["pair_factorization_reduction_vs_v15"],
)
print(
    "average graph sparsity:",
    data["adaptive"]["complexity"]["average_sparsity_fraction"],
)
