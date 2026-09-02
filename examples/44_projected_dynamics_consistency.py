from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v012_representation_consistent_campaign.json").read_text()
)
r=data["reference_case"]

print("v0.12 representation-consistent nine-Gaussian benchmark")
print("--------------------------------------------------------")
print("initial projection fidelity:          ",r["initial_projection_fidelity"])
print("initial reduced-density error:        ",r["initial_density_error"])
print("projected-state dynamics error:       ",r["projected_dynamics_density_error"])
print("original-target final density error:  ",r["target_density_error"])
print("population L2 error:                  ",r["target_population_error"])
print("trace distance:                       ",r["target_trace_distance"])
print("purity error:                         ",r["purity_error"])
print("coherence phase error / rad:          ",r["coherence_phase_error"])
print("maximum generalized norm drift:       ",r["max_norm_drift"])

print(
    "\nThe dynamics error compares Gaussian propagation with exact TDSE propagation "
    "from the identical projected initial state."
)
