from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v018_convergence_complete_campaign.json").read_text()
)
c=data["canonical"]

print("v0.18 full-wavefunction benchmark")
print("---------------------------------")
print("projected-reference fidelity:",
      c["wavefunction_projected"]["fidelity"])
print("phase-aligned L2:",
      c["wavefunction_projected"]["phase_aligned_l2_error"])
print("nuclear-density L2:",
      c["wavefunction_projected"]["nuclear_density_l2_error"])
print("nuclear-density TV:",
      c["wavefunction_projected"]["nuclear_density_total_variation"])
print("centroid error:",
      c["wavefunction_projected"]["mean_error_l2"])
print("covariance error:",
      c["wavefunction_projected"]["covariance_error_frobenius"])
print("reduced-density error:",
      c["reduced_density_projected"]["density_frobenius_error"])
