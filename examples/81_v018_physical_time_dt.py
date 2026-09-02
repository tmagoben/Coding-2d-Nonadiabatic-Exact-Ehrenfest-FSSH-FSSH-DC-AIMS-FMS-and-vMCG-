from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v018_convergence_complete_campaign.json").read_text()
)

print("v0.18 physical-time-normalized dt surface")
print("-----------------------------------------")
for row in data["axes"]["dt"]:
    r=row["resolved_control_steps"]
    print(
        f"dt={row['coordinates']['dt']:.4f} "
        f"defect_every={r['defect_interval']:3d} steps "
        f"audit_every={r['sampled_audit_interval']:3d} steps "
        f"fidelity={row['projected_wavefunction_fidelity']:.12f}"
    )
