from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v014_time_adaptive_defect_campaign.json").read_text()
)

c=data["adaptive"]["complexity"]
total=c["total_seconds"]

print("v0.14 algorithmic complexity audit")
print("----------------------------------")
for name in [
    "matrix_build_calls",
    "pair_matrix_evaluations",
    "ordered_pair_equivalent",
    "time_matrix_calls",
    "cayley_solve_calls",
    "defect_evaluations",
    "candidate_ranking_calls",
    "candidate_count_scored",
    "peak_basis_size",
    "peak_electronic_dimension",
    "peak_candidate_count",
]:
    print(f"{name:30s}: {c[name]}")

print("\nMeasured timing")
for name in [
    "matrix_build_seconds",
    "time_matrix_seconds",
    "cayley_solve_seconds",
    "defect_seconds",
    "candidate_ranking_seconds",
    "pruning_seconds",
]:
    value=c[name]
    print(
        f"{name:30s}: {value:.6f} s "
        f"({100*value/total:.2f}%)"
    )

print("\nSymbolic scaling")
for key,value in c["asymptotic"].items():
    print(f"\n{key}:")
    print(value)
