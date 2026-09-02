from pathlib import Path
import json

root=Path(__file__).resolve().parents[1]
data=json.loads(
    (root/"results"/"v015_cost_aware_cache_campaign.json").read_text()
)
c=data["adaptive"]["complexity"]

print("v0.15 complexity ledger")
print("-----------------------")
for key in [
    "full_matrix_builds",
    "incremental_expansions",
    "incremental_prunes",
    "pair_snapshots",
    "pair_requests",
    "pair_factorizations",
    "propagation_pair_factorizations",
    "candidate_pair_factorizations",
    "pair_direct_hits",
    "pair_reverse_views",
    "cache_hit_fraction",
    "v14_factorization_baseline",
    "factorization_avoided",
    "factorization_reduction_fraction",
    "cayley_solve_calls",
    "defect_evaluations",
    "candidate_count_scored",
]:
    print(f"{key:38s}: {c[key]}")

print("\nTiming")
for key in [
    "matrix_build_seconds",
    "time_matrix_seconds",
    "defect_seconds",
    "candidate_ranking_seconds",
    "cost_ranking_seconds",
    "cayley_solve_seconds",
    "total_seconds",
]:
    print(f"{key:38s}: {c[key]:.6f} s")
