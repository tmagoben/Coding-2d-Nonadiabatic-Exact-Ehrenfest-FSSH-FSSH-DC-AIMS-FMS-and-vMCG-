from gaussian_dynamics.fast_lvc_matrices_v14 import (
    ordered_pair_evaluation_count,
    hermitian_pair_evaluation_count,
    pair_evaluation_reduction,
)

print("v0.14 Hermitian S/H pair-count reduction")
print("-----------------------------------------")
for n in (2,5,10,11,25,100):
    full=ordered_pair_evaluation_count(n)
    half=hermitian_pair_evaluation_count(n)
    reduction=pair_evaluation_reduction(n)
    print(
        f"N={n:3d}  ordered={full:5d}  "
        f"half={half:5d}  "
        f"reduction={100*reduction:6.2f}%"
    )

print(
    "\nThe asymptotic class remains O(N^2); this optimization reduces the "
    "leading pair-evaluation count toward one half."
)
