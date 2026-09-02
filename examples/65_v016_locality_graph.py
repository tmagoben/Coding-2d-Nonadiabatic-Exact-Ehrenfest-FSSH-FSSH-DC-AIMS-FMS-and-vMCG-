import numpy as np

from gaussian_dynamics import (
    DynamicGraphTBF,
    LocalityGraphSettings,
    PersistentGaussianLocalityGraph,
)


def tbf(uid,x):
    return DynamicGraphTBF(
        uid=uid,state=0,
        q=np.array([x,0.0]),
        p=np.zeros(2),
        A=np.eye(2),
        node=("locality",uid),
    )


basis=[
    tbf(0,0.0),
    tbf(1,1.5),
    tbf(2,3.0),
    tbf(3,7.5),
]

graph=PersistentGaussianLocalityGraph(
    LocalityGraphSettings(
        enter_overlap=0.03,
        exit_overlap=0.015,
    )
)
update=graph.update(basis)

print("v0.16 persistent Gaussian locality graph")
print("----------------------------------------")
print("active edges:",update.active_edges)
print("all off-diagonal pairs:",update.total_offdiagonal_pairs)
print("KD-tree spatial candidates:",update.spatial_candidate_pairs)
print("globally screened:",update.globally_screened_pairs)
print("exact overlap checks:",update.exact_pair_checks)
print("sparsity fraction:",update.sparsity_fraction)
