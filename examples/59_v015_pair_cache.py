import numpy as np

from gaussian_dynamics import (
    AnalyticCI2DFrameProvider,
    DynamicGraphTBF,
    GaussianPairCache,
    build_cached_spinor_lvc_matrices,
    build_cached_spinor_time_matrix,
)

provider=AnalyticCI2DFrameProvider(
    nuclear_mass_au=8.0
)

basis=[
    DynamicGraphTBF(
        uid=0,state=1,
        q=np.array([-0.55,0.35]),
        p=np.array([0.7,-0.2]),
        A=np.array([[1.25,0.08],[0.08,0.85]]),
        node=("a",0),
    ),
    DynamicGraphTBF(
        uid=1,state=0,
        q=np.array([0.45,-0.25]),
        p=np.array([-0.15,0.55]),
        A=np.array([[0.75,-0.05],[-0.05,1.15]]),
        node=("b",1),
    ),
]

cache=GaussianPairCache(basis)
S,H,Snuc=build_cached_spinor_lvc_matrices(
    cache,provider
)

solves_after_sh=cache.stats.canonical_solves

qdots=np.array([
    [0.05,0.02],
    [-0.01,0.03],
])
pdots=np.array([
    [-0.02,0.01],
    [0.03,-0.01],
])
T=build_cached_spinor_time_matrix(
    cache,qdots,pdots
)

print("v0.15 shared Gaussian-pair cache")
print("--------------------------------")
print("canonical pair count:",cache.canonical_pair_count)
print("pair solves after S/H:",solves_after_sh)
print("pair solves after T:",cache.stats.canonical_solves)
print("pair requests:",cache.stats.requests)
print("direct hits:",cache.stats.direct_hits)
print("reverse views:",cache.stats.reverse_views)
print("||S-S^dag||:",np.linalg.norm(S-S.conj().T))
print("||H-H^dag||:",np.linalg.norm(H-H.conj().T))
print("T is ordered/non-Hermitian:",not np.allclose(T,T.conj().T))
