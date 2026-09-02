import numpy as np

from gaussian_dynamics import (
    AnalyticCI2DFrameProvider,
    DynamicGraphTBF,
    GaussianPairCache,
    build_cached_spinor_lvc_matrices,
    expand_cached_spinor_lvc_matrices,
)

provider=AnalyticCI2DFrameProvider(
    nuclear_mass_au=8.0
)

basis=[
    DynamicGraphTBF(
        uid=0,state=1,
        q=np.array([-0.5,0.3]),
        p=np.array([0.5,0.1]),
        A=np.eye(2),
        node=("a",0),
    ),
    DynamicGraphTBF(
        uid=1,state=0,
        q=np.array([0.4,-0.2]),
        p=np.array([-0.1,0.4]),
        A=1.3*np.eye(2),
        node=("b",1),
    ),
]
child=DynamicGraphTBF(
    uid=2,state=1,
    q=np.array([0.1,0.6]),
    p=np.array([0.2,-0.2]),
    A=0.8*np.eye(2),
    node=("child",2),
)

cache=GaussianPairCache(basis)
S,H,Snuc=build_cached_spinor_lvc_matrices(
    cache,provider
)

expanded=cache.expanded(child)
# Mimic candidate conditioning: compute the child pairs before acceptance.
for i in range(3):
    expanded.pair(i,2)

before=expanded.stats.canonical_solves
S2,H2,N2=expand_cached_spinor_lvc_matrices(
    S,H,Snuc,expanded,provider
)
after=expanded.stats.canonical_solves

print("v0.15 incremental accepted-child matrix expansion")
print("------------------------------------------------")
print("old basis:",2)
print("new basis:",3)
print("child pairs already cached:",3)
print("new pair solves during matrix expansion:",after-before)
print("new S dimension:",S2.shape)
print("new H Hermiticity error:",np.linalg.norm(H2-H2.conj().T))
