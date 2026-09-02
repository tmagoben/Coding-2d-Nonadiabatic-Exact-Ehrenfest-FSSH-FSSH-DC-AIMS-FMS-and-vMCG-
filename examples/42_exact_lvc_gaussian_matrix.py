import numpy as np

from gaussian_dynamics import (
    AnalyticCI2DFrameProvider,
    DynamicGraphTBF,
    build_exact_lvc_gaussian_matrices,
)

provider=AnalyticCI2DFrameProvider(nuclear_mass_au=8.0)

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

S,H=build_exact_lvc_gaussian_matrices(basis,provider)

print("v0.12 exact analytic LVC Gaussian matrices")
print("-------------------------------------------")
print("S =")
print(S)
print("\nH =")
print(H)
print("\n||S-S^dag|| =",np.linalg.norm(S-S.conj().T))
print("||H-H^dag|| =",np.linalg.norm(H-H.conj().T))
print("cond(S) =",np.linalg.cond(S))
