import numpy as np

from gaussian_dynamics import (
    AnalyticCI2DFrameProvider,
    DynamicGraphTBF,
    EdgeImportanceSettingsV17,
    ErrorControlledGaussianLocalityGraphV17,
)

provider=AnalyticCI2DFrameProvider(
    nuclear_mass_au=20.0
)

basis=[
    DynamicGraphTBF(
        uid=0,state=1,
        q=np.array([-0.7,0.25]),
        p=np.array([0.8,0.0]),
        A=1.2*np.eye(2),
        node=("v17",0),
    ),
    DynamicGraphTBF(
        uid=1,state=0,
        q=np.array([0.4,0.25]),
        p=np.array([-0.4,0.2]),
        A=1.0*np.eye(2),
        node=("v17",1),
    ),
]

graph=ErrorControlledGaussianLocalityGraphV17(
    provider,
    dt=0.005,
    settings=EdgeImportanceSettingsV17(
        enter_score=0.03,
        exit_score=0.015,
        search_overlap_floor=1e-5,
    ),
)
update=graph.update(basis)
info=update.importance[(0,1)]

print("v0.17 exact local S/H/T edge importance")
print("---------------------------------------")
print("overlap contribution:",info.overlap)
print("relative H contribution:",info.hamiltonian_relative)
print("dt-scaled T contribution:",info.time_connection_dt)
print("combined score:",info.score)
print("edge active:",bool(update.active_edges))
print("local omitted-score L2:",update.omitted_candidate_score_l2)
