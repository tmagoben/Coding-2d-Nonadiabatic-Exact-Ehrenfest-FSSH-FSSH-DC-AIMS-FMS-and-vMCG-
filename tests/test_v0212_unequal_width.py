import numpy as np

from gaussian_dynamics.synthetic_operator_provider_v21 import (
    SyntheticLinearOperatorConfigV21,SyntheticLinearOperatorProviderV21,
)
from gaussian_dynamics.complex_gauge_v21 import (
    PhaseMixingGaugeV21,GaugeTransformedOperatorProviderV21,random_unitary_v21,
)
from gaussian_dynamics.block_sparse_molecular_v21 import (
    BlockMolecularTBFV21,BlockSparseSettingsV21,build_dense_block_reference_v21,
)
from gaussian_dynamics.block_dynamics_v21 import gauge_block_matrices_v21,gauge_covariance_errors_v21


def _provider():
    return SyntheticLinearOperatorProviderV21(
        SyntheticLinearOperatorConfigV21(nstate=3,nq=2,seed=21201,mass=25.0)
    )


def test_unequal_width_complex_block_matrices_remain_gauge_covariant():
    basis=[
        BlockMolecularTBFV21(0,np.array([-0.6,0.2]),np.array([0.2,0.05]),np.array([[1.2,0.1],[0.1,1.7]])),
        BlockMolecularTBFV21(1,np.array([0.1,0.35]),np.array([-0.1,0.08]),np.array([[1.8,-0.12],[-0.12,1.1]])),
        BlockMolecularTBFV21(2,np.array([0.75,0.15]),np.array([0.15,-0.04]),np.diag([0.9,2.0])),
    ]
    qdot=np.array([[0.08,0.02],[-0.03,0.05],[0.06,-0.02]])
    pdot=np.zeros_like(qdot)
    gauge=PhaseMixingGaugeV21(
        random_unitary_v21(3,21202),
        np.array([[0.1,-0.04],[-0.06,0.08],[0.03,0.05]]),
        np.array([0.2,-0.1,0.3]),
    )
    settings=BlockSparseSettingsV21(
        enter_score=1e-14,exit_score=1e-14,search_overlap_floor=1e-14,
        local_omitted_score_l2_budget=0.0,use_kdtree=False,
    )
    base=_provider()
    transformed=GaugeTransformedOperatorProviderV21(_provider(),gauge)
    A=build_dense_block_reference_v21(basis,base,0.01,qdot,pdot,settings)
    B=build_dense_block_reference_v21(basis,transformed,0.01,qdot,pdot,settings)
    G,dG=gauge_block_matrices_v21(gauge,basis,qdot)
    err=gauge_covariance_errors_v21(A,B,G,dG)
    assert err["S_relative_error"]<2e-12
    assert err["H_relative_error"]<2e-12
    assert err["T_relative_error"]<2e-12
    assert all(np.isfinite(p.score) for p in A["pair_data"].values())
