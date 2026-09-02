import numpy as np

from gaussian_dynamics.synthetic_operator_provider_v21 import SyntheticLinearOperatorConfigV21,SyntheticLinearOperatorProviderV21
from gaussian_dynamics.complex_gauge_v21 import PhaseMixingGaugeV21,GaugeTransformedOperatorProviderV21,random_unitary_v21
from gaussian_dynamics.block_sparse_molecular_v21 import BlockMolecularTBFV21,BlockSparseSettingsV21,build_dense_block_reference_v21
from gaussian_dynamics.block_dynamics_v21 import gauge_block_matrices_v21
from gaussian_dynamics.electronic_observables_v212 import ElectronicObservableV212,build_electronic_observable_matrix_v212,observable_expectation_v212


def test_generic_electronic_observable_is_complex_gauge_invariant():
    def base(): return SyntheticLinearOperatorProviderV21(SyntheticLinearOperatorConfigV21(nstate=3,nq=2,seed=21220,mass=20.0))
    gauge=PhaseMixingGaugeV21(random_unitary_v21(3,21221),np.array([[.11,.02],[-.05,.07],[.04,-.08]]),np.array([.2,-.3,.1]))
    A=np.eye(2)
    basis=[
        BlockMolecularTBFV21(0,np.array([-.45,.2]),np.array([.1,0]),A),
        BlockMolecularTBFV21(1,np.array([.55,.25]),np.array([-.08,.03]),1.3*A),
    ]
    qdot=np.zeros((2,2)); pdot=np.zeros_like(qdot)
    settings=BlockSparseSettingsV21(enter_score=1e-14,exit_score=1e-14,search_overlap_floor=1e-14,local_omitted_score_l2_budget=0,use_kdtree=False)
    pb=base(); pg=GaugeTransformedOperatorProviderV21(base(),gauge)
    db=build_dense_block_reference_v21(basis,pb,.01,qdot,pdot,settings)
    dg=build_dense_block_reference_v21(basis,pg,.01,qdot,pdot,settings)
    observable=ElectronicObservableV212("physical dH/dq0",lambda snap:snap.point.dH_dq[0])
    Ob=build_electronic_observable_matrix_v212(basis,pb,observable)
    Og=build_electronic_observable_matrix_v212(basis,pg,observable)
    C=np.array([.6+.1j,.2-.2j,-.1+.3j,.3+.2j,.15-.1j,.25+.05j])
    G0,_=gauge_block_matrices_v21(gauge,basis,qdot)
    Cg=G0.conj().T@C
    eb=observable_expectation_v212(C,db["S"],Ob)
    eg=observable_expectation_v212(Cg,dg["S"],Og)
    assert abs(eb-eg)<2e-11
    assert np.max(np.abs(np.imag(pg.evaluate(basis[0].q).H)))>1e-6
