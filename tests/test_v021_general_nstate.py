import numpy as np
from gaussian_dynamics.synthetic_operator_provider_v21 import SyntheticLinearOperatorProviderV21, SyntheticLinearOperatorConfigV21
from gaussian_dynamics.block_sparse_molecular_v21 import BlockMolecularTBFV21, BlockSparseSettingsV21, build_dense_block_reference_v21

def test_block_framework_has_no_two_state_hardcoding():
    p=SyntheticLinearOperatorProviderV21(SyntheticLinearOperatorConfigV21(nstate=5,nq=2,seed=515)); A=1.1*np.eye(2); basis=[BlockMolecularTBFV21(i,np.array([-.5+.5*i,.2]),np.array([.1*(-1)**i,0]),A) for i in range(3)]; qd=np.array([[.1,0],[-.05,.02],[.03,-.04]]); d=build_dense_block_reference_v21(basis,p,.01,qd,np.zeros_like(qd),BlockSparseSettingsV21(1e-14,1e-14,1e-14,local_omitted_score_l2_budget=0,use_kdtree=False)); assert d['S'].shape==(15,15); assert np.allclose(d['H'],d['H'].conj().T,atol=1e-12)
