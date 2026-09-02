import numpy as np
from gaussian_dynamics.analytic_molecular_backend_v19 import AnalyticMolecularLVCBackendV19, default_diatomic_two_mode_map_v19
from gaussian_dynamics.indexed_molecular_provider_v20 import IndexedTrackedMolecularDirectProviderV20
from gaussian_dynamics.electronic_operator_v21 import ElectronicOperatorProviderAdapterV21
from gaussian_dynamics.complex_gauge_v21 import PhaseMixingGaugeV21, GaugeTransformedOperatorProviderV21, random_unitary_v21
from gaussian_dynamics.block_sparse_molecular_v21 import BlockMolecularTBFV21, BlockSparseSettingsV21, build_dense_block_reference_v21
from gaussian_dynamics.block_dynamics_v21 import gauge_block_matrices_v21, gauge_covariance_errors_v21

def _p():
    g=default_diatomic_two_mode_map_v19(); return ElectronicOperatorProviderAdapterV21(IndexedTrackedMolecularDirectProviderV20(AnalyticMolecularLVCBackendV19(g),g,rebuild_batch=4))

def test_block_matrices_are_complex_u2_covariant_and_score_invariant():
    A=1.3*np.eye(2); basis=[BlockMolecularTBFV21(0,np.array([-.7,.35]),np.array([1.2,.1]),A),BlockMolecularTBFV21(1,np.array([0,.42]),np.array([-.4,.2]),A),BlockMolecularTBFV21(2,np.array([.75,.33]),np.array([.7,-.1]),A)]; qd=np.array([[.12,.03],[-.05,.08],[.09,-.04]]); pd=np.zeros_like(qd); st=BlockSparseSettingsV21(1e-14,1e-14,1e-14,local_omitted_score_l2_budget=0,use_kdtree=False); gauge=PhaseMixingGaugeV21(random_unitary_v21(2,2121),np.array([[.4,-.15],[-.25,.3]]),np.array([.2,-.4])); a=build_dense_block_reference_v21(basis,_p(),.01,qd,pd,st); b=build_dense_block_reference_v21(basis,GaugeTransformedOperatorProviderV21(_p(),gauge),.01,qd,pd,st); G,dG=gauge_block_matrices_v21(gauge,basis,qd); err=gauge_covariance_errors_v21(a,b,G,dG)
    assert max(err.values())<1e-11
    assert max(abs(a['pair_data'][e].score-b['pair_data'][e].score) for e in a['pair_data'])<1e-11
