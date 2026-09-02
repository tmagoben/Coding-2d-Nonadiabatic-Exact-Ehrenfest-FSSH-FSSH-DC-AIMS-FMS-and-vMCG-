import numpy as np
from gaussian_dynamics.analytic_molecular_backend_v19 import AnalyticMolecularLVCBackendV19, default_diatomic_two_mode_map_v19
from gaussian_dynamics.indexed_molecular_provider_v20 import IndexedTrackedMolecularDirectProviderV20
from gaussian_dynamics.electronic_operator_v21 import ElectronicOperatorProviderAdapterV21
from gaussian_dynamics.complex_gauge_v21 import PhaseMixingGaugeV21, GaugeTransformedOperatorProviderV21, random_unitary_v21

def _providers():
    g=default_diatomic_two_mode_map_v19(); base=ElectronicOperatorProviderAdapterV21(IndexedTrackedMolecularDirectProviderV20(AnalyticMolecularLVCBackendV19(g),g,rebuild_batch=4)); gauge=PhaseMixingGaugeV21(random_unitary_v21(2,2101),np.array([[.4,-.15],[-.25,.3]]),np.array([.2,-.4])); return base,GaugeTransformedOperatorProviderV21(base,gauge),gauge

def test_complex_operator_contract_and_force_are_gauge_invariant():
    base,gp,gauge=_providers(); q=np.array([-.31,.42]); a=base.evaluate(q); b=gp.evaluate(q); G=gauge.matrix(q)
    assert np.allclose(b.H,G.conj().T@a.H@G,atol=1e-12)
    c=np.array([.8+.1j,-.2+.5j]); c/=np.linalg.norm(c); cg=G.conj().T@c
    assert np.allclose(a.force_expectation(c),b.force_expectation(cg),atol=1e-11)
