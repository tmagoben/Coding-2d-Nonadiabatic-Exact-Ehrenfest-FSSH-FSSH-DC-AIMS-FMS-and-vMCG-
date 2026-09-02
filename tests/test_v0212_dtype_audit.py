from pathlib import Path
import numpy as np

from gaussian_dynamics.complex_dtype_audit_v212 import audit_pre_soc_complex_core_v212
from gaussian_dynamics.synthetic_operator_provider_v21 import SyntheticLinearOperatorConfigV21,SyntheticLinearOperatorProviderV21
from gaussian_dynamics.complex_gauge_v21 import PhaseMixingGaugeV21,GaugeTransformedOperatorProviderV21,random_unitary_v21


def test_pre_soc_complex_core_has_no_unclassified_real_casts_and_preserves_imaginary_data():
    package=Path(__file__).resolve().parents[1]/"gaussian_dynamics"
    audit=audit_pre_soc_complex_core_v212(package)
    assert audit.passed, audit.suspicious_casts

    base=SyntheticLinearOperatorProviderV21(SyntheticLinearOperatorConfigV21(nstate=3,nq=2,seed=21240))
    gauge=PhaseMixingGaugeV21(random_unitary_v21(3,21241),np.array([[.12,.03],[-.07,.08],[.05,-.09]]),np.array([.2,-.1,.3]))
    point=GaugeTransformedOperatorProviderV21(base,gauge).evaluate(np.array([.2,-.15]))
    assert np.max(np.abs(np.imag(point.H)))>1e-6
    assert np.max(np.abs(np.imag(point.dH_dq)))>1e-6
    assert np.max(np.abs(np.imag(point.connection_q)))>1e-6
