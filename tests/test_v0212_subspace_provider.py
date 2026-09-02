import numpy as np

from gaussian_dynamics.synthetic_operator_provider_v21 import SyntheticLinearOperatorConfigV21,SyntheticLinearOperatorProviderV21
from gaussian_dynamics.complex_gauge_v21 import PhaseMixingGaugeV21,GaugeTransformedOperatorProviderV21,random_unitary_v21
from gaussian_dynamics.subspace_provider_v212 import SubspaceAwareOperatorProviderV212,SubspaceTrackingSettingsV212


def test_subspace_provider_checks_full_manifold_without_forcing_root_gauge():
    base=SyntheticLinearOperatorProviderV21(SyntheticLinearOperatorConfigV21(nstate=5,nq=2,seed=21230))
    gauge=PhaseMixingGaugeV21(random_unitary_v21(5,21231),np.column_stack([np.linspace(.02,.08,5),np.linspace(-.05,.04,5)]),np.linspace(-.2,.2,5))
    raw=GaugeTransformedOperatorProviderV21(base,gauge)
    p=SubspaceAwareOperatorProviderV212(raw,dimension=2,settings=SubspaceTrackingSettingsV212(minimum_singular_value=.999999,ambiguity_policy="raise",rebuild_batch=3))
    points=[np.array([x,.2+.01*np.sin(3*x)]) for x in np.linspace(-.5,.5,9)]
    snaps=[p.evaluate_snapshot(q) for q in points]
    diag=p.diagnostics_dict()["subspace"]
    assert diag["subspace_checks"]==8
    assert diag["subspace_ambiguities"]==0
    assert diag["minimum_seen_singular_value"]>1-1e-12
    # The wrapper returns the provider's frame unchanged; it only diagnoses/alignment-maps the subspace.
    assert np.allclose(snaps[-1].state_vectors,raw.evaluate_snapshot(points[-1]).state_vectors)
