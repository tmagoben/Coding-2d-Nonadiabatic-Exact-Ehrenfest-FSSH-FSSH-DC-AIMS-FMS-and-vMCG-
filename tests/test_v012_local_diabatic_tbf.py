import numpy as np

from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.local_diabatic_tbf_v12 import (
    from_adiabatic_guided_tbf,
    parallel_transport_spinor_full_space,
    reset_to_instantaneous_adiabatic_spinor,
)


def test_full_space_parallel_transport_preserves_global_spinor():
    provider=AnalyticCI2DFrameProvider()
    raw=DynamicGraphTBF(
        uid=0,
        state=1,
        q=np.array([0.8,0.5]),
        p=np.array([0.3,0.1]),
        A=np.eye(2),
        node=("seed",0),
    )
    b=from_adiabatic_guided_tbf(raw,provider)
    initial=b.spinor.copy()

    old_q=b.q.copy()
    new_q=np.array([0.55,0.8])
    b.q=new_q.copy()
    parallel_transport_spinor_full_space(
        b,old_q,new_q,provider
    )

    # In a complete two-state electronic space, overlap transport exactly represents
    # the same physical vector in the new local frame.
    assert np.allclose(b.spinor,initial,atol=1e-12)


def test_instantaneous_reset_changes_physical_spinor_along_path():
    provider=AnalyticCI2DFrameProvider()
    raw=DynamicGraphTBF(
        uid=0,
        state=1,
        q=np.array([0.8,0.5]),
        p=np.array([0.3,0.1]),
        A=np.eye(2),
        node=("seed",0),
    )
    b=from_adiabatic_guided_tbf(raw,provider)
    initial=b.spinor.copy()

    b.q=np.array([0.55,0.8])
    reset_to_instantaneous_adiabatic_spinor(b,provider)

    assert not np.allclose(b.spinor,initial)
