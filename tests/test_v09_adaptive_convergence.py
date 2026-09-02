import numpy as np

from gaussian_dynamics.adaptive_spawning import CouplingExposureTracker
from gaussian_dynamics.convergence import observed_order, scalar_refinement_study


def test_integrated_coupling_trigger_is_timestep_consistent_for_constant_rate():
    def trigger_time(dt):
        tr=CouplingExposureTracker(action_threshold=0.1,coupling_floor=0.0)
        t=0.0
        for _ in range(1000):
            t+=dt
            ready,_=tr.update(('p',1),2.0,dt)
            if ready:
                return t
        raise AssertionError('did not trigger')

    t1=trigger_time(0.01)
    t2=trigger_time(0.005)
    assert abs(t1-t2) <= 0.01
    assert abs(t1-0.05) <= 0.01


def test_observed_order_recovers_second_order_sequence():
    h=np.array([0.4,0.2,0.1,0.05])
    values=1.23+2.0*h*h
    study=scalar_refinement_study(h,values)
    assert np.allclose(study.observed_orders,2.0,atol=1e-10)
    assert abs(observed_order(0.04,0.01)-2.0) < 1e-12
