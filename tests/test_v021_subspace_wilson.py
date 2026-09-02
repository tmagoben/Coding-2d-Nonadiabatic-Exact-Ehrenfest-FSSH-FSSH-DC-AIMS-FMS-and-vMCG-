import numpy as np
from gaussian_dynamics.complex_gauge_v21 import random_unitary_v21
from gaussian_dynamics.subspace_tracking_v21 import procrustes_subspace_alignment_v21
from gaussian_dynamics.wilson_loop_v21 import gauge_transform_cycle_links_v21, sorted_wilson_eigenphases_v21

def test_full_subspace_procrustes_handles_arbitrary_degenerate_rotation():
    out=procrustes_subspace_alignment_v21(random_unitary_v21(4,2111)); assert out.antihermitian_residual<1e-12; assert np.allclose(out.aligned_overlap,np.eye(4),atol=1e-12)

def test_wilson_spectrum_is_invariant_under_local_complex_gauges():
    links=[random_unitary_v21(3,2200+k) for k in range(5)]; gauges=[random_unitary_v21(3,2300+k) for k in range(5)]
    assert np.allclose(sorted_wilson_eigenphases_v21(links),sorted_wilson_eigenphases_v21(gauge_transform_cycle_links_v21(links,gauges)),atol=1e-11)
