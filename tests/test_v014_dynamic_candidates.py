import numpy as np

from gaussian_dynamics.benchmark_campaign import CIPassageConfig
from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.born_huang_grid_v12 import build_born_huang_grid_2d
from gaussian_dynamics.spinor_complete_dynamics_v12 import (
    initialize_spinor_complete_coefficients,
)
from gaussian_dynamics.spinor_complete_lvc_v12 import (
    build_spinor_complete_lvc_matrices,
)
from gaussian_dynamics.adaptive_defect_dynamics_v14 import (
    compute_tdse_defect_with_matrices,
)
from gaussian_dynamics.defect_candidates_v14 import (
    generate_energy_conserving_defect_candidates,
    rank_dynamic_defect_candidates_prepared,
)
from gaussian_dynamics.tdse_defect_v13 import (
    defect_candidate_capture,
)


def _seed(c):
    return DynamicGraphTBF(
        uid=0,state=c.state,
        q=c.q_array(),p=c.p_array(),A=c.A_matrix(),
        node=("seed",0),
    )


def test_dynamic_candidates_are_energy_conserving_and_residual_rank_is_consistent():
    c=CIPassageConfig()
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=c.mass
    )
    grid=build_born_huang_grid_2d(
        grid_n=28,
        half_width=c.half_width,
        mass=c.mass,
    )
    basis=[_seed(c)]
    C=initialize_spinor_complete_coefficients(
        basis,provider
    )
    S,H,Snuc=build_spinor_complete_lvc_matrices(
        basis,provider
    )
    defect=compute_tdse_defect_with_matrices(
        C,basis,provider,grid,S,H,Snuc
    )

    dynamic=generate_energy_conserving_defect_candidates(
        basis,provider,
        position_shifts=(0.0,0.05,-0.05),
        width_scales=(0.8,1.2),
    )
    assert len(dynamic)>0
    assert max(abs(x.energy_residual) for x in dynamic)<1e-9

    ranked=rank_dynamic_defect_candidates_prepared(
        defect,basis,dynamic,grid,
        condition_limit=1e7,
        max_return=3,
    )
    assert ranked
    best=ranked[0]
    item=dynamic[best.candidate_index]

    slow=defect_candidate_capture(
        item.candidate,
        defect,
        basis,
        grid,
        condition_limit=1e7,
    )
    assert slow is not None
    assert np.isclose(
        best.capture_fraction,
        slow.capture_fraction,
        rtol=3e-3,
        atol=3e-5,
    )
