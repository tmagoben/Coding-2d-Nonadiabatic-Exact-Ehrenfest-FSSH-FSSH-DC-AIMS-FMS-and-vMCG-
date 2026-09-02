import numpy as np

from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.born_huang_grid_v12 import build_born_huang_grid_2d
from gaussian_dynamics.tdse_defect_v13 import compute_tdse_defect
from gaussian_dynamics.defect_candidates_v15 import (
    generate_energy_conserving_defect_candidates_v15,
    rank_dynamic_defect_candidates_cached,
)
from gaussian_dynamics.defect_candidates_v18 import (
    rank_dynamic_defect_candidates_batched_v18,
)
from gaussian_dynamics.pair_cache_v15 import GaussianPairCache
from gaussian_dynamics.residual_basis_v13 import nuclear_overlap_matrix


def _basis():
    return [
        DynamicGraphTBF(
            uid=0,state=1,
            q=np.array([-0.45,0.3]),
            p=np.array([0.8,0.1]),
            A=1.2*np.eye(2),
            node=("a",0),
        ),
        DynamicGraphTBF(
            uid=1,state=0,
            q=np.array([0.35,0.25]),
            p=np.array([-0.2,0.2]),
            A=0.9*np.eye(2),
            node=("b",1),
        ),
    ]


def test_batched_residual_ranking_matches_dense_v15():
    provider=AnalyticCI2DFrameProvider(
        nuclear_mass_au=20.0
    )
    basis=_basis()
    grid=build_born_huang_grid_2d(
        grid_n=20,half_width=3.0,mass=20.0
    )

    C=np.zeros(2*len(basis),complex)
    C[1]=1.0
    defect=compute_tdse_defect(
        C,basis,provider,grid
    )

    candidates=generate_energy_conserving_defect_candidates_v15(
        basis,provider,
        position_shifts=(0.0,0.05,-0.05),
        width_scales=(0.8,1.0,1.2),
    )
    cache=GaussianPairCache(basis)
    Snuc=nuclear_overlap_matrix(basis)

    dense=rank_dynamic_defect_candidates_cached(
        defect,basis,candidates,grid,cache,Snuc,
        exact_condition_top=8,max_return=6,
        condition_limit=1e10,
    )
    batched,diag=rank_dynamic_defect_candidates_batched_v18(
        defect,basis,candidates,grid,cache,Snuc,
        exact_condition_top=8,max_return=6,
        condition_limit=1e10,
        batch_size=3,
        return_diagnostics=True,
    )

    assert [x.candidate_index for x in batched] == [
        x.candidate_index for x in dense
    ]
    assert np.allclose(
        [x.capture_fraction for x in batched],
        [x.capture_fraction for x in dense],
        rtol=2e-12,atol=2e-12,
    )
    assert (
        diag.peak_candidate_grid_elements
        <=3*diag.grid_points
    )
    assert diag.peak_grid_element_reduction_fraction>0.0
