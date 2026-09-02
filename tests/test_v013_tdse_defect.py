import numpy as np

from gaussian_dynamics.benchmark_campaign import CIPassageConfig
from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.born_huang_grid_v12 import build_born_huang_grid_2d
from gaussian_dynamics.spinor_complete_dynamics_v12 import (
    initialize_spinor_complete_coefficients,
)
from gaussian_dynamics.residual_basis_v13 import (
    GaussianCandidate,
    cartesian_offsets_2d,
    generate_gaussian_dictionary,
)
from gaussian_dynamics.tdse_defect_v13 import (
    compute_tdse_defect,
    defect_candidate_capture,
    enrich_basis_from_tdse_defect,
)


def _seed(config):
    return DynamicGraphTBF(
        uid=0,
        state=config.state,
        q=config.q_array(),
        p=config.p_array(),
        A=config.A_matrix(),
        node=("seed",0),
    )


def test_tdse_defect_is_finite_and_mostly_outside_projected_basis():
    c=CIPassageConfig()
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=c.mass)
    grid=build_born_huang_grid_2d(
        grid_n=36,
        half_width=c.half_width,
        mass=c.mass,
    )
    basis=[_seed(c)]
    C=initialize_spinor_complete_coefficients(
        basis,provider
    )

    defect=compute_tdse_defect(
        C,basis,provider,grid
    )

    assert np.isfinite(defect.residual_norm)
    assert defect.residual_norm > 0.0
    assert np.isfinite(defect.relative_to_hpsi)
    assert defect.projected_residual_norm < 0.05*defect.residual_norm


def test_defect_candidate_reports_positive_capture_fraction():
    c=CIPassageConfig()
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=c.mass)
    grid=build_born_huang_grid_2d(
        grid_n=36,
        half_width=c.half_width,
        mass=c.mass,
    )
    basis=[_seed(c)]
    C=initialize_spinor_complete_coefficients(
        basis,provider
    )
    defect=compute_tdse_defect(
        C,basis,provider,grid
    )

    candidate=GaussianCandidate(
        q=c.q_array()+np.array([-0.2,0.4]),
        p=c.p_array(),
        A=2.0*c.A_matrix(),
        state=c.state,
        label="test",
    )

    score=defect_candidate_capture(
        candidate,
        defect,
        basis,
        grid,
        condition_limit=1e8,
    )

    assert score is not None
    assert score.captured_defect_norm > 0.0
    assert 0.0 < score.capture_fraction <= 1.0+1e-8
    assert score.expanded_condition_number > 1.0


def test_defect_enrichment_preserves_wavefunction_and_reduces_defect():
    c=CIPassageConfig()
    provider=AnalyticCI2DFrameProvider(nuclear_mass_au=c.mass)
    grid=build_born_huang_grid_2d(
        grid_n=32,
        half_width=c.half_width,
        mass=c.mass,
    )
    basis=[_seed(c)]
    C=initialize_spinor_complete_coefficients(
        basis,provider
    )

    candidates=generate_gaussian_dictionary(
        c.q_array(),
        c.p_array(),
        c.A_matrix(),
        c.state,
        cartesian_offsets_2d(radius=0.6,spacing=0.3),
        width_scales=(1.0,2.0,4.0),
    )

    result=enrich_basis_from_tdse_defect(
        C,basis,provider,grid,candidates,
        condition_limit=1e7,
    )

    assert result is not None
    assert len(result.basis)==2
    assert np.allclose(result.coefficients[-2:],0.0)

    # Zero coefficient insertion leaves Psi itself unchanged.
    assert np.allclose(
        result.defect_before.wavefunction,
        result.defect_after.wavefunction,
        atol=2e-12,
    )

    assert result.defect_after.residual_norm < result.defect_before.residual_norm
    assert result.actual_squared_defect_reduction > 0.0

    predicted=result.score.captured_defect_norm**2
    assert np.isclose(
        result.actual_squared_defect_reduction,
        predicted,
        rtol=5e-3,
        atol=5e-6,
    )
