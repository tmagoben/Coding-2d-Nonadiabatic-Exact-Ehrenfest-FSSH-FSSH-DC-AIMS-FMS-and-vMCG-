import numpy as np

from gaussian_dynamics.complexity_v14 import (
    asymptotic_complexity,
    dense_dimension_cost_proxy,
    candidate_ranking_cost_proxy,
)
from gaussian_dynamics.residual_pruning_v14 import (
    leave_one_out_projection_losses,
)
from gaussian_dynamics.paired_basis_management_v12 import (
    project_spinor_coefficients_to_subset,
    spinor_wavefunction_norm,
)


def test_complexity_model_exposes_dominant_scalings():
    c=asymptotic_complexity()
    assert "(sN)^3" in c.coefficient_solve
    assert "K G" in c.prepared_candidate_ranking
    assert "N^3" in c.pruning

    # Dense solve proxy must scale cubically in basis size.
    assert dense_dimension_cost_proxy(8)==8*dense_dimension_cost_proxy(4)

    a=candidate_ranking_cost_proxy(
        n_basis=4,n_candidates=10,grid_points=100
    )
    b=candidate_ranking_cost_proxy(
        n_basis=4,n_candidates=20,grid_points=100
    )
    assert b>a


def test_leave_one_out_formula_matches_direct_projection_loss():
    S=np.array([
        [1.0,0.3+0.05j,0.1],
        [0.3-0.05j,1.0,0.25],
        [0.1,0.25,1.0],
    ],dtype=complex)
    C=np.array([
        [0.5+0.1j,0.2],
        [0.15-0.05j,-0.1j],
        [0.08,0.03+0.02j],
    ],dtype=complex)

    scores=leave_one_out_projection_losses(
        C,S,uids=[10,11,12]
    )
    by_index={s.index:s for s in scores}
    old=spinor_wavefunction_norm(C,S)

    for j in range(3):
        keep=np.array([i for i in range(3) if i!=j])
        _,loss=project_spinor_coefficients_to_subset(
            C,S,keep
        )
        assert np.isclose(
            by_index[j].absolute_projection_loss,
            loss,
            rtol=1e-11,
            atol=1e-12,
        )
        assert np.isclose(
            by_index[j].fractional_projection_loss,
            loss/old,
            rtol=1e-11,
            atol=1e-12,
        )
