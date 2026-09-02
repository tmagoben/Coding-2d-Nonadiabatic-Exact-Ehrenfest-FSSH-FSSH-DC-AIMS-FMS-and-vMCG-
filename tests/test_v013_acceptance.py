from gaussian_dynamics.v13_benchmark import (
    V13AcceptanceThresholds,
    evaluate_v13_acceptance,
)


def _good_reference():
    return {
        "initial_density_error":0.02,
        "projected_dynamics_density_error":1e-4,
        "target_density_error":0.025,
        "target_population_error":0.02,
        "coherence_phase_error":0.001,
        "max_norm_drift":1e-6,
        "max_condition_number":1000.0,
    }


def _good_defect():
    return {
        "predicted_squared_reduction":0.01,
        "actual_squared_reduction":0.010001,
    }


def _history():
    return [
        {"relative_residual":0.8},
        {"relative_residual":0.5},
        {"relative_residual":0.3},
    ]


def test_v13_acceptance_requires_monotone_residual_and_defect_prediction():
    result=evaluate_v13_acceptance(
        _good_reference(),
        _good_defect(),
        _history(),
        V13AcceptanceThresholds(),
    )
    assert result["passed"]

    bad_history=[
        {"relative_residual":0.8},
        {"relative_residual":0.82},
    ]
    bad=evaluate_v13_acceptance(
        _good_reference(),
        _good_defect(),
        bad_history,
        V13AcceptanceThresholds(),
    )
    assert not bad["passed"]
    assert not bad["checks"]["monotone_residual_refinement"]


def test_v13_acceptance_detects_wrong_defect_gain_prediction():
    defect=_good_defect()
    defect["actual_squared_reduction"]=0.02

    result=evaluate_v13_acceptance(
        _good_reference(),
        defect,
        _history(),
        V13AcceptanceThresholds(),
    )
    assert not result["passed"]
    assert not result["checks"]["defect_gain_prediction"]
