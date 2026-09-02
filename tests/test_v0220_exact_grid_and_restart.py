import numpy as np
import pytest

from gaussian_dynamics.analytic_soc_models_v220 import (
    AnalyticDoubletSOCProviderV220,
    AnalyticSingletTripletSOCProviderV220,
    DoubletSOCConfigV220,
    SingletTripletSOCConfigV220,
)
from gaussian_dynamics.block_sparse_molecular_v21 import BlockMolecularTBFV21
from gaussian_dynamics.complex_gauge_v21 import (
    GaugeTransformedOperatorProviderV21,
    PhaseMixingGaugeV21,
    random_unitary_v21,
)
from gaussian_dynamics.electronic_contract_v213 import (
    ContractedElectronicOperatorProviderV213,
)
from gaussian_dynamics.checkpoint_restart_v214 import (
    SelfConsistentBlockSettingsV214,
    run_self_consistent_block_dynamics_v214,
)
from gaussian_dynamics.spinor_exact_grid_v220 import (
    SpinorGridSettingsV220,
    initial_gaussian_spinor_v220,
    phase_aligned_spinor_grid_error_v220,
    run_spinor_exact_grid_v220,
    spinor_split_operator_step_v220,
)


def _grid():
    return np.linspace(-8.0, 8.0, 256, endpoint=False)


def _initial_grid_state(x):
    return initial_gaussian_spinor_v220(
        x, np.asarray([1.0, 0.0, 0.0, 0.0]), center=-1.0, momentum=1.2, width=0.7
    )


@pytest.mark.parametrize(
    "config_type,provider_type,population_names",
    [
        (
            SingletTripletSOCConfigV220,
            AnalyticSingletTripletSOCProviderV220,
            ("singlet", "triplet"),
        ),
        (
            DoubletSOCConfigV220,
            AnalyticDoubletSOCProviderV220,
            ("doublet_1", "doublet_2"),
        ),
    ],
)
def test_enabled_zero_soc_grid_dynamics_exactly_matches_disabled_spin_free(
    config_type, provider_type, population_names
):
    x = _grid()
    psi0 = _initial_grid_state(x)
    settings = SpinorGridSettingsV220(dt=0.05, steps=30, store_every=10)
    enabled = run_spinor_exact_grid_v220(
        provider_type(config_type(soc_scale=0.0, soc_enabled=True)),
        x,
        psi0,
        settings=settings,
    )
    disabled = run_spinor_exact_grid_v220(
        provider_type(config_type(soc_scale=0.0, soc_enabled=False)),
        x,
        psi0,
        settings=settings,
    )

    assert np.array_equal(enabled["psi"], disabled["psi"])
    assert np.array_equal(enabled["norm"], disabled["norm"])
    for name in population_names:
        assert np.array_equal(enabled["populations"][name], disabled["populations"][name])
    assert enabled["populations"][population_names[0]][-1] == pytest.approx(1.0, abs=2.0e-14)


@pytest.mark.parametrize(
    "provider_type,transfer_population",
    [
        (AnalyticSingletTripletSOCProviderV220, "triplet"),
        (AnalyticDoubletSOCProviderV220, "doublet_2"),
    ],
)
def test_soc_exact_grid_conserves_norm_energy_and_transfers_population(
    provider_type, transfer_population
):
    x = _grid()
    output = run_spinor_exact_grid_v220(
        provider_type(),
        x,
        _initial_grid_state(x),
        settings=SpinorGridSettingsV220(dt=0.04, steps=100, store_every=20),
    )

    assert output["maximum_norm_drift"] < 5.0e-14
    assert output["maximum_energy_drift"] < 1.0e-12
    assert output["populations"][transfer_population][-1] > 5.0e-5
    population_sum = sum(output["populations"].values())
    assert np.max(np.abs(population_sum - 1.0)) < 2.0e-14


def test_exact_grid_strang_propagation_has_second_order_timestep_convergence():
    x = _grid()
    dx = x[1] - x[0]
    provider = AnalyticDoubletSOCProviderV220()
    psi0 = _initial_grid_state(x)
    outputs = []
    for dt in (0.08, 0.04, 0.02):
        steps = int(round(4.0 / dt))
        outputs.append(
            run_spinor_exact_grid_v220(
                provider,
                x,
                psi0,
                settings=SpinorGridSettingsV220(
                    dt=dt, steps=steps, store_every=steps
                ),
            )
        )
    coarse_difference = phase_aligned_spinor_grid_error_v220(
        outputs[0]["psi"][-1], outputs[1]["psi"][-1], dx
    )
    fine_difference = phase_aligned_spinor_grid_error_v220(
        outputs[1]["psi"][-1], outputs[2]["psi"][-1], dx
    )
    observed_order = np.log(coarse_difference / fine_difference) / np.log(2.0)

    assert coarse_difference > fine_difference
    assert observed_order > 1.99


def test_exact_grid_step_is_forward_backward_reversible_to_roundoff():
    x = _grid()
    dx = x[1] - x[0]
    provider = AnalyticSingletTripletSOCProviderV220()
    psi0 = _initial_grid_state(x)
    potential = np.asarray(
        [provider.evaluate_snapshot(np.asarray([coordinate])).point.H for coordinate in x]
    )
    psi = psi0.copy()
    for _ in range(20):
        psi = spinor_split_operator_step_v220(
            psi, dx, 0.05, provider.config.mass, potential
        )
    for _ in range(20):
        psi = spinor_split_operator_step_v220(
            psi, dx, -0.05, provider.config.mass, potential
        )

    assert phase_aligned_spinor_grid_error_v220(psi0, psi, dx) < 2.0e-13


def _basis():
    return [
        BlockMolecularTBFV21(
            3, np.asarray([-0.55]), np.asarray([0.20]), np.asarray([[1.15]])
        ),
        BlockMolecularTBFV21(
            8, np.asarray([0.45]), np.asarray([-0.10]), np.asarray([[1.35]])
        ),
    ]


def _coefficients():
    return np.asarray(
        [0.76 + 0.08j, 0.0, 0.0, 0.0, 0.22 - 0.11j, 0.0, 0.0, 0.0]
    )


def _dense_settings():
    return SelfConsistentBlockSettingsV214(
        use_dense_reference=True,
        corrector_iterations=3,
        momentum_tolerance=1.0e-12,
    )


def test_soc_active_gaussian_checkpoint_restart_matches_uninterrupted():
    provider = AnalyticSingletTripletSOCProviderV220()
    full = run_self_consistent_block_dynamics_v214(
        provider,
        provider.provenance,
        initial_basis=_basis(),
        C0=_coefficients(),
        dt=0.002,
        steps=8,
        store_every=2,
        settings=_dense_settings(),
    )
    first_provider = AnalyticSingletTripletSOCProviderV220()
    first = run_self_consistent_block_dynamics_v214(
        first_provider,
        first_provider.provenance,
        initial_basis=_basis(),
        C0=_coefficients(),
        dt=0.002,
        steps=3,
        store_every=1,
        settings=_dense_settings(),
    )
    resumed_provider = AnalyticSingletTripletSOCProviderV220()
    resumed = run_self_consistent_block_dynamics_v214(
        resumed_provider,
        resumed_provider.provenance,
        checkpoint=first["checkpoint"],
        steps=5,
        store_every=1,
        settings=_dense_settings(),
    )

    assert max(
        np.linalg.norm(left.q - right.q)
        for left, right in zip(full["final_basis"], resumed["final_basis"])
    ) == 0.0
    assert max(
        np.linalg.norm(left.p - right.p)
        for left, right in zip(full["final_basis"], resumed["final_basis"])
    ) == 0.0
    assert np.array_equal(full["final_coefficients"], resumed["final_coefficients"])
    assert resumed["checkpoint"].step == 8


def test_soc_checkpoint_rejects_changed_soc_parameters():
    provider = AnalyticDoubletSOCProviderV220()
    first = run_self_consistent_block_dynamics_v214(
        provider,
        provider.provenance,
        initial_basis=_basis(),
        C0=_coefficients(),
        dt=0.002,
        steps=1,
        store_every=1,
        settings=_dense_settings(),
    )
    changed = AnalyticDoubletSOCProviderV220(DoubletSOCConfigV220(soc_scale=0.0026))

    with pytest.raises(ValueError, match="provider provenance fingerprint mismatch"):
        run_self_consistent_block_dynamics_v214(
            changed,
            changed.provenance,
            checkpoint=first["checkpoint"],
            steps=1,
            store_every=1,
            settings=_dense_settings(),
        )


@pytest.mark.parametrize(
    "provider_type,transferred_population",
    [
        (AnalyticSingletTripletSOCProviderV220, "triplet"),
        (AnalyticDoubletSOCProviderV220, "doublet_2"),
    ],
)
def test_short_time_gaussian_soc_population_agrees_with_independent_exact_grid(
    provider_type, transferred_population
):
    provider = provider_type()
    basis = [
        BlockMolecularTBFV21(
            1, np.asarray([-1.0]), np.asarray([1.2]), np.asarray([[0.7]])
        )
    ]
    gaussian = run_self_consistent_block_dynamics_v214(
        provider,
        provider.provenance,
        initial_basis=basis,
        C0=np.asarray([1.0, 0.0, 0.0, 0.0], dtype=complex),
        dt=0.002,
        steps=100,
        store_every=100,
        settings=SelfConsistentBlockSettingsV214(
            use_dense_reference=True,
            corrector_iterations=4,
            momentum_tolerance=1.0e-12,
        ),
    )
    coefficients = gaussian["final_coefficients"]
    projector = provider.projectors[transferred_population]
    gaussian_population = float(
        np.real(np.vdot(coefficients, projector @ coefficients))
        / np.real(np.vdot(coefficients, coefficients))
    )

    x = np.linspace(-8.0, 8.0, 512, endpoint=False)
    exact = run_spinor_exact_grid_v220(
        provider_type(),
        x,
        _initial_grid_state(x),
        settings=SpinorGridSettingsV220(dt=0.01, steps=20, store_every=20),
    )
    exact_population = float(exact["populations"][transferred_population][-1])

    assert exact_population > 1.0e-7
    assert abs(gaussian_population - exact_population) < 1.0e-8


def test_soc_active_dynamics_is_covariant_in_a_moving_complex_frame():
    config = DoubletSOCConfigV220()
    base = AnalyticDoubletSOCProviderV220(config)
    gauge = PhaseMixingGaugeV21(
        random_unitary_v21(4, 22031),
        np.asarray([[0.11], [-0.08], [0.17], [-0.13]]),
        np.asarray([0.20, -0.31, 0.14, -0.09]),
    )
    local_provenance = config.provenance("local_general")
    local = ContractedElectronicOperatorProviderV213(
        GaugeTransformedOperatorProviderV21(
            AnalyticDoubletSOCProviderV220(config), gauge
        ),
        local_provenance,
    )
    basis = _basis()
    coefficients = _coefficients()
    local_coefficients = np.concatenate(
        [
            gauge.matrix(item.q).conj().T
            @ coefficients[4 * index : 4 * index + 4]
            for index, item in enumerate(basis)
        ]
    )
    common = dict(
        dt=0.002,
        steps=6,
        store_every=2,
        settings=_dense_settings(),
    )
    fixed = run_self_consistent_block_dynamics_v214(
        base,
        base.provenance,
        initial_basis=basis,
        C0=coefficients,
        **common,
    )
    moving = run_self_consistent_block_dynamics_v214(
        local,
        local_provenance,
        initial_basis=basis,
        C0=local_coefficients,
        **common,
    )
    mapped = np.concatenate(
        [
            gauge.matrix(item.q)
            @ moving["final_coefficients"][4 * index : 4 * index + 4]
            for index, item in enumerate(moving["final_basis"])
        ]
    )
    metric = fixed["final_S"]
    overlap = np.vdot(fixed["final_coefficients"], metric @ mapped)
    phase = np.exp(-1j * np.angle(overlap))
    difference = phase * mapped - fixed["final_coefficients"]
    coefficient_error = float(
        np.sqrt(max(np.real(np.vdot(difference, metric @ difference)), 0.0))
    )

    assert max(
        np.linalg.norm(left.q - right.q)
        for left, right in zip(fixed["final_basis"], moving["final_basis"])
    ) < 1.0e-13
    assert max(
        np.linalg.norm(left.p - right.p)
        for left, right in zip(fixed["final_basis"], moving["final_basis"])
    ) < 1.0e-13
    assert coefficient_error < 2.0e-12
