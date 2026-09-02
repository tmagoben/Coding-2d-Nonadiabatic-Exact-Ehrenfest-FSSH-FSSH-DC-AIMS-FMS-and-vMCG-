from dataclasses import replace

import numpy as np
import pytest

from gaussian_dynamics.analytic_soc_models_v220 import (
    AnalyticDoubletSOCProviderV220,
    AnalyticSingletTripletSOCProviderV220,
    DoubletSOCConfigV220,
    SOCOperatorComponentsV220,
    SingletTripletSOCConfigV220,
)
from gaussian_dynamics.complex_gauge_v21 import (
    GaugeTransformedOperatorProviderV21,
    PhaseMixingGaugeV21,
    random_unitary_v21,
)
from gaussian_dynamics.electronic_contract_v213 import (
    ContractedElectronicOperatorProviderV213,
    ElectronicModelSpaceV213,
    ElectronicStateDescriptorV213,
)
from gaussian_dynamics.physical_soc_validation_v220 import (
    audit_kramers_degeneracy_v220,
    audit_physical_soc_provider_v220,
    projector_population_v220,
    time_reversal_residual_v220,
    transform_projector_v220,
    transform_time_reversal_matrix_v220,
)
from gaussian_dynamics.provider_differential_audit_v214 import (
    require_provider_differential_contract_v214,
)


def _gauge():
    return PhaseMixingGaugeV21(
        random_unitary_v21(4, 22031),
        np.asarray(
            [
                [0.11],
                [-0.08],
                [0.17],
                [-0.13],
            ]
        ),
        np.asarray([0.20, -0.31, 0.14, -0.09]),
    )


def test_singlet_triplet_physical_soc_contract_passes():
    provider = AnalyticSingletTripletSOCProviderV220()
    report = audit_physical_soc_provider_v220(
        provider, np.asarray([0.17]), fermionic=False
    )

    assert report.passed
    assert all(report.checks.values())
    assert report.H_composition_error == 0.0
    assert report.K_composition_error == 0.0
    assert report.maximum_time_reversal_residual == 0.0
    assert report.time_reversal_square_residual == 0.0
    assert report.soc_force_error < 1.0e-12


def test_doublet_physical_soc_contract_and_kramers_degeneracy_pass():
    provider = AnalyticDoubletSOCProviderV220()
    report = audit_physical_soc_provider_v220(
        provider, np.asarray([-0.23]), fermionic=True
    )
    kramers = audit_kramers_degeneracy_v220(
        provider,
        [np.asarray([-1.2]), np.asarray([0.0]), np.asarray([1.1])],
    )

    assert report.passed
    assert kramers.passed
    assert kramers.time_reversal_square_residual == 0.0
    assert kramers.maximum_pair_splitting < 2.0e-17


def test_constant_soc_has_nonzero_H_soc_and_exactly_zero_K_soc():
    config = SingletTripletSOCConfigV220(
        lambda_real_gradient=0.0,
        lambda_imag_gradient=0.0,
        lambda_zero_gradient=0.0,
    )
    provider = AnalyticSingletTripletSOCProviderV220(config)
    components = provider.components(np.asarray([0.31]))
    report = audit_physical_soc_provider_v220(
        provider, np.asarray([0.31]), fermionic=False
    )

    assert np.linalg.norm(components.H_soc) > 0.0
    assert np.linalg.norm(components.K_soc) == 0.0
    assert report.passed
    assert report.soc_force_analytic == 0.0
    assert report.soc_force_finite_difference == pytest.approx(0.0, abs=1.0e-14)


@pytest.mark.parametrize(
    "config_type,provider_type,fermionic",
    [
        (SingletTripletSOCConfigV220, AnalyticSingletTripletSOCProviderV220, False),
        (DoubletSOCConfigV220, AnalyticDoubletSOCProviderV220, True),
    ],
)
def test_enabled_zero_soc_and_disabled_spin_free_operators_are_exactly_equal(
    config_type, provider_type, fermionic
):
    enabled = provider_type(config_type(soc_scale=0.0, soc_enabled=True))
    disabled = provider_type(config_type(soc_scale=0.0, soc_enabled=False))
    q = np.asarray([0.29])
    enabled_point = enabled.evaluate_snapshot(q).point
    disabled_point = disabled.evaluate_snapshot(q).point

    assert np.array_equal(enabled_point.H, disabled_point.H)
    assert np.array_equal(enabled_point.dH_dq, disabled_point.dH_dq)
    assert np.array_equal(enabled_point.connection_q, disabled_point.connection_q)
    assert audit_physical_soc_provider_v220(enabled, q, fermionic=fermionic).passed
    assert audit_physical_soc_provider_v220(disabled, q, fermionic=fermionic).passed


class _WrongSOCDerivativeProvider:
    def __init__(self, base):
        self.base = base
        self.config = base.config
        self.provenance = base.provenance

    @property
    def time_reversal_matrix(self):
        return self.base.time_reversal_matrix

    @property
    def projectors(self):
        return self.base.projectors

    def components(self, q):
        components = self.base.components(q)
        wrong = components.K_soc.copy()
        wrong[0] += 1.0e-3 * np.eye(4)
        return SOCOperatorComponentsV220(
            components.q,
            components.H_spin_free,
            components.K_spin_free,
            components.H_soc,
            wrong,
        ).validate()

    def evaluate_snapshot(self, q):
        return self.base.evaluate_snapshot(q)

    def snapshot_overlap(self, left, right):
        return self.base.snapshot_overlap(left, right)


def test_wrong_but_hermitian_soc_derivative_is_detected():
    provider = _WrongSOCDerivativeProvider(AnalyticSingletTripletSOCProviderV220())
    report = audit_physical_soc_provider_v220(
        provider, np.asarray([0.17]), fermionic=False
    )

    assert not report.passed
    assert not report.checks["K_decomposition"]
    assert not report.checks["SOC_force_derivative"]
    assert report.checks["cross_geometry_differentials"]


def test_general_complex_gauge_transforms_projectors_and_time_reversal_operator():
    provider = AnalyticDoubletSOCProviderV220()
    q = np.asarray([0.21])
    gauge = _gauge()
    G = gauge.matrix(q)
    point = provider.evaluate_snapshot(q).point
    transformed_H = G.conj().T @ point.H @ G
    transformed_J = transform_time_reversal_matrix_v220(
        provider.time_reversal_matrix, G
    )

    assert time_reversal_residual_v220(transformed_H, transformed_J) < 1.0e-12
    assert time_reversal_residual_v220(
        transformed_H, provider.time_reversal_matrix
    ) > 1.0e-6

    vector = np.asarray([0.55 + 0.1j, -0.12j, 0.31 - 0.2j, 0.42 + 0.08j])
    transformed_vector = G.conj().T @ vector
    for projector in provider.projectors.values():
        base_population = projector_population_v220(vector, projector)
        transformed_population = projector_population_v220(
            transformed_vector, transform_projector_v220(projector, G)
        )
        assert transformed_population == pytest.approx(base_population, abs=1.0e-14)


def test_coordinate_dependent_complex_gauge_passes_soc_H_K_D_differentials():
    config = DoubletSOCConfigV220()
    base = AnalyticDoubletSOCProviderV220(config)
    provenance = config.provenance("local_general")
    provider = ContractedElectronicOperatorProviderV213(
        GaugeTransformedOperatorProviderV21(base, _gauge()), provenance
    )

    report = require_provider_differential_contract_v214(
        provider, np.asarray([0.19]), provenance
    )
    assert report.maximum_hamiltonian_derivative_scaled_error < 1.0e-10
    assert report.maximum_connection_scaled_error < 2.0e-9
    assert report.maximum_overlap_isometry_residual < 1.0e-12


def test_incomplete_doublet_and_cross_parity_model_spaces_are_rejected():
    with pytest.raises(ValueError, match="multiplicity 2 requires 2"):
        ElectronicModelSpaceV213(
            name="incomplete doublet",
            representation="fixed_spin_diabatic",
            states=(
                ElectronicStateDescriptorV213(
                    "D(+1/2)", "D", 2, "M=+1/2", 0
                ),
            ),
            complete_multiplets=True,
        ).validate()

    singlet_space = SingletTripletSOCConfigV220().model_space()
    doublet_space = DoubletSOCConfigV220().model_space()
    assert {state.multiplicity % 2 for state in singlet_space.states} == {1}
    assert {state.multiplicity % 2 for state in doublet_space.states} == {0}


def test_broken_kramers_pair_is_detected_by_time_reversal_and_splitting():
    provider = AnalyticDoubletSOCProviderV220()
    H = provider.evaluate_snapshot(np.asarray([0.1])).point.H.copy()
    H[0, 0] += 1.0e-3

    assert time_reversal_residual_v220(H, provider.time_reversal_matrix) > 1.0e-4
    energies = np.linalg.eigvalsh(H)
    assert max(abs(energies[1::2] - energies[0::2])) > 1.0e-5

