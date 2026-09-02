import numpy as np
import pytest

from gaussian_dynamics.electronic_contract_v213 import (
    ElectronicModelSpaceV213,
    ElectronicOperatorProvenanceV213,
    ElectronicStateDescriptorV213,
    compose_electronic_operator_v213,
    hartree_to_wavenumber_v213,
    validate_electronic_contract_v213,
    wavenumber_to_hartree_v213,
)
from gaussian_dynamics.electronic_operator_v21 import ElectronicOperatorPointV21


def _complete_spin_model_space():
    return ElectronicModelSpaceV213(
        name="S0-plus-T1-components",
        representation="fixed_spin_diabatic",
        states=(
            ElectronicStateDescriptorV213("S0", "S0", 1, "M=0", 0),
            ElectronicStateDescriptorV213("T1(M=-1)", "T1", 3, "M=-1", 0),
            ElectronicStateDescriptorV213("T1(M=0)", "T1", 3, "M=0", 0),
            ElectronicStateDescriptorV213("T1(M=+1)", "T1", 3, "M=+1", 0),
        ),
        complete_multiplets=True,
    ).validate()


def _provenance(**parameters):
    return ElectronicOperatorProvenanceV213(
        model_name="four-state-zero-SOC fixture",
        model_version="1",
        model_space=_complete_spin_model_space(),
        spin_free_method="analytic fixture",
        parameters=parameters,
    ).validate()


def test_strict_validator_rejects_defect_hidden_by_default_allclose_rtol():
    point = ElectronicOperatorPointV21(
        q=np.array([0.0]),
        H=np.diag([1.0 + 1.0e-6j, 2.0]),
        dH_dq=np.zeros((1, 2, 2), dtype=complex),
        connection_q=np.zeros((1, 2, 2), dtype=complex),
        mass_matrix_q_au=np.eye(1),
    )
    with pytest.raises(ValueError, match="Hamiltonian Hermiticity residual"):
        point.validate(atol=1.0e-12)


def test_zero_soc_composition_freezes_total_H_K_and_units():
    provenance = _provenance(coupling=0.0125)
    H0 = np.diag([0.0, 0.02, 0.02, 0.02]).astype(complex)
    K0 = np.zeros((2, 4, 4), dtype=complex)
    K0[0] = np.diag([0.01, -0.02, -0.02, -0.02])
    K0[1, 0, 1] = 0.003 + 0.004j
    K0[1, 1, 0] = K0[1, 0, 1].conjugate()

    point = compose_electronic_operator_v213(
        q=np.array([0.1, -0.2]),
        H_spin_free=H0,
        dH_spin_free_dq=K0,
        connection_q=np.zeros_like(K0),
        mass_matrix_q_au=np.diag([1500.0, 1800.0]),
        provenance=provenance,
    )

    assert np.array_equal(point.H, H0)
    assert np.array_equal(point.hamiltonian_derivative_operator_q, K0)
    assert point.metadata["operator_decomposition"]["soc_enabled"] is False
    assert point.metadata["v213_provenance_fingerprint"] == provenance.fingerprint()
    values = np.array([0.0, 125.0, 1000.0])
    assert np.allclose(
        hartree_to_wavenumber_v213(wavenumber_to_hartree_v213(values)),
        values,
        rtol=0.0,
        atol=2.0e-13,
    )


def test_nonzero_soc_terms_require_explicit_soc_provenance():
    provenance = _provenance()
    with pytest.raises(ValueError, match="nonzero SOC terms"):
        compose_electronic_operator_v213(
            q=np.array([0.0]),
            H_spin_free=np.zeros((4, 4)),
            dH_spin_free_dq=np.zeros((1, 4, 4)),
            H_soc=np.diag([0.0, 1.0e-4, 0.0, 0.0]),
            connection_q=np.zeros((1, 4, 4)),
            mass_matrix_q_au=np.eye(1),
            provenance=provenance,
        )


def test_fixed_frame_contract_rejects_nonzero_connection():
    provenance = _provenance()
    point = ElectronicOperatorPointV21(
        q=np.array([0.0]),
        H=np.zeros((4, 4), dtype=complex),
        dH_dq=np.zeros((1, 4, 4), dtype=complex),
        connection_q=np.asarray([1.0j * np.eye(4)]),
        mass_matrix_q_au=np.eye(1),
    ).validate()
    with pytest.raises(ValueError, match="fixed electronic frame"):
        validate_electronic_contract_v213(point, provenance)


def test_complete_multiplet_contract_rejects_missing_component():
    incomplete = ElectronicModelSpaceV213(
        name="incomplete triplet",
        representation="fixed_spin_diabatic",
        states=(
            ElectronicStateDescriptorV213("T(M=-1)", "T", 3, "M=-1"),
            ElectronicStateDescriptorV213("T(M=0)", "T", 3, "M=0"),
        ),
        complete_multiplets=True,
    )
    with pytest.raises(ValueError, match="requires 3"):
        incomplete.validate()


def test_provenance_fingerprint_is_order_stable_and_parameter_sensitive():
    first = _provenance(alpha=1, nested={"z": 2, "a": [3, 4]})
    reordered = _provenance(nested={"a": [3, 4], "z": 2}, alpha=1)
    changed = _provenance(alpha=2, nested={"z": 2, "a": [3, 4]})
    assert first.fingerprint() == reordered.fingerprint()
    assert first.fingerprint() != changed.fingerprint()
    with pytest.raises(ValueError, match="non-finite"):
        _provenance(bad=np.nan)


def test_operator_contract_never_silently_discards_complex_coordinates_or_mass():
    common = dict(
        H=np.zeros((2, 2)),
        dH_dq=np.zeros((1, 2, 2)),
        connection_q=np.zeros((1, 2, 2)),
    )
    with pytest.raises(ValueError, match="q must be real"):
        ElectronicOperatorPointV21(
            q=np.asarray([0.1 + 1.0e-9j]),
            mass_matrix_q_au=np.eye(1),
            **common,
        ).validate()
    with pytest.raises(ValueError, match="mass_matrix_q_au must be real"):
        ElectronicOperatorPointV21(
            q=np.asarray([0.1]),
            mass_matrix_q_au=np.asarray([[1.0 + 1.0e-9j]]),
            **common,
        ).validate()
