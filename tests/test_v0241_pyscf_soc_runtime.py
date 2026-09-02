import numpy as np
import pytest

from gaussian_dynamics.pyscf_soc_runtime_v241 import (
    run_pyscf_oh_static_soc_evidence_v241,
)
from gaussian_dynamics.pyscf_state_interaction_soc_v241 import (
    SpinFreeRootV241,
    assemble_state_interaction_soc_v241,
    probe_pyscf_static_soc_runtime_v241,
    wigner_reduced_transition_density_from_pyscf_ci_v241,
)


_PROBE = probe_pyscf_static_soc_runtime_v241()
pytestmark = pytest.mark.skipif(
    not _PROBE.usable,
    reason="real pinned PySCF 2.13.1 static SOC integration test",
)


@pytest.fixture(scope="module")
def evidence():
    return run_pyscf_oh_static_soc_evidence_v241()


def test_real_oh_bp_somf_returns_direct_complex_soc_elements(evidence):
    matrices = evidence.result.matrices
    H_soc = matrices.H_soc

    assert evidence.passed
    assert H_soc.shape == (6, 6)
    assert np.iscomplexobj(H_soc)
    assert np.linalg.norm(H_soc) > 1.0e-4
    assert np.max(np.abs(H_soc)) > 1.0e-5
    assert matrices.hermiticity_residual < 1.0e-12
    assert matrices.state_order == (
        "D1(M=+1/2)",
        "D1(M=-1/2)",
        "D2(M=+1/2)",
        "D2(M=-1/2)",
        "D3(M=+1/2)",
        "D3(M=-1/2)",
    )


def test_real_oh_integrals_and_state_interaction_have_independent_gates(evidence):
    checks = evidence.audit.checks
    metrics = evidence.audit.metrics

    assert len(checks) == 39
    assert all(type(value) is bool for value in checks.values())
    assert all(checks.values())
    assert checks["explicit_somf_matches_independent_pyscf_jk_path"]
    assert metrics["somf_jk_crosscheck_max_abs_error"] < 1.0e-12
    assert metrics["H_soc_frobenius_norm_cm_inverse"] > 100.0
    assert metrics["H_soc_max_abs_cm_inverse"] > 40.0


def test_real_oh_doublet_obeys_time_reversal_and_kramers_pairing(evidence):
    matrices = evidence.result.matrices

    assert matrices.time_reversal_residual < 1.0e-12
    assert matrices.time_reversal_square_residual < 1.0e-12
    assert matrices.maximum_kramers_pair_splitting_hartree < 1.0e-10
    assert evidence.audit.checks["kramers_or_even_sector_condition"]


def test_static_result_cannot_be_used_as_a_trajectory_provider(evidence):
    result = evidence.result

    assert result.capabilities.tier == "static_soc"
    assert result.static_soc_admitted
    assert not result.trajectory_ready
    assert not result.molecular_soc_contract.real_backend_admission_ready
    assert not evidence.claims["live_molecular_SOC_backend_admitted"]
    assert not evidence.claims["physical_SOC_derivatives_validated"]
    assert not evidence.claims["cross_geometry_SOC_tracking_validated"]
    assert not evidence.claims["ab_initio_SOC_accuracy_validated"]


def test_mixed_singlet_triplet_wigner_reduction_uses_a_nonzero_lifted_component():
    # Two orthogonal, exactly spin-pure 2e/3o CI roots in a common real orbital
    # basis.  The triplet M=0 diagonal q=0 Clebsch-Gordan coefficient vanishes,
    # so extraction must raise that root to M=+1 instead of dividing by zero.
    singlet = np.zeros((3, 3))
    singlet[0, 1] = singlet[1, 0] = 1.0 / np.sqrt(2.0)
    triplet = np.zeros((3, 3))
    triplet[0, 2] = 1.0 / np.sqrt(2.0)
    triplet[2, 0] = -1.0 / np.sqrt(2.0)
    roots = (
        SpinFreeRootV241("S", 0.0, 0, 0, 0.0),
        SpinFreeRootV241("T", 0.01, 2, 0, 2.0),
    )
    reduced, sample_ms, coefficients = (
        wigner_reduced_transition_density_from_pyscf_ci_v241(
            (singlet, triplet),
            roots,
            ncore=0,
            ncas=3,
            nmo=3,
            nelecas=(1, 1),
        )
    )
    antisymmetric = np.zeros((3, 3))
    antisymmetric[1, 2] = 1.0
    antisymmetric[2, 1] = -1.0
    soc_integrals = -1j * np.stack(
        (antisymmetric, 2.0 * antisymmetric, 3.0 * antisymmetric)
    ) * 1.0e-3
    matrices = assemble_state_interaction_soc_v241(
        roots, reduced, soc_integrals
    )

    assert sample_ms[0, 1] == 0
    assert sample_ms[1, 1] == 2
    assert coefficients[0, 1] != 0.0
    assert coefficients[1, 1] != 0.0
    assert matrices.state_order == (
        "S(M=+0)",
        "T(M=+1)",
        "T(M=+0)",
        "T(M=-1)",
    )
    assert np.linalg.norm(matrices.H_soc) > 1.0e-3
    assert matrices.hermiticity_residual < 1.0e-12
    assert matrices.time_reversal_residual < 1.0e-12
