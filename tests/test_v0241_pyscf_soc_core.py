from dataclasses import replace
import math

import numpy as np
import pytest

from gaussian_dynamics.pyscf_state_interaction_soc_v241 import (
    BPSOMFIntegralsV241,
    SpinFreeRootV241,
    assemble_state_interaction_soc_v241,
    clebsch_gordan_twice_v241,
    complete_spin_microstates_v241,
    probe_pyscf_static_soc_runtime_v241,
    root_projectors_v241,
    time_reversal_matrix_v241,
)


def test_integer_twice_quantum_number_clebsch_gordan_values():
    assert clebsch_gordan_twice_v241(1, 1, 2, 0, 1, 1) == pytest.approx(
        1.0 / math.sqrt(3.0), abs=1.0e-15
    )
    assert clebsch_gordan_twice_v241(2, 0, 2, 0, 2, 0) == 0.0
    assert clebsch_gordan_twice_v241(0, 0, 2, 0, 2, 0) == 1.0
    assert clebsch_gordan_twice_v241(2, 0, 2, 0, 0, 0) == pytest.approx(
        -1.0 / math.sqrt(3.0), abs=1.0e-15
    )
    assert clebsch_gordan_twice_v241(1, 1, 2, 2, 1, 1) == 0.0
    with pytest.raises(TypeError, match="integer twice"):
        clebsch_gordan_twice_v241(0.5, 0, 2, 0, 2, 0)


def test_doublet_quartet_microstate_order_is_complete_and_time_reversal_is_fermionic():
    roots = (
        SpinFreeRootV241("D", -1.0, 1, 1, 0.75),
        SpinFreeRootV241("Q", -0.8, 3, 1, 3.75),
    )
    states = complete_spin_microstates_v241(roots)
    labels = tuple(state.label for state in states)
    J = time_reversal_matrix_v241(states)

    assert labels == (
        "D(M=+1/2)",
        "D(M=-1/2)",
        "Q(M=+3/2)",
        "Q(M=+1/2)",
        "Q(M=-1/2)",
        "Q(M=-3/2)",
    )
    assert np.allclose(J.conj().T @ J, np.eye(6), atol=1.0e-15)
    assert np.allclose(J @ J.conj(), -np.eye(6), atol=1.0e-15)
    assert set(root_projectors_v241(states)) == {"D", "Q"}


def test_singlet_triplet_microstate_space_is_bosonic_and_static_assembly_is_dimension_safe():
    roots = (
        SpinFreeRootV241("S", -1.0, 0, 0, 0.0),
        SpinFreeRootV241("T", -0.9, 2, 0, 2.0),
    )
    states = complete_spin_microstates_v241(roots)
    reduced = np.zeros((2, 2, 3, 3), dtype=complex)
    integrals = np.zeros((3, 3, 3), dtype=complex)
    matrices = assemble_state_interaction_soc_v241(roots, reduced, integrals)
    J = matrices.time_reversal_matrix

    assert matrices.state_order == (
        "S(M=+0)",
        "T(M=+1)",
        "T(M=+0)",
        "T(M=-1)",
    )
    assert matrices.H_soc.shape == (4, 4)
    assert np.array_equal(matrices.H_soc, np.zeros((4, 4)))
    assert np.allclose(J @ J.conj(), np.eye(4), atol=1.0e-15)
    assert matrices.maximum_kramers_pair_splitting_hartree is None


def test_cross_electron_parity_multiplets_are_rejected():
    roots = (
        SpinFreeRootV241("S", 0.0, 0, 0),
        SpinFreeRootV241("D", 0.1, 1, 1),
    )
    with pytest.raises(ValueError, match="even- and odd-electron"):
        complete_spin_microstates_v241(roots)


def test_nonhermitian_direct_soc_input_is_not_silently_symmetrized():
    roots = (
        SpinFreeRootV241("D1", 0.0, 1, 1),
        SpinFreeRootV241("D2", 0.1, 1, 1),
    )
    reduced = np.zeros((2, 2, 2, 2), dtype=complex)
    reduced[0, 1, 0, 1] = 1.0
    orbital = np.asarray([[0.0, -1j], [1j, 0.0]])
    integrals = np.stack((orbital, 0.5 * orbital, 0.25 * orbital))

    with pytest.raises(ValueError, match="not Hermitian"):
        assemble_state_interaction_soc_v241(roots, reduced, integrals)


def test_bp_somf_integral_prefactor_and_antisymmetry_fail_closed():
    antisymmetric = np.asarray([[0.0, 1.0], [-1.0, 0.0]])
    one = np.stack((antisymmetric, 2.0 * antisymmetric, -antisymmetric))
    two = 0.25 * one
    c = 137.035999
    prefactor = 0.5 / c**2
    effective_ao = prefactor * (one - two)
    effective_mo = -1j * effective_ao
    density = np.eye(2)
    valid = BPSOMFIntegralsV241(
        one,
        two,
        effective_ao,
        effective_mo,
        density,
        density,
        c,
        prefactor,
    ).validate()

    assert valid.one_electron_antisymmetry_residual == 0.0
    with pytest.raises(ValueError, match=r"0.5/c\^2"):
        replace(valid, prefactor=2.0 * prefactor).validate()
    corrupted = one.copy()
    corrupted[0, 0, 0] = 1.0
    with pytest.raises(ValueError, match="not antisymmetric"):
        replace(valid, one_electron_ao_cartesian=corrupted).validate()


def test_static_runtime_probe_never_infers_trajectory_capabilities():
    probe = probe_pyscf_static_soc_runtime_v241()

    assert type(probe.usable) is bool
    assert probe.required_version == "2.13.1"
    if probe.installed:
        assert probe.module_version is not None
    else:
        assert not probe.usable
