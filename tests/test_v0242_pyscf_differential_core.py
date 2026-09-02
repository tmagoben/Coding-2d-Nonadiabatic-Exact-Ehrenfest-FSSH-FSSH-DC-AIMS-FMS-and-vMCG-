from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

import gaussian_dynamics.pyscf_differential_soc_v242 as differential
from gaussian_dynamics.pyscf_state_interaction_soc_v241 import (
    SpinFreeRootV241,
    complete_spin_microstates_v241,
    time_reversal_matrix_v241,
)


def _two_doublet_roots():
    return (
        SpinFreeRootV241("D1", -75.0, 1, 1, 0.75).validate(),
        SpinFreeRootV241("D2", -74.9, 1, 1, 0.75).validate(),
    )


def test_complete_multiplet_overlap_lifts_root_matrix_without_mixing_ms(monkeypatch):
    roots = _two_doublet_roots()
    states = complete_spin_microstates_v241(roots)
    root_overlap = np.array([[0.8, 0.3], [-0.3, 0.8]], dtype=complex)
    monkeypatch.setattr(
        differential,
        "casscf_state_overlap_matrix",
        lambda left, right: root_overlap.copy(),
    )

    class FakeSnapshot:
        def __init__(self):
            self.roots = roots
            self.matrices = SimpleNamespace(microstates=states)
            self.wavefunction_snapshot = object()

        def validate(self):
            return self

    observed = differential.complete_multiplet_overlap_v242(
        FakeSnapshot(), FakeSnapshot()
    )
    expected = np.array(
        [
            [0.8, 0.0, 0.3, 0.0],
            [0.0, 0.8, 0.0, 0.3],
            [-0.3, 0.0, 0.8, 0.0],
            [0.0, -0.3, 0.0, 0.8],
        ],
        dtype=complex,
    )
    assert np.allclose(observed, expected)


def test_root_phase_alignment_shares_one_phase_across_each_complete_doublet():
    states = complete_spin_microstates_v241(_two_doublet_roots())
    phases = np.diag([1j, 1j, -1.0, -1.0]).astype(complex)
    aligned, correction = differential.phase_align_complete_multiplet_overlap_v242(
        phases, states
    )
    assert np.allclose(aligned, np.eye(4))
    assert np.allclose(correction[:2], -1j)
    assert np.allclose(correction[2:], -1.0)


def test_root_phase_alignment_rejects_degenerate_rotation_as_ambiguous():
    states = complete_spin_microstates_v241(_two_doublet_roots())
    rotation = np.kron(
        np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
        np.eye(2),
    )
    with pytest.raises(ValueError, match="phase is ambiguous"):
        differential.phase_align_complete_multiplet_overlap_v242(rotation, states)


def test_polar_transported_component_differences_remove_degenerate_endpoint_gauges(
    monkeypatch,
):
    roots = _two_doublet_roots()
    states = complete_spin_microstates_v241(roots)
    J = time_reversal_matrix_v241(states)
    identity = np.eye(2)
    h = 0.02
    H0_sf_root = np.diag([-75.0, -74.9])
    H0_soc_root = np.array([[0.0, 2.0e-4], [2.0e-4, 0.0]])
    K_sf_root = np.array([[0.3, 0.04], [0.04, -0.2]])
    K_soc_root = np.array([[0.0, 7.0e-5], [7.0e-5, 0.0]])

    def lift(root_matrix):
        return np.kron(root_matrix, identity).astype(complex)

    angle_minus, angle_plus = 0.71, -0.43
    root_gauges = {
        "minus": np.array(
            [
                [np.cos(angle_minus), -np.sin(angle_minus)],
                [np.sin(angle_minus), np.cos(angle_minus)],
            ]
        ),
        "plus": np.array(
            [
                [np.cos(angle_plus), -np.sin(angle_plus)],
                [np.sin(angle_plus), np.cos(angle_plus)],
            ]
        ),
    }
    full_gauges = {name: lift(value) for name, value in root_gauges.items()}

    class FakeMatrices:
        def __init__(self, H_sf, H_soc):
            self.microstates = states
            self.H_spin_free = np.asarray(H_sf, dtype=complex)
            self.H_soc = np.asarray(H_soc, dtype=complex)
            self.H_total = self.H_spin_free + self.H_soc
            self.time_reversal_matrix = J

    class FakeSnapshot:
        def __init__(self, tag, matrices):
            self.tag = tag
            self.roots = roots
            self.matrices = matrices

        def validate(self):
            return self

        def fingerprint(self):
            return (self.tag * 64)[:64]

    center = FakeSnapshot("c", FakeMatrices(lift(H0_sf_root), lift(H0_soc_root)))
    endpoints = {}
    for tag, fingerprint_tag, sign in (("minus", "b", -1.0), ("plus", "d", 1.0)):
        gauge = full_gauges[tag]
        physical_sf = lift(H0_sf_root + sign * h * K_sf_root)
        physical_soc = lift(H0_soc_root + sign * h * K_soc_root)
        endpoints[tag] = FakeSnapshot(
            fingerprint_tag,
            FakeMatrices(
                gauge.conj().T @ physical_sf @ gauge,
                gauge.conj().T @ physical_soc @ gauge,
            ),
        )

    monkeypatch.setattr(
        differential,
        "complete_multiplet_overlap_v242",
        lambda left, right: full_gauges[right.tag == "b" and "minus" or "plus"],
    )
    record = differential.transported_soc_central_difference_v242(
        center,
        endpoints["minus"],
        endpoints["plus"],
        displacement_bohr=h,
        coordinate_label="synthetic",
    )
    assert np.allclose(record.K_spin_free, lift(K_sf_root), atol=1.0e-11)
    assert np.allclose(record.K_soc, lift(K_soc_root), atol=1.0e-11)
    assert np.allclose(record.K_total, lift(K_sf_root + K_soc_root), atol=1.0e-11)
    assert np.allclose(record.derivative_connection, 0.0, atol=1.0e-12)


def test_transport_record_rejects_nonhermitian_component_without_symmetrizing():
    matrix = np.zeros((2, 2), dtype=complex)
    record = differential.TransportedSOCDerivativeV242(
        coordinate_label="q",
        displacement_bohr=0.1,
        center_fingerprint="a" * 64,
        minus_fingerprint="b" * 64,
        plus_fingerprint="c" * 64,
        overlap_center_minus=np.eye(2),
        overlap_center_plus=np.eye(2),
        transport_minus_to_center=np.eye(2),
        transport_plus_to_center=np.eye(2),
        H_spin_free_minus_to_center=matrix,
        H_spin_free_plus_to_center=matrix,
        H_soc_minus_to_center=matrix,
        H_soc_plus_to_center=matrix,
        K_spin_free=matrix,
        K_soc=matrix,
        K_total=matrix,
        derivative_connection=matrix,
        overlap_metrics={},
        residuals={},
    ).validate()
    broken = replace(
        record,
        K_soc=np.array([[0.0, 1.0], [0.0, 0.0]], dtype=complex),
        K_total=np.array([[0.0, 1.0], [0.0, 0.0]], dtype=complex),
    )
    with pytest.raises(ValueError, match="K_soc is not Hermitian"):
        broken.validate()


def test_transport_record_binds_each_component_to_its_endpoint_matrices():
    matrix = np.zeros((2, 2), dtype=complex)
    record = differential.TransportedSOCDerivativeV242(
        coordinate_label="q",
        displacement_bohr=0.1,
        center_fingerprint="a" * 64,
        minus_fingerprint="b" * 64,
        plus_fingerprint="c" * 64,
        overlap_center_minus=np.eye(2),
        overlap_center_plus=np.eye(2),
        transport_minus_to_center=np.eye(2),
        transport_plus_to_center=np.eye(2),
        H_spin_free_minus_to_center=matrix,
        H_spin_free_plus_to_center=matrix,
        H_soc_minus_to_center=matrix,
        H_soc_plus_to_center=matrix,
        K_spin_free=matrix,
        K_soc=matrix,
        K_total=matrix,
        derivative_connection=matrix,
        overlap_metrics={},
        residuals={},
    ).validate()
    compensating_shift = np.diag([0.2, -0.2]).astype(complex)
    broken = replace(
        record,
        K_spin_free=record.K_spin_free + compensating_shift,
        K_soc=record.K_soc - compensating_shift,
    )
    with pytest.raises(ValueError, match="spin-free derivative disagrees"):
        broken.validate()


def test_oh_scan_rejects_nonpositive_or_nonhalving_step_lists_before_runtime():
    with pytest.raises(ValueError, match="positive finite steps"):
        differential.run_pyscf_oh_bond_differential_soc_v242(
            steps_bohr=(0.08, 0.04, 0.0)
        )
    with pytest.raises(ValueError, match="strictly decreasing"):
        differential.run_pyscf_oh_bond_differential_soc_v242(
            steps_bohr=(0.08, 0.04, 0.04)
        )
