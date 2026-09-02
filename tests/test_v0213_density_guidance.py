import numpy as np
import pytest

from gaussian_dynamics.block_sparse_molecular_v21 import BlockMolecularTBFV21
from gaussian_dynamics.complex_gauge_v21 import (
    GaugeTransformedOperatorProviderV21,
    PhaseMixingGaugeV21,
)
from gaussian_dynamics.density_guidance_v213 import BlockDensityMatrixGuidanceV213
from gaussian_dynamics.electronic_operator_v21 import (
    ElectronicOperatorPointV21,
    ElectronicOperatorSnapshotV21,
)
from gaussian_dynamics.self_consistent_block_v212 import MeanFieldGuidanceSettingsV212
from gaussian_dynamics.self_consistent_block_v213 import (
    SelfConsistentBlockSettingsV213,
    run_self_consistent_block_dynamics_v213,
)
from gaussian_dynamics.synthetic_operator_provider_v21 import (
    SyntheticLinearOperatorConfigV21,
    SyntheticLinearOperatorProviderV21,
)


class _ExactlyDegenerateProvider:
    def evaluate_snapshot(self, q):
        q = np.asarray(q, dtype=float)
        K = np.asarray([np.diag([1.0, -1.0])], dtype=complex)
        point = ElectronicOperatorPointV21(
            q=q.copy(),
            H=np.zeros((2, 2), dtype=complex),
            dH_dq=K,
            connection_q=np.zeros_like(K),
            mass_matrix_q_au=np.asarray([[20.0]]),
            metadata={"exact_degeneracy": True},
        ).validate()
        return ElectronicOperatorSnapshotV21(
            point=point,
            state_vectors=np.eye(2, dtype=complex),
        ).validate()

    @staticmethod
    def snapshot_overlap(left, right):
        return left.state_vectors.conj().T @ right.state_vectors


def _gauge():
    hadamard = np.asarray([[1.0, 1.0], [1.0, -1.0]], dtype=complex) / np.sqrt(2.0)
    return PhaseMixingGaugeV21(
        U0=hadamard,
        phase_gradient=np.asarray([[0.23], [-0.17]]),
        phase_offset=np.asarray([0.31, -0.22]),
    )


def _tbf(uid, q):
    return BlockMolecularTBFV21(uid, np.asarray([q]), np.asarray([0.0]), np.asarray([[1.2]]))


def test_density_guidance_is_gauge_invariant_at_exact_degeneracy_and_zero_amplitude():
    base_provider = _ExactlyDegenerateProvider()
    gauge = _gauge()
    gauge_provider = GaugeTransformedOperatorProviderV21(base_provider, gauge)
    base_guide = BlockDensityMatrixGuidanceV213()
    gauge_guide = BlockDensityMatrixGuidanceV213()

    initial = _tbf(7, 0.0)
    c = np.asarray([1.0, 0.0], dtype=complex)
    c_gauge = gauge.matrix(initial.q).conj().T @ c
    f_base, _, _ = base_guide.forces_and_masses([initial], c, base_provider, 2)
    f_gauge, _, _ = gauge_guide.forces_and_masses(
        [initial], c_gauge, gauge_provider, 2
    )
    assert np.allclose(f_base, [[-1.0]], rtol=0.0, atol=1.0e-13)
    assert np.allclose(f_gauge, f_base, rtol=0.0, atol=1.0e-13)

    moved = _tbf(7, 0.4)
    zeros = np.zeros(2, dtype=complex)
    retained_base, _, _ = base_guide.forces_and_masses(
        [moved], zeros, base_provider, 2
    )
    retained_gauge, _, _ = gauge_guide.forces_and_masses(
        [moved], zeros, gauge_provider, 2
    )
    assert np.allclose(retained_base, [[-1.0]], rtol=0.0, atol=1.0e-13)
    assert np.allclose(retained_gauge, retained_base, rtol=0.0, atol=1.0e-13)
    assert gauge_guide.diagnostics_dict()["density_transports"] == 1
    assert gauge_guide.diagnostics_dict()["retained_density_uses"] == 1


def test_unseeded_zero_block_has_zero_force_and_parent_density_can_be_inherited():
    provider = _ExactlyDegenerateProvider()
    guide = BlockDensityMatrixGuidanceV213()
    parent = _tbf(1, -0.2)
    child = _tbf(2, 0.3)

    unseeded_force, _, _ = guide.forces_and_masses(
        [parent], np.zeros(2, dtype=complex), provider, 2
    )
    assert np.array_equal(unseeded_force, np.zeros((1, 1)))
    assert guide.density(parent.uid) is None

    guide.forces_and_masses(
        [parent], np.asarray([0.0, 1.0], dtype=complex), provider, 2
    )
    guide.on_insert(child, provider, parent_uid=parent.uid)
    assert np.allclose(guide.density(child.uid), guide.density(parent.uid))
    child_force, _, _ = guide.forces_and_masses(
        [child], np.zeros(2, dtype=complex), provider, 2
    )
    assert np.allclose(child_force, [[1.0]], rtol=0.0, atol=1.0e-13)
    assert guide.diagnostics_dict()["parent_inheritances"] == 1
    guide.on_prune(child.uid)
    assert guide.density(child.uid) is None


def test_historical_lowest_eigenvector_fallback_is_explicitly_unavailable():
    with pytest.raises(ValueError, match="gauge-dependent lowest-eigenvector"):
        MeanFieldGuidanceSettingsV212(
            low_amplitude_policy="lowest_eigenvector"
        ).validate()


def test_v0213_runner_installs_density_guidance_in_the_propagation_path():
    provider = SyntheticLinearOperatorProviderV21(
        SyntheticLinearOperatorConfigV21(
            nstate=2, nq=1, seed=21309, mass=25.0
        )
    )
    basis = [_tbf(11, -0.15)]
    output = run_self_consistent_block_dynamics_v213(
        basis,
        np.asarray([0.8 + 0.1j, -0.2 + 0.3j]),
        provider,
        dt=0.002,
        steps=3,
        store_every=1,
        settings=SelfConsistentBlockSettingsV213(
            corrector_iterations=2,
            momentum_tolerance=1.0e-12,
        ),
    )
    assert output["release_path"] == "v0.21.3"
    assert output["settings"]["guidance_contract"] == (
        "transported electronic density matrix"
    )
    assert output["guidance_diagnostics"]["coefficient_refreshes"] > 0
    assert output["guidance_trial_state_rollbacks"] > 0
    assert output["maximum_norm_drift"] < 1.0e-11
