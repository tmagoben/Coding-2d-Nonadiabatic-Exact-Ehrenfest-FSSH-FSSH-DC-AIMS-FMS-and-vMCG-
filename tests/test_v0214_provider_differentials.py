import numpy as np
import pytest

from gaussian_dynamics.complex_gauge_v21 import (
    GaugeTransformedOperatorProviderV21,
    PhaseMixingGaugeV21,
    random_unitary_v21,
)
from gaussian_dynamics.electronic_contract_v213 import (
    ContractedElectronicOperatorProviderV213,
    ElectronicModelSpaceV213,
    ElectronicOperatorProvenanceV213,
    ElectronicStateDescriptorV213,
)
from gaussian_dynamics.electronic_operator_v21 import (
    ElectronicOperatorPointV21,
    ElectronicOperatorSnapshotV21,
)
from gaussian_dynamics.provider_differential_audit_v214 import (
    ProviderDifferentialAuditSettingsV214,
    audit_provider_differentials_v214,
    require_provider_differential_contract_v214,
)
from gaussian_dynamics.synthetic_operator_provider_v21 import (
    SyntheticLinearOperatorConfigV21,
    SyntheticLinearOperatorProviderV21,
)


def _provenance(representation):
    space = ElectronicModelSpaceV213(
        name=f"three-state {representation} differential fixture",
        representation=representation,
        states=tuple(ElectronicStateDescriptorV213(f"state-{i}") for i in range(3)),
    )
    return ElectronicOperatorProvenanceV213(
        model_name="v0.21.4 differential fixture",
        model_version="1",
        model_space=space,
        spin_free_method="analytic complex linear model",
        derivative_method="analytic physical operator derivative",
        parameters={"seed": 21401},
    )


def _base():
    return SyntheticLinearOperatorProviderV21(
        SyntheticLinearOperatorConfigV21(
            nstate=3,
            nq=2,
            mass=30.0,
            seed=21401,
            base_scale=0.04,
            derivative_scale=0.015,
        )
    )


def _gauge():
    return PhaseMixingGaugeV21(
        U0=random_unitary_v21(3, 21402),
        phase_gradient=np.asarray([[0.23, -0.11], [-0.17, 0.08], [0.05, 0.19]]),
        phase_offset=np.asarray([0.31, -0.22, 0.14]),
    )


def test_fixed_frame_provider_passes_cross_geometry_H_K_D_certification():
    provenance = _provenance("fixed_general")
    provider = ContractedElectronicOperatorProviderV213(_base(), provenance)
    report = audit_provider_differentials_v214(
        provider,
        np.asarray([0.17, -0.28]),
        provenance,
    )

    assert report.passed
    assert report.maximum_hamiltonian_derivative_scaled_error < 1.0e-10
    assert report.maximum_connection_scaled_error == 0.0
    assert report.maximum_overlap_isometry_residual == 0.0
    assert len(report.rows) == 2


def test_differential_audit_requires_an_explicit_provider_fingerprint():
    provenance = _provenance("fixed_general")
    report = audit_provider_differentials_v214(
        _base(), np.asarray([0.17, -0.28]), provenance
    )

    assert not report.passed
    assert not report.checks["provenance_fingerprint_consistency"]


def test_coordinate_dependent_complex_gauge_passes_K_and_D_certification():
    provenance = _provenance("local_general")
    provider = ContractedElectronicOperatorProviderV213(
        GaugeTransformedOperatorProviderV21(_base(), _gauge()), provenance
    )
    report = require_provider_differential_contract_v214(
        provider,
        np.asarray([0.17, -0.28]),
        provenance,
        settings=ProviderDifferentialAuditSettingsV214(
            default_step=5.0e-5,
            hamiltonian_derivative_tolerance=2.0e-9,
            connection_tolerance=2.0e-9,
        ),
    )

    assert report.maximum_hamiltonian_derivative_scaled_error < 1.0e-10
    assert report.maximum_connection_scaled_error < 2.0e-10
    assert report.maximum_overlap_isometry_residual < 1.0e-12


class _PerturbedOperatorProvider:
    def __init__(self, base, *, perturb_K=False, erase_D=False):
        self.base = base
        self.perturb_K = bool(perturb_K)
        self.erase_D = bool(erase_D)

    def evaluate_snapshot(self, q):
        snapshot = self.base.evaluate_snapshot(q)
        K = snapshot.point.dH_dq.copy()
        D = snapshot.point.connection_q.copy()
        if self.perturb_K:
            K[0] += 1.0e-3 * np.eye(snapshot.point.nstate)
        if self.erase_D:
            D[:] = 0.0
        point = ElectronicOperatorPointV21(
            q=snapshot.point.q.copy(),
            H=snapshot.point.H.copy(),
            dH_dq=K,
            connection_q=D,
            mass_matrix_q_au=snapshot.point.mass_matrix_q_au.copy(),
            metadata=dict(snapshot.point.metadata),
        ).validate()
        return ElectronicOperatorSnapshotV21(
            point=point,
            state_vectors=snapshot.state_vectors.copy(),
            parent_snapshot=snapshot,
        ).validate()

    def snapshot_overlap(self, left, right):
        return self.base.snapshot_overlap(left.parent_snapshot, right.parent_snapshot)


def test_differential_audit_rejects_pointwise_valid_but_inconsistent_K():
    provenance = _provenance("fixed_general")
    contracted = ContractedElectronicOperatorProviderV213(_base(), provenance)
    bad = _PerturbedOperatorProvider(contracted, perturb_K=True)

    report = audit_provider_differentials_v214(
        bad, np.asarray([0.17, -0.28]), provenance
    )
    assert not report.passed
    assert report.checks["structural_invariants"]
    assert not report.checks["physical_H_derivatives"]
    with pytest.raises(ValueError, match="physical_H_derivatives"):
        require_provider_differential_contract_v214(
            bad, np.asarray([0.17, -0.28]), provenance
        )


def test_differential_audit_rejects_pointwise_valid_but_inconsistent_D():
    provenance = _provenance("local_general")
    contracted = ContractedElectronicOperatorProviderV213(
        GaugeTransformedOperatorProviderV21(_base(), _gauge()), provenance
    )
    bad = _PerturbedOperatorProvider(contracted, erase_D=True)

    report = audit_provider_differentials_v214(
        bad, np.asarray([0.17, -0.28]), provenance
    )
    assert not report.passed
    assert report.checks["structural_invariants"]
    assert report.checks["physical_H_derivatives"]
    assert not report.checks["derivative_connections"]
