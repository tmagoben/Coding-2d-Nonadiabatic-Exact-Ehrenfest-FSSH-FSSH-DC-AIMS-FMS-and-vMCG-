import numpy as np
import pytest

from gaussian_dynamics.checkpoint_restart_v214 import (
    SelfConsistentBlockSettingsV214,
    run_self_consistent_block_dynamics_v214,
)
from gaussian_dynamics.electronic_contract_v213 import (
    ContractedElectronicOperatorProviderV213,
    ElectronicModelSpaceV213,
    ElectronicOperatorProvenanceV213,
    ElectronicStateDescriptorV213,
)
from gaussian_dynamics.provider_differential_audit_v214 import (
    require_provider_differential_contract_v214,
)
from gaussian_dynamics.synthetic_operator_provider_v21 import (
    SyntheticLinearOperatorConfigV21,
    SyntheticLinearOperatorProviderV21,
)
from gaussian_dynamics.zero_soc_rehearsal_v214 import (
    ZeroSOCRehearsalProviderV214,
    audit_zero_soc_equivalence_v214,
)
from tests.test_v0214_checkpoint_restart import (
    _basis,
    _coefficients,
    _phase_aligned_metric_error,
)


def _provenance(*, soc_enabled=False):
    return ElectronicOperatorProvenanceV213(
        model_name="v0.21.4 explicit zero-SOC rehearsal",
        model_version="1",
        model_space=ElectronicModelSpaceV213(
            name="two-state fixed zero-SOC fixture",
            representation="fixed_general",
            states=(
                ElectronicStateDescriptorV213("state-0"),
                ElectronicStateDescriptorV213("state-1"),
            ),
        ),
        spin_free_method="analytic linear fixture",
        soc_enabled=soc_enabled,
        soc_method="analytic test" if soc_enabled else "none",
        derivative_method="analytic physical operator derivative",
        parameters={"seed": 21430},
    )


def _base():
    return SyntheticLinearOperatorProviderV21(
        SyntheticLinearOperatorConfigV21(
            nstate=2,
            nq=1,
            mass=28.0,
            seed=21430,
            base_scale=0.025,
            derivative_scale=0.008,
        )
    )


def test_explicit_zero_soc_path_is_exact_for_H_K_D_mass_and_overlaps():
    provenance = _provenance()
    base = _base()
    rehearsal = ZeroSOCRehearsalProviderV214(_base(), provenance)
    report = audit_zero_soc_equivalence_v214(
        base,
        rehearsal,
        [np.asarray([-0.3]), np.asarray([0.0]), np.asarray([0.4])],
        tolerance=0.0,
    )

    assert report.passed
    assert report.maximum_H_error == 0.0
    assert report.maximum_K_error == 0.0
    assert report.maximum_D_error == 0.0
    assert report.maximum_mass_error == 0.0
    assert report.maximum_overlap_error == 0.0
    require_provider_differential_contract_v214(
        rehearsal, np.asarray([0.17]), provenance
    )


def test_zero_soc_rehearsal_and_spin_free_dynamics_are_identical():
    provenance = _provenance()
    settings = SelfConsistentBlockSettingsV214(
        use_dense_reference=True,
        corrector_iterations=3,
        momentum_tolerance=1.0e-12,
    )
    base = run_self_consistent_block_dynamics_v214(
        ContractedElectronicOperatorProviderV213(_base(), provenance),
        provenance,
        initial_basis=_basis(),
        C0=_coefficients(),
        dt=0.002,
        steps=6,
        store_every=2,
        settings=settings,
    )
    rehearsal = run_self_consistent_block_dynamics_v214(
        ZeroSOCRehearsalProviderV214(_base(), provenance),
        provenance,
        initial_basis=_basis(),
        C0=_coefficients(),
        dt=0.002,
        steps=6,
        store_every=2,
        settings=settings,
    )

    assert max(
        np.linalg.norm(left.q - right.q)
        for left, right in zip(base["final_basis"], rehearsal["final_basis"])
    ) == 0.0
    assert max(
        np.linalg.norm(left.p - right.p)
        for left, right in zip(base["final_basis"], rehearsal["final_basis"])
    ) == 0.0
    assert _phase_aligned_metric_error(
        base["final_coefficients"],
        rehearsal["final_coefficients"],
        base["final_S"],
    ) < 1.0e-15
    assert rehearsal["checkpoint"].provider_fingerprint == provenance.fingerprint()


def test_zero_soc_rehearsal_refuses_soc_enabled_provenance():
    with pytest.raises(ValueError, match="explicitly disabled SOC"):
        ZeroSOCRehearsalProviderV214(_base(), _provenance(soc_enabled=True))
