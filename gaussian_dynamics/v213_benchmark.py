"""Canonical v0.21.3 SOC-contract-freeze acceptance campaign."""

from dataclasses import asdict, dataclass
import importlib.util
from pathlib import Path
import tempfile

import numpy as np

from .block_sparse_molecular_v21 import BlockMolecularTBFV21
from .complex_gauge_v21 import (
    GaugeTransformedOperatorProviderV21,
    PhaseMixingGaugeV21,
)
from .complex_operator_cache_v213 import FixedFrameComplexOperatorCacheV213
from .density_guidance_v213 import BlockDensityMatrixGuidanceV213
from .electronic_contract_v213 import (
    ElectronicModelSpaceV213,
    ElectronicOperatorProvenanceV213,
    ElectronicStateDescriptorV213,
    compose_electronic_operator_v213,
    hartree_to_wavenumber_v213,
    validate_electronic_contract_v213,
    wavenumber_to_hartree_v213,
)
from .electronic_operator_v21 import (
    ElectronicOperatorPointV21,
    ElectronicOperatorSnapshotV21,
)
from .gaussian_nd import gaussian_nd
from .initial_projection_v213 import project_grid_wavefunction_fixed_frame_v213
from .self_consistent_block_v212 import MeanFieldGuidanceSettingsV212
from .self_consistent_block_v213 import (
    SelfConsistentBlockSettingsV213,
    run_self_consistent_block_dynamics_v213,
)
from .synthetic_operator_provider_v21 import (
    SyntheticLinearOperatorConfigV21,
    SyntheticLinearOperatorProviderV21,
)
from .v212_benchmark import run_v0212_release_benchmark


@dataclass(frozen=True)
class V213AcceptanceThresholds:
    max_structural_residual: float = 1.0e-12
    max_zero_soc_composition_error: float = 1.0e-15
    max_unit_roundtrip_error: float = 5.0e-13
    max_degenerate_gauge_force_error: float = 5.0e-13
    max_retained_force_error: float = 5.0e-13
    max_unseeded_force: float = 1.0e-15
    min_projection_fidelity: float = 1.0 - 1.0e-12
    max_projection_relative_residual: float = 1.0e-11
    max_cache_roundtrip_error: float = 0.0
    min_cached_imaginary_signal: float = 1.0e-8
    max_integrated_norm_drift: float = 1.0e-10
    require_inherited_v0212: bool = True


def _spin_model_space_v213():
    return ElectronicModelSpaceV213(
        name="complete S0 plus T1 spin-component fixture",
        representation="fixed_spin_diabatic",
        states=(
            ElectronicStateDescriptorV213("S0", "S0", 1, "M=0", 0),
            ElectronicStateDescriptorV213("T1(M=-1)", "T1", 3, "M=-1", 0),
            ElectronicStateDescriptorV213("T1(M=0)", "T1", 3, "M=0", 0),
            ElectronicStateDescriptorV213("T1(M=+1)", "T1", 3, "M=+1", 0),
        ),
        complete_multiplets=True,
    ).validate()


def _provenance_v213(parameter=1.0):
    return ElectronicOperatorProvenanceV213(
        model_name="v0.21.3 zero-SOC contract fixture",
        model_version="1",
        model_space=_spin_model_space_v213(),
        spin_free_method="analytic linear fixture",
        soc_enabled=False,
        soc_method="none",
        derivative_method="analytic physical operator derivative",
        parameters={"linear_parameter": float(parameter)},
    ).validate()


class _ExactlyDegenerateProviderV213:
    def evaluate_snapshot(self, q):
        q = np.asarray(q, dtype=float)
        derivative = np.asarray([np.diag([1.0, -1.0])], dtype=complex)
        point = ElectronicOperatorPointV21(
            q=q.copy(),
            H=np.zeros((2, 2), dtype=complex),
            dH_dq=derivative,
            connection_q=np.zeros_like(derivative),
            mass_matrix_q_au=np.asarray([[20.0]]),
            metadata={"fixture": "exact electronic degeneracy"},
        ).validate()
        return ElectronicOperatorSnapshotV21(
            point=point, state_vectors=np.eye(2, dtype=complex)
        ).validate()

    @staticmethod
    def snapshot_overlap(left, right):
        return left.state_vectors.conj().T @ right.state_vectors


def _strict_invariant_campaign_v213():
    defective = ElectronicOperatorPointV21(
        q=np.asarray([0.0]),
        H=np.diag([1.0 + 1.0e-6j, 2.0]),
        dH_dq=np.zeros((1, 2, 2), dtype=complex),
        connection_q=np.zeros((1, 2, 2), dtype=complex),
        mass_matrix_q_au=np.eye(1),
    )
    historical_allclose_accepts = bool(
        np.allclose(defective.H, defective.H.conj().T, atol=1.0e-12)
    )
    strict_rejected = False
    rejection = ""
    try:
        defective.validate(atol=1.0e-12)
    except ValueError as exc:
        strict_rejected = True
        rejection = str(exc)
    return {
        "historical_default_allclose_accepts_fixture": historical_allclose_accepts,
        "strict_validator_rejected_fixture": strict_rejected,
        "strict_rejection": rejection,
    }


def _operator_contract_campaign_v213():
    provenance = _provenance_v213()
    H0 = np.diag([0.0, 0.02, 0.02, 0.02]).astype(complex)
    K0 = np.zeros((2, 4, 4), dtype=complex)
    K0[0] = np.diag([0.01, -0.02, -0.02, -0.02])
    K0[1, 0, 1] = 0.003 + 0.004j
    K0[1, 1, 0] = K0[1, 0, 1].conjugate()
    point = compose_electronic_operator_v213(
        q=np.asarray([0.1, -0.2]),
        H_spin_free=H0,
        dH_spin_free_dq=K0,
        connection_q=np.zeros_like(K0),
        mass_matrix_q_au=np.diag([1500.0, 1800.0]),
        provenance=provenance,
    )
    residuals = point.structural_residuals_v213()
    structural_maximum = max(
        residuals["H_hermiticity"],
        *residuals["dH_hermiticity"],
        *residuals["connection_antihermiticity"],
        residuals["mass_symmetry"],
    )
    values = np.asarray([0.0, 125.0, 1000.0])
    unit_roundtrip_error = float(
        np.max(
            np.abs(
                hartree_to_wavenumber_v213(wavenumber_to_hartree_v213(values))
                - values
            )
        )
    )

    incomplete_multiplet_rejected = False
    try:
        ElectronicModelSpaceV213(
            name="incomplete triplet fixture",
            representation="fixed_spin_diabatic",
            states=(
                ElectronicStateDescriptorV213("T(-1)", "T", 3, "M=-1"),
                ElectronicStateDescriptorV213("T(0)", "T", 3, "M=0"),
            ),
            complete_multiplets=True,
        ).validate()
    except ValueError:
        incomplete_multiplet_rejected = True

    nonzero_fixed_connection_rejected = False
    connection_point = ElectronicOperatorPointV21(
        q=np.asarray([0.0]),
        H=np.zeros((4, 4), dtype=complex),
        dH_dq=np.zeros((1, 4, 4), dtype=complex),
        connection_q=np.asarray([1.0j * np.eye(4)]),
        mass_matrix_q_au=np.eye(1),
    ).validate()
    try:
        validate_electronic_contract_v213(connection_point, provenance)
    except ValueError:
        nonzero_fixed_connection_rejected = True

    nonzero_soc_without_provenance_rejected = False
    try:
        compose_electronic_operator_v213(
            q=np.asarray([0.0]),
            H_spin_free=np.zeros((4, 4)),
            dH_spin_free_dq=np.zeros((1, 4, 4)),
            H_soc=np.diag([0.0, 1.0e-4, 0.0, 0.0]),
            connection_q=np.zeros((1, 4, 4)),
            mass_matrix_q_au=np.eye(1),
            provenance=provenance,
        )
    except ValueError:
        nonzero_soc_without_provenance_rejected = True

    return {
        "nstate": provenance.model_space.nstate,
        "complete_multiplets": provenance.model_space.complete_multiplets,
        "provenance_fingerprint": provenance.fingerprint(),
        "zero_soc_H_error": float(np.linalg.norm(point.H - H0, ord="fro")),
        "zero_soc_K_error": float(np.linalg.norm(point.dH_dq - K0)),
        "maximum_structural_residual": float(structural_maximum),
        "unit_roundtrip_error_cm_inverse": unit_roundtrip_error,
        "incomplete_multiplet_rejected": incomplete_multiplet_rejected,
        "nonzero_fixed_connection_rejected": nonzero_fixed_connection_rejected,
        "nonzero_soc_without_provenance_rejected": (
            nonzero_soc_without_provenance_rejected
        ),
    }


def _degenerate_guidance_campaign_v213():
    base_provider = _ExactlyDegenerateProviderV213()
    hadamard = np.asarray([[1.0, 1.0], [1.0, -1.0]], dtype=complex) / np.sqrt(2.0)
    gauge = PhaseMixingGaugeV21(
        U0=hadamard,
        phase_gradient=np.asarray([[0.23], [-0.17]]),
        phase_offset=np.asarray([0.31, -0.22]),
    )
    gauge_provider = GaugeTransformedOperatorProviderV21(base_provider, gauge)
    base_guide = BlockDensityMatrixGuidanceV213()
    gauge_guide = BlockDensityMatrixGuidanceV213()
    initial = BlockMolecularTBFV21(
        7, np.asarray([0.0]), np.asarray([0.0]), np.asarray([[1.2]])
    )
    c = np.asarray([1.0, 0.0], dtype=complex)
    c_gauge = gauge.matrix(initial.q).conj().T @ c
    current_base = base_guide.forces_and_masses([initial], c, base_provider, 2)[0]
    current_gauge = gauge_guide.forces_and_masses(
        [initial], c_gauge, gauge_provider, 2
    )[0]

    moved = BlockMolecularTBFV21(
        7, np.asarray([0.4]), np.asarray([0.0]), np.asarray([[1.2]])
    )
    zeros = np.zeros(2, dtype=complex)
    retained_base = base_guide.forces_and_masses(
        [moved], zeros, base_provider, 2
    )[0]
    retained_gauge = gauge_guide.forces_and_masses(
        [moved], zeros, gauge_provider, 2
    )[0]
    empty_guide = BlockDensityMatrixGuidanceV213()
    unseeded = empty_guide.forces_and_masses(
        [moved], zeros, base_provider, 2
    )[0]

    fallback_rejected = False
    try:
        MeanFieldGuidanceSettingsV212(
            low_amplitude_policy="lowest_eigenvector"
        ).validate()
    except ValueError:
        fallback_rejected = True
    integrated_provider = SyntheticLinearOperatorProviderV21(
        SyntheticLinearOperatorConfigV21(
            nstate=2, nq=1, seed=21309, mass=25.0
        )
    )
    integrated = run_self_consistent_block_dynamics_v213(
        [
            BlockMolecularTBFV21(
                11,
                np.asarray([-0.15]),
                np.asarray([0.0]),
                np.asarray([[1.2]]),
            )
        ],
        np.asarray([0.8 + 0.1j, -0.2 + 0.3j]),
        integrated_provider,
        dt=0.002,
        steps=3,
        store_every=1,
        settings=SelfConsistentBlockSettingsV213(
            corrector_iterations=2,
            momentum_tolerance=1.0e-12,
        ),
    )
    return {
        "current_base_force": current_base.tolist(),
        "current_gauge_force": current_gauge.tolist(),
        "retained_base_force": retained_base.tolist(),
        "retained_gauge_force": retained_gauge.tolist(),
        "current_gauge_error": float(np.max(np.abs(current_gauge - current_base))),
        "retained_gauge_error": float(np.max(np.abs(retained_gauge - retained_base))),
        "retained_expected_force_error": float(
            np.max(np.abs(retained_base - np.asarray([[-1.0]])))
        ),
        "unseeded_zero_force_norm": float(np.linalg.norm(unseeded)),
        "unsafe_eigenvector_fallback_rejected": fallback_rejected,
        "gauge_diagnostics": gauge_guide.diagnostics_dict(),
        "integrated_runner": {
            "release_path": integrated["release_path"],
            "maximum_norm_drift": integrated["maximum_norm_drift"],
            "guidance_trial_state_rollbacks": integrated[
                "guidance_trial_state_rollbacks"
            ],
            "coefficient_refreshes": integrated["guidance_diagnostics"][
                "coefficient_refreshes"
            ],
        },
    }


def _projection_campaign_v213():
    x = np.linspace(-9.0, 9.0, 4001)
    points = x[:, None]
    basis = [
        BlockMolecularTBFV21(
            0,
            np.asarray([0.35]),
            np.asarray([0.42]),
            np.asarray([[1.4]]),
        )
    ]
    electronic = np.asarray(
        [0.55 + 0.10j, -0.20 + 0.35j, 0.30 - 0.15j, -0.12 - 0.28j]
    )
    electronic /= np.linalg.norm(electronic)
    target = gaussian_nd(
        points, basis[0].q, basis[0].p, basis[0].A
    )[..., None] * electronic
    result = project_grid_wavefunction_fixed_frame_v213(
        target, points, x[1] - x[0], basis
    )
    return {
        "nstate": result.nstate,
        "nuclear_dimension": result.nuclear_dimension,
        "fidelity": result.fidelity,
        "relative_residual": result.relative_residual,
        "condition_number": result.condition_number,
        "coefficient_error": float(np.linalg.norm(result.coefficients - electronic)),
    }


def _cache_campaign_v213():
    model_space = ElectronicModelSpaceV213(
        name="three-state fixed complex cache fixture",
        representation="fixed_general",
        states=tuple(ElectronicStateDescriptorV213(f"state-{i}") for i in range(3)),
    )

    def provenance(parameter):
        return ElectronicOperatorProvenanceV213(
            model_name="synthetic complex cache fixture",
            model_version="1",
            model_space=model_space,
            spin_free_method="analytic linear fixture",
            parameters={"parameter": float(parameter)},
        )

    base = SyntheticLinearOperatorProviderV21(
        SyntheticLinearOperatorConfigV21(nstate=3, nq=2, seed=21301)
    )
    q = np.asarray([0.17, -0.28])
    with tempfile.TemporaryDirectory(prefix="v213-complex-cache-") as directory:
        first = FixedFrameComplexOperatorCacheV213(
            base, Path(directory), provenance(1.0), namespace="acceptance"
        )
        miss = first.evaluate_snapshot(q)
        hit = first.evaluate_snapshot(q)
        second = FixedFrameComplexOperatorCacheV213(
            base, Path(directory), provenance(2.0), namespace="acceptance"
        )
        second.evaluate_snapshot(q)
        files = sorted(path.name for path in Path(directory).iterdir())
    roundtrip_error = max(
        float(np.max(np.abs(hit.point.H - miss.point.H))),
        float(np.max(np.abs(hit.point.dH_dq - miss.point.dH_dq))),
        float(np.max(np.abs(hit.point.connection_q - miss.point.connection_q))),
    )
    imaginary_signal = max(
        float(np.max(np.abs(np.imag(hit.point.H)))),
        float(np.max(np.abs(np.imag(hit.point.dH_dq)))),
    )
    return {
        "roundtrip_error": roundtrip_error,
        "imaginary_signal": imaginary_signal,
        "base_provider_calls": int(base.calls),
        "first_cache_hits": int(first.hits),
        "first_cache_misses": int(first.misses),
        "second_cache_misses": int(second.misses),
        "fingerprints_differ": (
            first.provider_fingerprint != second.provider_fingerprint
        ),
        "entry_files": files,
        "entry_file_count": len(files),
    }


def run_v0213_release_benchmark():
    strict = _strict_invariant_campaign_v213()
    contract = _operator_contract_campaign_v213()
    guidance = _degenerate_guidance_campaign_v213()
    projection = _projection_campaign_v213()
    cache = _cache_campaign_v213()
    inherited = bool(run_v0212_release_benchmark()["acceptance"]["passed"])
    thresholds = V213AcceptanceThresholds()
    checks = {
        "strict_validator_closes_allclose_rtol_gap": (
            strict["historical_default_allclose_accepts_fixture"]
            and strict["strict_validator_rejected_fixture"]
        ),
        "complete_multiplet_model_space": (
            contract["nstate"] == 4
            and contract["complete_multiplets"]
            and contract["incomplete_multiplet_rejected"]
        ),
        "zero_soc_H_composition": (
            contract["zero_soc_H_error"]
            <= thresholds.max_zero_soc_composition_error
        ),
        "zero_soc_K_composition": (
            contract["zero_soc_K_error"]
            <= thresholds.max_zero_soc_composition_error
        ),
        "strict_operator_structure": (
            contract["maximum_structural_residual"]
            <= thresholds.max_structural_residual
        ),
        "explicit_internal_units": (
            contract["unit_roundtrip_error_cm_inverse"]
            <= thresholds.max_unit_roundtrip_error
        ),
        "fixed_frame_connection_contract": contract[
            "nonzero_fixed_connection_rejected"
        ],
        "soc_requires_explicit_provenance": contract[
            "nonzero_soc_without_provenance_rejected"
        ],
        "degenerate_current_force_covariance": (
            guidance["current_gauge_error"]
            <= thresholds.max_degenerate_gauge_force_error
        ),
        "degenerate_retained_force_covariance": (
            guidance["retained_gauge_error"]
            <= thresholds.max_degenerate_gauge_force_error
        ),
        "retained_density_force": (
            guidance["retained_expected_force_error"]
            <= thresholds.max_retained_force_error
        ),
        "unseeded_zero_block_force": (
            guidance["unseeded_zero_force_norm"]
            <= thresholds.max_unseeded_force
        ),
        "unsafe_eigenvector_fallback_retired": guidance[
            "unsafe_eigenvector_fallback_rejected"
        ],
        "transactional_corrector_guidance": (
            guidance["integrated_runner"]["release_path"] == "v0.21.3"
            and guidance["integrated_runner"]["guidance_trial_state_rollbacks"] > 0
            and guidance["integrated_runner"]["coefficient_refreshes"] > 0
            and guidance["integrated_runner"]["maximum_norm_drift"]
            <= thresholds.max_integrated_norm_drift
        ),
        "arbitrary_state_projection_fidelity": (
            projection["nstate"] == 4
            and projection["nuclear_dimension"] == 1
            and projection["fidelity"] >= thresholds.min_projection_fidelity
        ),
        "arbitrary_state_projection_residual": (
            projection["relative_residual"]
            <= thresholds.max_projection_relative_residual
        ),
        "complex_cache_roundtrip": (
            cache["roundtrip_error"] <= thresholds.max_cache_roundtrip_error
            and cache["first_cache_hits"] == 1
            and cache["first_cache_misses"] == 1
        ),
        "complex_cache_preserves_imaginary_data": (
            cache["imaginary_signal"] >= thresholds.min_cached_imaginary_signal
        ),
        "cache_provenance_separation": (
            cache["fingerprints_differ"]
            and cache["base_provider_calls"] == 2
            and cache["entry_file_count"] == 4
        ),
        "inherited_v0212": (
            inherited if thresholds.require_inherited_v0212 else True
        ),
    }
    return {
        "release": "v0.21.3",
        "theme": "SOC-contract freeze and degeneracy-safe pre-integration procedures",
        "strict_matrix_invariants": strict,
        "electronic_operator_contract": contract,
        "degeneracy_safe_density_guidance": guidance,
        "arbitrary_state_projection": projection,
        "fingerprinted_complex_cache": cache,
        "inherited_v0212_acceptance": inherited,
        "soc": {
            "physical_hamiltonian_introduced": False,
            "spin_free_mode_permanent": True,
            "first_physical_soc_target": "v0.22",
        },
        "pyscf": {
            "installed_in_build_environment": bool(
                importlib.util.find_spec("pyscf") is not None
            ),
            "runtime_validated": False,
            "note": "No ab-initio SOC claim is made in v0.21.3.",
        },
        "acceptance": {
            "passed": bool(all(checks.values())),
            "checks": checks,
            "thresholds": asdict(thresholds),
        },
    }
