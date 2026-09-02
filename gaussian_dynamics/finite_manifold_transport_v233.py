"""Certified separation of physical overlaps and unitary electronic transport.

For finite retained manifolds the raw cross-geometry matrix

    O_lr[i,j] = <Phi_i(q_l) | Phi_j(q_r)>

is a contraction and need not be unitary.  Coefficients and operators must be
transported with its unitary polar factor, never with the contraction itself.
v0.23.3 centralizes that distinction and adds an independent trajectory-quality
policy on top of the necessary Hilbert-space consistency checks.
"""

from dataclasses import asdict, dataclass
import math

import numpy as np


OVERLAP_CONTRACT_ID_V233 = "finite-orthonormal-manifold-contraction-v2"
TRANSPORT_CONTRACT_ID_V233 = "right-to-left-unitary-polar-factor-v1"


@dataclass(frozen=True)
class FiniteManifoldOverlapPolicyV233:
    """Physical and numerical requirements for one cross-geometry overlap."""

    contraction_tolerance: float = 1.0e-10
    minimum_retained_singular_value: float = 0.5
    maximum_condition_number: float = 1.0e6
    maximum_principal_angle_radians: float = math.pi / 3.0
    transport_unitarity_tolerance: float = 1.0e-10
    polar_hermiticity_tolerance: float = 1.0e-10

    def validate(self):
        finite_nonnegative = (
            "contraction_tolerance",
            "transport_unitarity_tolerance",
            "polar_hermiticity_tolerance",
        )
        for name in finite_nonnegative:
            value = float(getattr(self, name))
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and nonnegative.")
        minimum = float(self.minimum_retained_singular_value)
        if not np.isfinite(minimum) or not 0.0 < minimum <= 1.0:
            raise ValueError(
                "minimum_retained_singular_value must lie in (0,1]."
            )
        condition = float(self.maximum_condition_number)
        if not np.isfinite(condition) or condition < 1.0:
            raise ValueError("maximum_condition_number must be finite and >= 1.")
        angle = float(self.maximum_principal_angle_radians)
        if not np.isfinite(angle) or not 0.0 <= angle < math.pi / 2.0:
            raise ValueError(
                "maximum_principal_angle_radians must lie in [0,pi/2)."
            )
        return self

    def as_dict(self):
        return asdict(self)


CONSUMER_OVERLAP_POLICY_V233 = FiniteManifoldOverlapPolicyV233(
    minimum_retained_singular_value=1.0e-8,
    maximum_condition_number=1.0e8,
    maximum_principal_angle_radians=math.acos(1.0e-8),
)


@dataclass(frozen=True)
class FiniteManifoldTransportV233:
    """Raw overlap, its right-to-left unitary polar transport, and diagnostics."""

    overlap: np.ndarray
    right_to_left_transport: np.ndarray
    singular_values: np.ndarray
    minimum_singular_value: float
    maximum_singular_value: float
    condition_number: float
    maximum_principal_angle_radians: float
    contraction_excess: float
    transport_unitarity_residual: float
    polar_hermiticity_residual: float
    polar_minimum_eigenvalue: float
    physically_consistent: bool
    trajectory_ready: bool
    failed_quality_checks: tuple

    def as_dict(self):
        return {
            "overlap_contract": OVERLAP_CONTRACT_ID_V233,
            "transport_contract": TRANSPORT_CONTRACT_ID_V233,
            "singular_values": np.asarray(self.singular_values, dtype=float).tolist(),
            "minimum_singular_value": float(self.minimum_singular_value),
            "maximum_singular_value": float(self.maximum_singular_value),
            "condition_number": float(self.condition_number),
            "maximum_principal_angle_radians": float(
                self.maximum_principal_angle_radians
            ),
            "contraction_excess": float(self.contraction_excess),
            "transport_unitarity_residual": float(
                self.transport_unitarity_residual
            ),
            "polar_hermiticity_residual": float(self.polar_hermiticity_residual),
            "polar_minimum_eigenvalue": float(self.polar_minimum_eigenvalue),
            "physically_consistent": bool(self.physically_consistent),
            "trajectory_ready": bool(self.trajectory_ready),
            "failed_quality_checks": list(self.failed_quality_checks),
        }


def _square_finite_overlap_v233(overlap):
    matrix = np.asarray(overlap, dtype=complex)
    if matrix.ndim != 2 or matrix.shape[0] < 1 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("finite-manifold overlap must be a nonempty square matrix.")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("finite-manifold overlap contains non-finite data.")
    return matrix


def analyze_finite_manifold_overlap_v233(
    overlap,
    *,
    policy=FiniteManifoldOverlapPolicyV233(),
):
    """Analyze a physical overlap without silently using it as transport."""
    policy = policy.validate()
    matrix = _square_finite_overlap_v233(overlap)
    left_vectors, singular_values, right_vectors_h = np.linalg.svd(
        matrix, full_matrices=False
    )
    transport = left_vectors @ right_vectors_h
    identity = np.eye(matrix.shape[0], dtype=complex)
    transport_unitarity_residual = float(
        max(
            np.linalg.norm(transport.conj().T @ transport - identity, ord="fro"),
            np.linalg.norm(transport @ transport.conj().T - identity, ord="fro"),
        )
    )
    positive_factor = transport.conj().T @ matrix
    polar_hermiticity_residual = float(
        np.linalg.norm(positive_factor - positive_factor.conj().T, ord="fro")
    )
    polar_eigenvalues = np.linalg.eigvalsh(
        0.5 * (positive_factor + positive_factor.conj().T)
    )
    minimum = float(np.min(singular_values))
    maximum = float(np.max(singular_values))
    condition = float(np.inf if minimum == 0.0 else maximum / minimum)
    principal_angle = float(math.acos(float(np.clip(minimum, 0.0, 1.0))))
    contraction_excess = float(max(0.0, maximum - 1.0))

    physical_failures = []
    if contraction_excess > policy.contraction_tolerance:
        physical_failures.append("spectral_expansion")
    if transport_unitarity_residual > policy.transport_unitarity_tolerance:
        physical_failures.append("transport_nonunitarity")
    if polar_hermiticity_residual > policy.polar_hermiticity_tolerance:
        physical_failures.append("polar_factor_nonhermitian")
    if float(np.min(polar_eigenvalues)) < -policy.polar_hermiticity_tolerance:
        physical_failures.append("polar_factor_not_positive")

    quality_failures = []
    if minimum < policy.minimum_retained_singular_value:
        quality_failures.append("insufficient_manifold_retention")
    if condition > policy.maximum_condition_number:
        quality_failures.append("ill_conditioned_overlap")
    if principal_angle > policy.maximum_principal_angle_radians:
        quality_failures.append("principal_angle_too_large")

    physically_consistent = not physical_failures
    failures = tuple(physical_failures + quality_failures)
    return FiniteManifoldTransportV233(
        overlap=matrix.copy(),
        right_to_left_transport=transport,
        singular_values=np.asarray(singular_values, dtype=float),
        minimum_singular_value=minimum,
        maximum_singular_value=maximum,
        condition_number=condition,
        maximum_principal_angle_radians=principal_angle,
        contraction_excess=contraction_excess,
        transport_unitarity_residual=transport_unitarity_residual,
        polar_hermiticity_residual=polar_hermiticity_residual,
        polar_minimum_eigenvalue=float(np.min(polar_eigenvalues)),
        physically_consistent=physically_consistent,
        trajectory_ready=bool(physically_consistent and not quality_failures),
        failed_quality_checks=failures,
    )


def certified_transport_from_overlap_v233(
    overlap,
    *,
    policy=CONSUMER_OVERLAP_POLICY_V233,
):
    """Return the unitary polar transport only after all policy checks pass."""
    result = analyze_finite_manifold_overlap_v233(overlap, policy=policy)
    if not result.physically_consistent:
        raise ValueError(
            "finite-manifold overlap is physically inconsistent: "
            + ", ".join(result.failed_quality_checks)
        )
    if not result.trajectory_ready:
        raise ValueError(
            "finite-manifold overlap is not trajectory ready: "
            + ", ".join(result.failed_quality_checks)
        )
    return result


@dataclass(frozen=True)
class ReciprocalTransportPairV233:
    left_to_right_block: FiniteManifoldTransportV233
    right_to_left_block: FiniteManifoldTransportV233
    overlap_reciprocity_residual: float
    transport_reciprocity_residual: float

    def as_dict(self):
        return {
            "left_to_right_block": self.left_to_right_block.as_dict(),
            "right_to_left_block": self.right_to_left_block.as_dict(),
            "overlap_reciprocity_residual": float(
                self.overlap_reciprocity_residual
            ),
            "transport_reciprocity_residual": float(
                self.transport_reciprocity_residual
            ),
        }


def certify_reciprocal_transport_pair_v233(
    overlap_left_right,
    overlap_right_left,
    *,
    policy=FiniteManifoldOverlapPolicyV233(),
    reciprocity_tolerance=1.0e-10,
):
    """Certify raw-overlap and polar-transport adjoint reciprocity."""
    tolerance = float(reciprocity_tolerance)
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("reciprocity_tolerance must be finite and nonnegative.")
    forward = certified_transport_from_overlap_v233(
        overlap_left_right, policy=policy
    )
    reverse = certified_transport_from_overlap_v233(
        overlap_right_left, policy=policy
    )
    overlap_residual = float(
        np.linalg.norm(
            forward.overlap - reverse.overlap.conj().T, ord="fro"
        )
    )
    transport_residual = float(
        np.linalg.norm(
            forward.right_to_left_transport
            - reverse.right_to_left_transport.conj().T,
            ord="fro",
        )
    )
    if overlap_residual > tolerance:
        raise ValueError("finite-manifold overlaps violate adjoint reciprocity.")
    if transport_residual > tolerance:
        raise ValueError("finite-manifold transports violate adjoint reciprocity.")
    return ReciprocalTransportPairV233(
        left_to_right_block=forward,
        right_to_left_block=reverse,
        overlap_reciprocity_residual=overlap_residual,
        transport_reciprocity_residual=transport_residual,
    )
