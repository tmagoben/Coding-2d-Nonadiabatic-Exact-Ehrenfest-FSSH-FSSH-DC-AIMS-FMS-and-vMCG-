"""Explicit zero-SOC provider path used to rehearse the future decomposition.

No physical spin-orbit matrix is introduced.  The wrapper composes the total
operators from spin-free H/K and exact complex zero arrays through the frozen
v0.21.3 composition contract, then preserves the base provider's overlap path.
"""

from dataclasses import asdict, dataclass
import numpy as np

from .electronic_contract_v213 import compose_electronic_operator_v213
from .electronic_operator_v21 import ElectronicOperatorPointV21, ElectronicOperatorSnapshotV21


class ZeroSOCRehearsalProviderV214:
    def __init__(self, spin_free_provider, provenance):
        self.spin_free_provider = spin_free_provider
        self.provenance = provenance.validate()
        if self.provenance.soc_enabled or self.provenance.soc_method != "none":
            raise ValueError("the v0.21.4 rehearsal requires explicitly disabled SOC.")
        self.calls = 0

    def evaluate_snapshot(self, q):
        self.calls += 1
        base = self.spin_free_provider.evaluate_snapshot(q).validate()
        point = compose_electronic_operator_v213(
            q=base.point.q,
            H_spin_free=base.point.H,
            dH_spin_free_dq=base.point.hamiltonian_derivative_operator_q,
            H_soc=np.zeros_like(base.point.H, dtype=complex),
            dH_soc_dq=np.zeros_like(
                base.point.hamiltonian_derivative_operator_q, dtype=complex
            ),
            connection_q=base.point.connection_q,
            mass_matrix_q_au=base.point.mass_matrix_q_au,
            provenance=self.provenance,
        )
        point = ElectronicOperatorPointV21(
            q=point.q.copy(),
            H=point.H.copy(),
            dH_dq=point.dH_dq.copy(),
            connection_q=point.connection_q.copy(),
            mass_matrix_q_au=point.mass_matrix_q_au.copy(),
            metadata={
                **dict(base.point.metadata),
                **dict(point.metadata),
                "v214_zero_soc_rehearsal": True,
            },
        ).validate()
        return ElectronicOperatorSnapshotV21(
            point=point,
            state_vectors=(
                None
                if base.state_vectors is None
                else np.asarray(base.state_vectors, dtype=complex).copy()
            ),
            wavefunction_snapshot=base.wavefunction_snapshot,
            parent_snapshot=base,
            frame_from_parent=np.eye(point.nstate, dtype=complex),
            metadata={
                **dict(base.metadata),
                "provider": "ZeroSOCRehearsalProviderV214",
                "zero_soc": True,
                "provenance_fingerprint": self.provenance.fingerprint(),
            },
        ).validate()

    def evaluate(self, q):
        return self.evaluate_snapshot(q).point

    def snapshot_overlap(self, left, right):
        return self.spin_free_provider.snapshot_overlap(
            left.parent_snapshot, right.parent_snapshot
        )

    def diagnostics_dict(self):
        base = (
            self.spin_free_provider.diagnostics_dict()
            if hasattr(self.spin_free_provider, "diagnostics_dict")
            else {}
        )
        return {
            "provider": "ZeroSOCRehearsalProviderV214",
            "calls": int(self.calls),
            "physical_soc": False,
            "provenance_fingerprint": self.provenance.fingerprint(),
            "base": base,
        }


@dataclass(frozen=True)
class ZeroSOCEquivalenceReportV214:
    geometries: int
    maximum_H_error: float
    maximum_K_error: float
    maximum_D_error: float
    maximum_mass_error: float
    maximum_overlap_error: float
    provenance_fingerprint: str
    passed: bool
    tolerance: float

    def as_dict(self):
        return asdict(self)


def audit_zero_soc_equivalence_v214(
    spin_free_provider,
    rehearsal_provider,
    geometries,
    *,
    tolerance=0.0,
):
    tolerance = float(tolerance)
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("zero-SOC equivalence tolerance must be finite and nonnegative.")
    geometries = tuple(np.asarray(q, dtype=float) for q in geometries)
    if not geometries:
        raise ValueError("zero-SOC equivalence requires at least one geometry.")
    base_snapshots = []
    rehearsal_snapshots = []
    H_errors = []
    K_errors = []
    D_errors = []
    mass_errors = []
    for q in geometries:
        if q.ndim != 1 or not np.all(np.isfinite(q)):
            raise ValueError("zero-SOC audit geometries must be finite vectors.")
        base = spin_free_provider.evaluate_snapshot(q).validate()
        rehearsal = rehearsal_provider.evaluate_snapshot(q).validate()
        base_snapshots.append(base)
        rehearsal_snapshots.append(rehearsal)
        H_errors.append(float(np.max(np.abs(base.point.H - rehearsal.point.H))))
        K_errors.append(
            float(np.max(np.abs(base.point.dH_dq - rehearsal.point.dH_dq)))
        )
        D_errors.append(
            float(
                np.max(
                    np.abs(base.point.connection_q - rehearsal.point.connection_q)
                )
            )
        )
        mass_errors.append(
            float(
                np.max(
                    np.abs(
                        base.point.mass_matrix_q_au
                        - rehearsal.point.mass_matrix_q_au
                    )
                )
            )
        )

    overlap_errors = []
    for index in range(len(geometries) - 1):
        base_overlap = np.asarray(
            spin_free_provider.snapshot_overlap(
                base_snapshots[index], base_snapshots[index + 1]
            ),
            dtype=complex,
        )
        rehearsal_overlap = np.asarray(
            rehearsal_provider.snapshot_overlap(
                rehearsal_snapshots[index], rehearsal_snapshots[index + 1]
            ),
            dtype=complex,
        )
        if base_overlap.shape != rehearsal_overlap.shape:
            raise ValueError("zero-SOC overlap paths returned incompatible shapes.")
        overlap_errors.append(
            float(np.max(np.abs(base_overlap - rehearsal_overlap)))
        )
    maximum_H = max(H_errors, default=0.0)
    maximum_K = max(K_errors, default=0.0)
    maximum_D = max(D_errors, default=0.0)
    maximum_mass = max(mass_errors, default=0.0)
    maximum_overlap = max(overlap_errors, default=0.0)
    passed = max(
        maximum_H,
        maximum_K,
        maximum_D,
        maximum_mass,
        maximum_overlap,
    ) <= tolerance
    return ZeroSOCEquivalenceReportV214(
        geometries=len(geometries),
        maximum_H_error=maximum_H,
        maximum_K_error=maximum_K,
        maximum_D_error=maximum_D,
        maximum_mass_error=maximum_mass,
        maximum_overlap_error=maximum_overlap,
        provenance_fingerprint=rehearsal_provider.provenance.fingerprint(),
        passed=bool(passed),
        tolerance=tolerance,
    )
