from dataclasses import dataclass, field
from itertools import permutations
import numpy as np


@dataclass
class StateTrackingResult:
    """Result of maximum-overlap state assignment.

    permutation[i] is the raw-current state assigned to previous/tracked state i.
    phase_factors[i] multiplies that raw-current ket so the assigned overlap is
    positive in the selected gauge.
    """
    permutation: np.ndarray
    phase_factors: np.ndarray
    assigned_overlaps: np.ndarray
    best_score: float
    second_best_score: float
    score_margin: float
    ambiguous: bool
    reasons: tuple = field(default_factory=tuple)

    def as_metadata(self):
        return {
            "permutation_tracked_to_raw": self.permutation.astype(int).tolist(),
            "phase_factors": [
                [float(np.real(z)), float(np.imag(z))]
                for z in self.phase_factors
            ],
            "assigned_overlap_magnitudes": np.abs(self.assigned_overlaps).tolist(),
            "assigned_overlaps": [
                [float(np.real(z)), float(np.imag(z))]
                for z in self.assigned_overlaps
            ],
            "best_score": float(self.best_score),
            "second_best_score": float(self.second_best_score),
            "score_margin": float(self.score_margin),
            "ambiguous": bool(self.ambiguous),
            "ambiguity_reasons": list(self.reasons),
        }


def _assignment_scores(overlap):
    overlap = np.asarray(overlap, dtype=complex)
    if overlap.ndim != 2 or overlap.shape[0] != overlap.shape[1]:
        raise ValueError("State-overlap matrix must be square.")

    n = overlap.shape[0]
    scored = []

    for perm in permutations(range(n)):
        idx = np.asarray(perm, dtype=int)
        assigned = overlap[np.arange(n), idx]
        # Squared magnitudes give a natural probability-like continuity score.
        score = float(np.sum(np.abs(assigned) ** 2))
        scored.append((score, idx, assigned.copy()))

    scored.sort(key=lambda item: item[0], reverse=True)
    return scored


def maximum_overlap_assignment(
    overlap,
    minimum_overlap=0.50,
    minimum_score_margin=0.05,
    real_gauge=True,
    imaginary_tolerance=1e-8,
):
    """Assign current roots to previous roots by maximum many-electron overlap.

    Parameters
    ----------
    overlap
        O_ij = <Psi_i(previous) | Psi_j(current)>.
    minimum_overlap
        Minimum magnitude accepted for every assigned pair.
    minimum_score_margin
        Required gap between the best and second-best global permutation score.
        This is deliberately conservative near ambiguous degeneracies.
    real_gauge
        If True, restrict phase correction to +/-1 and reject assigned overlaps with
        appreciable imaginary components. This matches the real RHF/ROHF SA-CASSCF
        backend used in v0.6.
    """
    overlap = np.asarray(overlap, dtype=complex)
    scored = _assignment_scores(overlap)

    best_score, perm, assigned = scored[0]
    second_score = scored[1][0] if len(scored) > 1 else -np.inf
    margin = np.inf if len(scored) == 1 else best_score - second_score

    reasons = []

    mags = np.abs(assigned)
    if np.min(mags) < minimum_overlap:
        reasons.append(
            f"minimum assigned overlap {np.min(mags):.6g} "
            f"is below threshold {minimum_overlap:.6g}"
        )

    if np.isfinite(margin) and margin < minimum_score_margin:
        reasons.append(
            f"assignment score margin {margin:.6g} "
            f"is below threshold {minimum_score_margin:.6g}"
        )

    phases = np.ones(len(assigned), dtype=complex)

    for i, z in enumerate(assigned):
        if abs(z) < 1e-15:
            phases[i] = 1.0
            continue

        if real_gauge:
            if abs(np.imag(z)) > imaginary_tolerance:
                reasons.append(
                    f"assigned overlap for tracked state {i} has imaginary part "
                    f"{np.imag(z):.6g}, incompatible with real-gauge tracking"
                )
            phases[i] = 1.0 if np.real(z) >= 0.0 else -1.0
        else:
            # If |current'> = phase |current>, then
            # <previous|current'> = phase <previous|current>.
            phases[i] = np.conj(z) / abs(z)

    # The overlap *after* the phase correction is useful for diagnostics.
    assigned_after = assigned * phases

    return StateTrackingResult(
        permutation=perm,
        phase_factors=phases,
        assigned_overlaps=assigned_after,
        best_score=best_score,
        second_best_score=second_score,
        score_margin=margin,
        ambiguous=bool(reasons),
        reasons=tuple(reasons),
    )


def transform_state_properties(
    energies,
    gradients,
    nac,
    tracking_result,
):
    """Reorder and gauge-transform state-resolved quantities.

    If |phi_i'> = p_i |phi_raw[perm_i]>, then

        d'_ij = p_i^* p_j d_raw[perm_i, perm_j].

    Energies and diagonal state gradients are phase invariant and only reorder.
    """
    energies = np.asarray(energies)
    gradients = np.asarray(gradients)
    nac = np.asarray(nac)

    perm = np.asarray(tracking_result.permutation, dtype=int)
    phase = np.asarray(tracking_result.phase_factors, dtype=complex)

    if energies.shape[0] != len(perm):
        raise ValueError("Energy array is inconsistent with tracking permutation.")
    if gradients.shape[0] != len(perm):
        raise ValueError("Gradient array is inconsistent with tracking permutation.")
    if nac.shape[0] != len(perm) or nac.shape[1] != len(perm):
        raise ValueError("NAC array is inconsistent with tracking permutation.")

    E = energies[perm].copy()
    G = gradients[perm].copy()

    D = nac[np.ix_(perm, perm)].astype(complex, copy=True)

    phase_matrix = np.conj(phase)[:, None] * phase[None, :]
    extra_axes = (1,) * (D.ndim - 2)
    D = D * phase_matrix.reshape(phase_matrix.shape + extra_axes)

    if np.isrealobj(nac):
        if np.max(np.abs(np.imag(D))) > 1e-10:
            raise ValueError(
                "Complex gauge transformation generated complex NACs from a real "
                "backend; use a complex-valued data contract for that case."
            )
        D = np.real(D)

    return E, G, D


def reorder_and_phase_ci_roots(ci_roots, tracking_result):
    """Return CI roots in tracked order with the selected ket phase."""
    roots = list(ci_roots)
    perm = tracking_result.permutation
    phases = tracking_result.phase_factors

    out = []
    for i, raw_index in enumerate(perm):
        out.append(np.asarray(roots[int(raw_index)], dtype=complex) * phases[i])
    return out


def energy_degeneracy_clusters(energies, tolerance=1e-4):
    """Group adjacent states whose energy gaps are <= tolerance.

    Energies are assumed to be in whatever state order the caller wants diagnosed.
    """
    E = np.asarray(energies, dtype=float)
    if E.ndim != 1:
        raise ValueError("energies must be one-dimensional.")

    clusters = []
    current = [0]

    for i in range(1, len(E)):
        if abs(E[i] - E[i - 1]) <= tolerance:
            current.append(i)
        else:
            clusters.append(tuple(current))
            current = [i]

    clusters.append(tuple(current))
    return tuple(clusters)


def subspace_overlap_singular_values(overlap, previous_indices, current_indices):
    """Principal-overlap singular values for selected previous/current subspaces."""
    O = np.asarray(overlap, dtype=complex)
    block = O[np.ix_(list(previous_indices), list(current_indices))]
    return np.linalg.svd(block, compute_uv=False)
