from dataclasses import dataclass
from pathlib import Path
import json
import numpy as np

from .molecular_backend import GeneralizedElectronicStructurePoint


@dataclass
class TrackedScanResult:
    q: np.ndarray
    energies: np.ndarray
    gradients_q: np.ndarray
    nac_q: np.ndarray
    mass_matrix_q_au: np.ndarray
    metadata: list

    def validate(self):
        q = np.asarray(self.q, dtype=float)
        E = np.asarray(self.energies, dtype=float)
        G = np.asarray(self.gradients_q, dtype=float)
        D = np.asarray(self.nac_q, dtype=float)
        M = np.asarray(self.mass_matrix_q_au, dtype=float)

        if q.ndim != 2:
            raise ValueError("q must have shape (npoint,nq).")

        npoint, nq = q.shape
        if E.shape[0] != npoint:
            raise ValueError("energies do not match scan length.")
        nstate = E.shape[1]

        if G.shape != (npoint, nstate, nq):
            raise ValueError("gradients_q has incompatible shape.")
        if D.shape != (npoint, nstate, nstate, nq):
            raise ValueError("nac_q has incompatible shape.")
        if M.shape != (npoint, nq, nq):
            raise ValueError("mass_matrix_q_au has incompatible shape.")
        if len(self.metadata) != npoint:
            raise ValueError("metadata length does not match scan length.")

        self.q = q
        self.energies = E
        self.gradients_q = G
        self.nac_q = D
        self.mass_matrix_q_au = M
        return self


def run_tracked_scan(q_path, generalized_provider):
    """Evaluate a stateful tracked provider in the exact order supplied."""
    q_path = np.asarray(q_path, dtype=float)
    if q_path.ndim == 1:
        q_path = q_path[:, None]
    if q_path.ndim != 2:
        raise ValueError("q_path must be one- or two-dimensional.")

    points = [generalized_provider.evaluate(q) for q in q_path]

    return TrackedScanResult(
        q=q_path.copy(),
        energies=np.asarray([p.energies for p in points]),
        gradients_q=np.asarray([p.gradients_q for p in points]),
        nac_q=np.asarray([p.nac_q for p in points]),
        mass_matrix_q_au=np.asarray([p.mass_matrix_q_au for p in points]),
        metadata=[dict(p.metadata) for p in points],
    ).validate()


def save_tracked_scan(scan, prefix):
    scan = scan.validate()
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    npz_path = prefix.with_suffix(".npz")
    json_path = prefix.with_suffix(".json")

    np.savez_compressed(
        npz_path,
        q=scan.q,
        energies=scan.energies,
        gradients_q=scan.gradients_q,
        nac_q=scan.nac_q,
        mass_matrix_q_au=scan.mass_matrix_q_au,
    )
    json_path.write_text(
        json.dumps(scan.metadata, indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )
    return npz_path, json_path


class TrackedScan1DProvider:
    """Order-independent linear interpolation of an already tracked 1D scan."""

    def __init__(self, scan):
        scan = scan.validate()
        if scan.q.shape[1] != 1:
            raise ValueError("TrackedScan1DProvider requires nq=1.")

        q = scan.q[:, 0]
        if np.any(np.diff(q) <= 0.0):
            raise ValueError("Tracked 1D scan coordinates must be strictly increasing.")

        self.scan = scan
        self.grid = q

    def _interp_vector(self, arr, x):
        return np.array([
            np.interp(x, self.grid, arr[:, j])
            for j in range(arr.shape[1])
        ])

    def evaluate(self, q):
        q = np.asarray(q, dtype=float)
        if q.shape not in {(1,), ()}:
            raise ValueError("TrackedScan1DProvider expects one coordinate.")

        x = float(q.reshape(-1)[0])
        if x < self.grid[0] or x > self.grid[-1]:
            raise ValueError("Requested point lies outside the tracked scan domain.")

        E = self._interp_vector(self.scan.energies, x)

        nstate = self.scan.energies.shape[1]
        G = np.zeros((nstate, 1), dtype=float)
        D = np.zeros((nstate, nstate, 1), dtype=float)

        for i in range(nstate):
            G[i, 0] = np.interp(
                x, self.grid, self.scan.gradients_q[:, i, 0]
            )
            for j in range(nstate):
                D[i, j, 0] = np.interp(
                    x, self.grid, self.scan.nac_q[:, i, j, 0]
                )

        # Restore exact antisymmetry after interpolation.
        D[:, :, 0] = 0.5 * (D[:, :, 0] - D[:, :, 0].T)

        nq = 1
        M = np.zeros((nq, nq), dtype=float)
        M[0, 0] = np.interp(
            x,
            self.grid,
            self.scan.mass_matrix_q_au[:, 0, 0],
        )

        return GeneralizedElectronicStructurePoint(
            q=np.array([x]),
            energies=E,
            gradients_q=G,
            nac_q=D,
            mass_matrix_q_au=M,
            metadata={
                "provider": "tracked_scan_1d_linear_interpolation",
                "scan_domain": [float(self.grid[0]), float(self.grid[-1])],
            },
        ).validate()
