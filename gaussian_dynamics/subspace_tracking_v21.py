from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class SubspaceTrackingResultV21:
    transform: np.ndarray
    singular_values: np.ndarray
    aligned_overlap: np.ndarray
    antihermitian_residual: float
    minimum_singular_value: float
    principal_angle_max: float
    def as_dict(self):
        return {"singular_values": self.singular_values.tolist(), "antihermitian_residual": float(self.antihermitian_residual), "minimum_singular_value": float(self.minimum_singular_value), "principal_angle_max": float(self.principal_angle_max)}

def procrustes_subspace_alignment_v21(overlap):
    O = np.asarray(overlap, dtype=complex)
    U, s, Vh = np.linalg.svd(O, full_matrices=False)
    W = Vh.conj().T @ U.conj().T
    aligned = O @ W
    anti = float(np.linalg.norm(aligned - aligned.conj().T, ord="fro"))
    angle = float(np.max(np.arccos(np.clip(s, 0.0, 1.0))))
    return SubspaceTrackingResultV21(W, s, aligned, anti, float(np.min(s)), angle)

def transform_subspace_operator_v21(matrix, transform):
    A = np.asarray(matrix, dtype=complex); W = np.asarray(transform, dtype=complex)
    return W.conj().T @ A @ W
