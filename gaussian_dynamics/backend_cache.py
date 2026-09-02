from pathlib import Path
import hashlib
import json
import numpy as np

from .molecular_backend import GeneralizedElectronicStructurePoint


class DiskCachedGeneralizedProvider:
    """NPZ/JSON cache for generalized-coordinate electronic points."""

    def __init__(self, provider, directory, namespace="default", decimals=12):
        self.provider = provider
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.namespace = str(namespace)
        self.decimals = int(decimals)
        self.hits = 0
        self.misses = 0

    def _key(self, q):
        payload = {
            "namespace": self.namespace,
            "q": np.asarray(q, dtype=float).round(self.decimals).tolist(),
        }
        raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
        return hashlib.sha256(raw).hexdigest()

    def evaluate(self, q):
        q = np.asarray(q, dtype=float)
        key = self._key(q)
        npz_path = self.directory / f"{key}.npz"
        json_path = self.directory / f"{key}.json"

        if npz_path.exists() and json_path.exists():
            self.hits += 1
            data = np.load(npz_path)
            metadata = json.loads(json_path.read_text(encoding="utf-8"))

            return GeneralizedElectronicStructurePoint(
                q=data["q"],
                energies=data["energies"],
                gradients_q=data["gradients_q"],
                nac_q=data["nac_q"],
                mass_matrix_q_au=data["mass_matrix_q_au"],
                metadata=metadata,
            ).validate()

        self.misses += 1
        point = self.provider.evaluate(q).validate()

        np.savez_compressed(
            npz_path,
            q=point.q,
            energies=point.energies,
            gradients_q=point.gradients_q,
            nac_q=point.nac_q,
            mass_matrix_q_au=point.mass_matrix_q_au,
        )
        json_path.write_text(
            json.dumps(point.metadata, sort_keys=True, indent=2, default=str),
            encoding="utf-8",
        )

        return point
