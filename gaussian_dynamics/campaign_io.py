from pathlib import Path
import json
import numpy as np
from dataclasses import is_dataclass, asdict


def _convert(obj):
    if is_dataclass(obj):
        return {k: _convert(v) for k, v in asdict(obj).items()}
    if isinstance(obj, np.ndarray):
        return _convert(obj.tolist())
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, complex):
        return [float(obj.real), float(obj.imag)]
    if isinstance(obj, dict):
        return {str(k): _convert(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_convert(v) for v in obj]
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    return str(obj)


def save_campaign_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_convert(payload), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return path
