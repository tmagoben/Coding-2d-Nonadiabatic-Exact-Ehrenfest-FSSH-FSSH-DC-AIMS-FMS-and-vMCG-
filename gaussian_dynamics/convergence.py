from dataclasses import dataclass
import json
from pathlib import Path
import numpy as np


@dataclass
class RefinementResult:
    resolutions: np.ndarray
    values: np.ndarray
    successive_errors: np.ndarray
    observed_orders: np.ndarray


def phase_aligned_state_error(state, reference):
    a = np.asarray(state, dtype=complex).reshape(-1)
    b = np.asarray(reference, dtype=complex).reshape(-1)
    inner = np.vdot(b, a)
    if abs(inner) > 0:
        a = a * np.conj(inner) / abs(inner)
    return float(np.linalg.norm(a - b))


def observed_order(e_coarse, e_fine, refinement_ratio=2.0):
    e_coarse = float(e_coarse)
    e_fine = float(e_fine)
    if e_coarse <= 0.0 or e_fine <= 0.0:
        return np.nan
    return float(np.log(e_coarse / e_fine) / np.log(refinement_ratio))


def scalar_refinement_study(resolutions, values):
    r = np.asarray(resolutions, dtype=float)
    v = np.asarray(values, dtype=float)
    if r.ndim != 1 or v.ndim != 1 or len(r) != len(v):
        raise ValueError("resolutions and values must be equal-length vectors")
    if len(r) < 3:
        raise ValueError("at least three refinement levels are required")

    errors = np.abs(np.diff(v))
    orders = []
    for i in range(len(errors)-1):
        ratio = r[i] / r[i+1]
        orders.append(observed_order(errors[i], errors[i+1], ratio))

    return RefinementResult(
        resolutions=r,
        values=v,
        successive_errors=errors,
        observed_orders=np.asarray(orders),
    )


def vector_refinement_study(resolutions, values):
    r = np.asarray(resolutions, dtype=float)
    values = [np.asarray(v, dtype=float) for v in values]
    if len(r) != len(values) or len(r) < 3:
        raise ValueError("at least three matching refinement levels are required")

    errors = np.asarray([
        np.linalg.norm(values[i] - values[i+1])
        for i in range(len(values)-1)
    ])
    orders = []
    for i in range(len(errors)-1):
        ratio = r[i] / r[i+1]
        orders.append(observed_order(errors[i], errors[i+1], ratio))

    return RefinementResult(
        resolutions=r,
        values=np.asarray(values),
        successive_errors=errors,
        observed_orders=np.asarray(orders),
    )


def converged(last_error, atol, rtol=0.0, scale=1.0):
    return bool(float(last_error) <= float(atol) + float(rtol) * abs(float(scale)))


def save_convergence_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    def convert(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.floating, np.integer)):
            return obj.item()
        if isinstance(obj, complex):
            return [obj.real, obj.imag]
        raise TypeError(type(obj).__name__)

    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=convert), encoding="utf-8")
    return path
