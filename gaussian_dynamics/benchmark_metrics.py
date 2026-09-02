from dataclasses import dataclass, asdict
import numpy as np


def population_l2_error(candidate, reference):
    candidate = np.asarray(candidate, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if candidate.shape != reference.shape:
        raise ValueError("population vectors must have equal shape.")
    return float(np.linalg.norm(candidate - reference))


def population_l1_error(candidate, reference):
    candidate = np.asarray(candidate, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if candidate.shape != reference.shape:
        raise ValueError("population vectors must have equal shape.")
    return float(np.sum(np.abs(candidate-reference)))


def population_sum_error(populations):
    return float(abs(np.sum(np.asarray(populations, dtype=float)) - 1.0))


def norm_error(norm):
    return float(abs(float(norm) - 1.0))


@dataclass(frozen=True)
class ManagedRunMetrics:
    final_populations: np.ndarray
    final_norm: float
    max_norm_error: float
    max_condition_number: float
    max_basis_size: int
    spawn_count: int
    prune_count: int
    total_pruning_loss: float
    max_spa1_relative_correction: float
    population_l2_vs_reference: float | None = None

    def to_dict(self):
        out = asdict(self)
        out["final_populations"] = np.asarray(self.final_populations).tolist()
        return out


def summarize_managed_run(run, reference_populations=None):
    records = run["records"]
    if not records:
        raise ValueError("managed run has no records.")

    final = records[-1]
    final_pop = np.asarray(final["state_populations"], dtype=float)

    l2 = None
    if reference_populations is not None:
        l2 = population_l2_error(final_pop, reference_populations)

    spawn_count = sum(e.get("kind") == "spawn" for e in run.get("events", []))
    prune_events = [e for e in run.get("events", []) if e.get("kind") == "prune"]

    return ManagedRunMetrics(
        final_populations=final_pop,
        final_norm=float(final["norm"]),
        max_norm_error=max(norm_error(r["norm"]) for r in records),
        max_condition_number=max(float(r["condition_number"]) for r in records),
        max_basis_size=max(int(r["basis_size"]) for r in records),
        spawn_count=int(spawn_count),
        prune_count=len(prune_events),
        total_pruning_loss=float(sum(e.get("projection_loss", 0.0) for e in prune_events)),
        max_spa1_relative_correction=max(
            float(r.get("spa1_relative_correction", 0.0)) for r in records
        ),
        population_l2_vs_reference=l2,
    )
