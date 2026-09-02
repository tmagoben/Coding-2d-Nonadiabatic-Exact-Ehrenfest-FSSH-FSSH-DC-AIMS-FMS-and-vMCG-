"""v0.21.3 self-consistent propagation with density-matrix guidance."""

from dataclasses import dataclass, asdict
import numpy as np

from .block_sparse_molecular_v21 import BlockSparseSettingsV21
from .density_guidance_v213 import (
    BlockDensityMatrixGuidanceV213,
    DensityMatrixGuidanceSettingsV213,
)
from .self_consistent_block_v212 import (
    MeanFieldGuidanceSettingsV212,
    SelfConsistentBlockSettingsV212,
    run_self_consistent_block_dynamics_v212,
)


@dataclass(frozen=True)
class SelfConsistentBlockSettingsV213:
    graph: BlockSparseSettingsV21 = BlockSparseSettingsV21()
    guidance: DensityMatrixGuidanceSettingsV213 = DensityMatrixGuidanceSettingsV213()
    use_dense_reference: bool = True
    corrector_iterations: int = 2
    momentum_tolerance: float = 1.0e-10

    def validate(self):
        self.graph.validate()
        self.guidance.validate()
        if (
            self.corrector_iterations < 1
            or int(self.corrector_iterations) != self.corrector_iterations
        ):
            raise ValueError("corrector_iterations must be >=1.")
        if not np.isfinite(self.momentum_tolerance) or self.momentum_tolerance <= 0.0:
            raise ValueError("momentum_tolerance must be positive.")
        return self

    def legacy_control_settings(self):
        """Reuse the validated v0.21.2 integrator while replacing its guidance layer."""
        return SelfConsistentBlockSettingsV212(
            graph=self.graph,
            guidance=MeanFieldGuidanceSettingsV212(
                minimum_local_amplitude=self.guidance.minimum_local_amplitude,
                low_amplitude_policy="zero_force",
            ),
            use_dense_reference=self.use_dense_reference,
            corrector_iterations=self.corrector_iterations,
            momentum_tolerance=self.momentum_tolerance,
        )


def run_self_consistent_block_dynamics_v213(
    initial_basis,
    C0,
    provider,
    *,
    dt=0.002,
    steps=20,
    settings=SelfConsistentBlockSettingsV213(),
    store_every=5,
    adaptation_policy=None,
):
    settings = settings.validate()
    guidance = BlockDensityMatrixGuidanceV213(settings.guidance)
    output = run_self_consistent_block_dynamics_v212(
        initial_basis,
        C0,
        provider,
        dt=dt,
        steps=steps,
        settings=settings.legacy_control_settings(),
        store_every=store_every,
        adaptation_policy=adaptation_policy,
        guidance_engine=guidance,
    )
    output["settings"] = {
        "dt": float(dt),
        "steps": int(steps),
        "control": asdict(settings),
        "guidance_contract": "transported electronic density matrix",
    }
    output["release_path"] = "v0.21.3"
    return output
