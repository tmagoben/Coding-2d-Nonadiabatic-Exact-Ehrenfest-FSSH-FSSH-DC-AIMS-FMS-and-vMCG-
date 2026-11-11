"""Small analytic free-packet campaign for the v0.30.2 development integrator."""
import json
import platform

import numpy as np
import scipy
from scipy.sparse.linalg import expm_multiply

from .correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from .electronic_contract_v213 import (
    ElectronicModelSpaceV213, ElectronicOperatorProvenanceV213, ElectronicStateDescriptorV213,
)
from .molecular_frame_bridge_v300 import MolecularFrameBridgeV300
from .molecular_hamiltonian_v301 import (
    build_molecular_grid_hamiltonian_v301, build_molecular_tdvp_snapshot_v301,
)
from .molecular_midpoint_v302 import (
    molecular_midpoint_step_v302, run_molecular_midpoint_v302,
    serialize_molecular_trajectory_v302, replay_molecular_trajectory_v302,
)
from .synthetic_operator_provider_v21 import (
    SyntheticLinearOperatorConfigV21, SyntheticLinearOperatorProviderV21,
)


def free_packet_fixture_v302():
    provider = SyntheticLinearOperatorProviderV21(SyntheticLinearOperatorConfigV21(
        nstate=1, nq=1, mass=2.0, base_scale=0.0, derivative_scale=0.0,
    ))
    space = ElectronicModelSpaceV213(
        "scalar fixture", "fixed_general", (ElectronicStateDescriptorV213("scalar"),),
    )
    provider.provenance = ElectronicOperatorProvenanceV213(
        "zero-potential grid fixture", "1", space, "analytic zero potential",
        parameters={"mass": 2.0, "ndim": 1},
    ).validate()
    bridge = MolecularFrameBridgeV300(provider, [0.0])
    grid = build_molecular_grid_hamiltonian_v301(bridge, (np.linspace(-3, 3, 13),))
    state = CorrelatedGaussianSpinorStateV270(
        [[0.1]], [[0.3]], [[[1.1]]], [[[0.08]]], [[1.0+0j]],
    ).normalized()
    return state, bridge, grid


def run_v302_midpoint_campaign():
    """Separate time refinement from fixed-grid/basis projection error."""
    state, bridge, grid = free_packet_fixture_v302()
    initial = build_molecular_tdvp_snapshot_v301(state, bridge, grid).psi
    endpoints, rows, trajectories = [], [], []
    final_time = 0.16
    exact = expm_multiply(-1j*final_time*grid.matrix, initial.reshape(-1)).reshape(initial.shape)
    weighted_norm = lambda a: float(np.sqrt(np.sum(grid.quadrature.weights[:, None]*abs(a)**2)))
    for steps in (2, 4, 8):
        dt = final_time/steps
        trajectory = run_molecular_midpoint_v302(state, bridge, grid, dt=dt, steps=steps)
        psi = build_molecular_tdvp_snapshot_v301(trajectory[-1].after, bridge, grid).psi
        endpoints.append(psi)
        trajectories.append(trajectory)
        first = trajectory[0].receipt["before_observables"]
        last = trajectory[-1].receipt["after_observables"]
        rows.append({
            "dt": dt, "steps": steps,
            "norm_drift": abs(last["norm"]-first["norm"]),
            "energy_drift_hartree": abs(last["energy_hartree"]-first["energy_hartree"]),
            "grid_exact_wavefunction_error": weighted_norm(psi-exact),
            "maximum_scaled_midpoint_residual": max(s.receipt["scaled_midpoint_residual"] for s in trajectory),
            "maximum_svd_residual": max(s.receipt["midpoint_svd_residual"] for s in trajectory),
        })
    differences = [weighted_norm(endpoints[i]-endpoints[i+1]) for i in range(2)]
    order = float(np.log2(differences[0]/differences[1]))
    forward = trajectories[0][0]
    reverse = molecular_midpoint_step_v302(forward.after, bridge, grid, -forward.receipt["dt"])
    reverse_psi = build_molecular_tdvp_snapshot_v301(reverse.after, bridge, grid).psi
    reversal_error = weighted_norm(reverse_psi-initial)
    serialized = serialize_molecular_trajectory_v302(trajectories[0])
    replayed = replay_molecular_trajectory_v302(serialized, bridge, grid)
    replay_exact = serialize_molecular_trajectory_v302(replayed) == serialized
    gates = {
        "second_order_time_refinement": 1.8 < order < 2.2,
        "accepted_nonlinear_residuals": all(row["maximum_scaled_midpoint_residual"] <= 1 for row in rows),
        "svd_equation_residuals": all(row["maximum_svd_residual"] < 2e-8 for row in rows),
        "signed_step_reversal": reversal_error < 2e-8,
        "same_runtime_exact_replay": replay_exact,
        "fine_grid_norm_drift": rows[-1]["norm_drift"] < 1e-7,
        "energy_drift_refines": rows[0]["energy_drift_hartree"] > 8*rows[-1]["energy_drift_hartree"],
    }
    return {
        "schema": "gnd-molecular-midpoint-evidence-v0.30.2-dev",
        "source_kind": "analytic scalar zero-potential model",
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "passed": all(gates.values()), "gates": gates,
        "operator": grid.evidence(), "refinement": rows,
        "successive_wavefunction_differences": differences, "observed_time_order": order,
        "reversal_wavefunction_error": reversal_error,
        "trajectory": json.loads(serialized),
        "live_molecular_soc_validated": False, "dataset_certified": False,
        "scope": "Fixed basis, fixed finite grid; temporal self-convergence does not remove projection or spatial error.",
    }
