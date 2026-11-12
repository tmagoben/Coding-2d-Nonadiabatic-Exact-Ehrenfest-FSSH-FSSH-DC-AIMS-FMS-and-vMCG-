"""Recompute v0.30.1 analytic-model evidence (no ab-initio trajectory claim)."""
from pathlib import Path
import argparse
import json
import platform
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import scipy
from scipy.linalg import block_diag
from scipy.sparse.linalg import expm_multiply

from gaussian_dynamics.analytic_soc_models_v220 import AnalyticSingletTripletSOCProviderV220
from gaussian_dynamics.complex_gauge_v21 import (
    GaugeTransformedOperatorProviderV21, PhaseMixingGaugeV21, random_unitary_v21,
)
from gaussian_dynamics.correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from gaussian_dynamics.molecular_frame_bridge_v300 import MolecularFrameBridgeV300
from gaussian_dynamics.molecular_hamiltonian_v301 import (
    build_molecular_grid_hamiltonian_v301, build_molecular_tdvp_snapshot_v301,
)


def run():
    provider = AnalyticSingletTripletSOCProviderV220()
    gauge = PhaseMixingGaugeV21(
        random_unitary_v21(4, seed=301), np.array([[0.31], [-0.22], [0.17], [0.28]]),
        np.array([0.1, 0.2, -0.1, -0.3]),
    )
    bridge = MolecularFrameBridgeV300(provider, [0.0])
    moving_bridge = MolecularFrameBridgeV300(GaugeTransformedOperatorProviderV21(provider, gauge), [0.0])
    axes = (np.linspace(-0.8, 0.8, 9),)
    fixed = build_molecular_grid_hamiltonian_v301(bridge, axes)
    moving = build_molecular_grid_hamiltonian_v301(moving_bridge, axes)
    transform = block_diag(*[gauge.matrix(q) for q in fixed.quadrature.points])
    covariance = float(np.linalg.norm(
        moving.matrix.toarray() - transform.conj().T @ fixed.matrix.toarray() @ transform,
    ))
    state = CorrelatedGaussianSpinorStateV270(
        [[0.13]], [[0.2]], [[[1.1]]], [[[0.08]]], [[0.8, 0.3j, 0.1, -0.2j]],
    ).normalized()
    evaluation = build_molecular_tdvp_snapshot_v301(state, bridge, fixed)
    psi = evaluation.psi.reshape(-1)
    psi = psi / np.linalg.norm(psi)
    propagated = expm_multiply(-0.7j * fixed.matrix, psi)
    moved = expm_multiply(-0.7j * moving.matrix, transform.conj().T @ psi)
    propagation_covariance = float(np.linalg.norm(moved - transform.conj().T @ propagated))
    norm_drift = float(abs(np.vdot(propagated, propagated) - 1.0))
    energy_drift = float(abs(np.vdot(propagated, fixed.matrix @ propagated) - np.vdot(psi, fixed.matrix @ psi)))
    mass = float(provider.evaluate([0.0]).mass_matrix_q_au[0, 0])
    refinement = []
    for spacing in (0.4, 0.2, 0.1):
        axis = np.arange(-8 + spacing, 8 - spacing/2, spacing)
        grid = build_molecular_grid_hamiltonian_v301(bridge, (axis,))
        gaussian = np.exp(-0.5*axis**2 + 0.3j*axis)
        vector = np.zeros((len(axis), 4), complex)
        vector[:, 0] = gaussian
        kinetic = (grid.kinetic @ vector.reshape(-1)).reshape(vector.shape)[:, 0]
        exact = (1 - (-axis + 0.3j)**2) * gaussian / (2*mass)
        error = float(np.linalg.norm(kinetic-exact) / np.linalg.norm(exact))
        refinement.append({"spacing_bohr": spacing, "relative_kinetic_error": error})
    errors = np.array([row["relative_kinetic_error"] for row in refinement])
    orders = np.log2(errors[:-1] / errors[1:])
    evidence = moving.evidence()
    gates = {
        "hermitian_full_hamiltonian": evidence["hermiticity_residual"] < 1e-13,
        "coordinate_dependent_gauge_covariance": covariance < 1e-11,
        "second_order_kinetic_refinement": bool(np.all((orders > 1.9) & (orders < 2.1))),
        "discrete_grid_propagation_gauge_covariance": propagation_covariance < 1e-11,
        "discrete_grid_norm_preservation": norm_drift < 1e-12,
        "discrete_grid_energy_preservation": energy_drift < 1e-12,
        "tdvp_svd_equation_residual": evaluation.snapshot.solve_receipt.residual_relative < 2e-8,
        "tdvp_null_rhs_compatibility": evaluation.snapshot.solve_receipt.null_rhs_relative < 2e-8,
    }
    return {
        "schema": "gnd-molecular-hamiltonian-evidence-v0.30.1-dev",
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "source_kind": "analytic singlet/triplet SOC model",
        "passed": bool(all(gates.values())), "gates": gates,
        "operator": evidence,
        "full_hamiltonian_gauge_covariance_error": covariance,
        "grid_propagation_gauge_covariance_error": propagation_covariance,
        "grid_norm_drift": norm_drift, "grid_energy_drift_hartree": energy_drift,
        "kinetic_refinement": refinement, "observed_orders": orders.tolist(),
        "tdvp_metric_rank": evaluation.snapshot.solve_receipt.rank,
        "tdvp_metric_condition": evaluation.snapshot.solve_receipt.condition_number,
        "tdvp_svd_residual": evaluation.snapshot.solve_receipt.residual_relative,
        "tdvp_null_rhs_relative": evaluation.snapshot.solve_receipt.null_rhs_relative,
        "tdse_defect_norm": float(np.linalg.norm(evaluation.snapshot.residual_weighted)),
        "gaussian_time_trajectory_validated": False,
        "live_molecular_soc_validated": False,
        "dataset_certified": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/v0301_molecular_hamiltonian_evidence.json")
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    raise SystemExit(0 if result["passed"] else 1)
