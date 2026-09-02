import numpy as np

from gaussian_dynamics.v18_benchmark import (
    V18AcceptanceThresholds,
    evaluate_v18_acceptance,
)


def test_v18_acceptance_contract_passes_good_synthetic_case():
    canonical={
        "wavefunction_projected":{
            "fidelity":0.99,
            "phase_aligned_l2_error":0.10,
            "nuclear_density_l2_error":0.04,
            "mean_error_l2":0.001,
            "covariance_error_frobenius":0.01,
        },
        "reduced_density_projected":{
            "density_frobenius_error":1e-4,
        },
        "reduced_density_target":{
            "density_frobenius_error":0.03,
        },
        "maximum_norm_drift":1e-6,
        "maximum_condition_number":1000.0,
        "complexity":{
            "sampled_audit_failures":0,
            "sentinel_dense_audits":2,
            "sentinel_pair_factorizations":100,
            "candidate_peak_memory_reduction_fraction":0.97,
        },
        "sentinel_audits":[
            {"passed":True},
            {"passed":True},
        ],
        "final_density_matrix":
            np.eye(2,dtype=complex)/2,
    }
    trajectory=[
        {"fidelity":0.99},
        {"fidelity":0.985},
    ]
    exact_overlap={
        "maximum_fidelity_drift":1e-12,
    }
    basis=[
        {"wavefunction_projected":{"phase_aligned_l2_error":0.20}},
        {"wavefunction_projected":{"phase_aligned_l2_error":0.18}},
        {"wavefunction_projected":{"phase_aligned_l2_error":0.16}},
        {"wavefunction_projected":{"phase_aligned_l2_error":0.14}},
    ]
    dt_axis={"observed_self_order":2.0}
    dt_self=[
        {"phase_aligned_l2_error":0.004},
        {"phase_aligned_l2_error":0.0005},
    ]
    edge=[
        {"wavefunction_projected":{"phase_aligned_l2_error":0.14}},
        {"wavefunction_projected":{"phase_aligned_l2_error":0.139}},
        {"wavefunction_projected":{"phase_aligned_l2_error":0.138}},
    ]
    growth=[
        {"wavefunction_projected":{"phase_aligned_l2_error":0.20}},
        {"wavefunction_projected":{"phase_aligned_l2_error":0.16}},
        {"wavefunction_projected":{"phase_aligned_l2_error":0.15}},
    ]
    v17={
        "complexity":{"audit_pair_factorizations":400},
        "final_density_matrix":
            np.eye(2,dtype=complex)/2,
    }

    out=evaluate_v18_acceptance(
        canonical,trajectory,exact_overlap,
        basis,dt_axis,dt_self,edge,growth,v17,
        V18AcceptanceThresholds(),
    )
    assert out["passed"]


def test_v18_acceptance_rejects_nonconvergent_basis_ladder():
    canonical={
        "wavefunction_projected":{
            "fidelity":0.99,
            "phase_aligned_l2_error":0.10,
            "nuclear_density_l2_error":0.04,
            "mean_error_l2":0.001,
            "covariance_error_frobenius":0.01,
        },
        "reduced_density_projected":{"density_frobenius_error":1e-4},
        "reduced_density_target":{"density_frobenius_error":0.03},
        "maximum_norm_drift":1e-6,
        "maximum_condition_number":1000.0,
        "complexity":{
            "sampled_audit_failures":0,
            "sentinel_dense_audits":2,
            "sentinel_pair_factorizations":100,
            "candidate_peak_memory_reduction_fraction":0.97,
        },
        "sentinel_audits":[{"passed":True},{"passed":True}],
        "final_density_matrix":np.eye(2,dtype=complex)/2,
    }
    trajectory=[{"fidelity":0.99}]
    exact_overlap={"maximum_fidelity_drift":0.0}
    basis=[
        {"wavefunction_projected":{"phase_aligned_l2_error":0.20}},
        {"wavefunction_projected":{"phase_aligned_l2_error":0.21}},
    ]
    dt_axis={"observed_self_order":2.0}
    dt_self=[
        {"phase_aligned_l2_error":0.004},
        {"phase_aligned_l2_error":0.0005},
    ]
    edge=[
        {"wavefunction_projected":{"phase_aligned_l2_error":0.14}},
        {"wavefunction_projected":{"phase_aligned_l2_error":0.13}},
    ]
    growth=[
        {"wavefunction_projected":{"phase_aligned_l2_error":0.20}},
        {"wavefunction_projected":{"phase_aligned_l2_error":0.14}},
    ]

    out=evaluate_v18_acceptance(
        canonical,trajectory,exact_overlap,
        basis,dt_axis,dt_self,edge,growth,None,
    )
    assert not out["passed"]
    assert not out["checks"]["basis_ladder_monotone"]
