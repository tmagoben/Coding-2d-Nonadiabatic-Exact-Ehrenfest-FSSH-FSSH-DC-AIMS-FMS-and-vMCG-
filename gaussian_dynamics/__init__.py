from .grids import uniform_grid, inner_product, normalize

__version__ = "0.27.0"
__development_version__ = "0.28.0.dev0"
from .gaussian import (
    frozen_gaussian,
    analytic_overlap,
    kinetic_on_gaussian,
    gaussian_moments,
)
from .exact import split_operator_step, run_split_operator
from .heller import (
    initial_heller_parameters,
    heller_wavefunction,
    run_thawed_gaussian,
    run_frozen_gaussian,
)
from .moving_basis import run_moving_gaussian_basis
from .variational import (
    pack_parameters,
    unpack_parameters,
    variational_wavefunction,
    tdvp_velocity,
    run_variational_dynamics,
)
from .spawning import (
    TrajectoryBasisFunction,
    spawn_child,
    coupled_basis_matrices,
)

__all__ = [
    "uniform_grid",
    "inner_product",
    "normalize",
    "frozen_gaussian",
    "analytic_overlap",
    "kinetic_on_gaussian",
    "gaussian_moments",
    "split_operator_step",
    "run_split_operator",
    "initial_heller_parameters",
    "heller_wavefunction",
    "run_thawed_gaussian",
    "run_frozen_gaussian",
    "run_moving_gaussian_basis",
    "pack_parameters",
    "unpack_parameters",
    "variational_wavefunction",
    "tdvp_velocity",
    "run_variational_dynamics",
    "TrajectoryBasisFunction",
    "spawn_child",
    "coupled_basis_matrices",
]

from .exact_multistate import run_multistate_exact, multistate_split_step
from .adiabatic import adiabatic_point, adiabatic_grid, adiabatic_hamiltonian_action
from .adiabatic_spawning import (
    AdiabaticTBF, coupling_indicator, spawn_child_energy_conserving,
    maybe_spawn, adiabatic_gaussian_basis_matrices
)

from .electronic_structure import (
    ElectronicStructurePoint, AnalyticAvoidedCrossingProvider,
    TabulatedElectronicStructureProvider, CachedProvider,
    project_cartesian_vector_to_coordinate, point_fingerprint
)
from .provider_dynamics import (
    ProviderTBF, velocity_verlet_tbf, coupling_indicator as provider_coupling_indicator,
    energy_conserving_child as provider_energy_conserving_child, provider_grid
)

from .ci2d import (
    LVC2DParameters, diabatic_potential_2d, adiabatic_energies_2d,
    adiabatic_gradients_2d, analytic_adiabatic_vectors, vector_nac_2d,
    branching_plane_vectors, circle_path, berry_line_integral,
    parallel_transport_real_state
)
from .gauge_transport import align_subspace, subspace_projector, projector_distance
from .gaussian_nd import (
    gaussian_nd, gaussian_nd_gradient, gaussian_nd_laplacian,
    kinetic_on_gaussian_nd, analytic_overlap_equal_width,
    gaussian_nd_time_derivative
)
from .heller_nd import initial_heller_nd, run_thawed_gaussian_nd
from .exact2d import run_exact_2d, split_operator_2d_step
from .spawned_basis_2d import (
    AdiabaticTBF2D, midpoint_grid_2d, adiabatic_fields_2d,
    basis_matrices_2d, nac_coupling_indicator, energy_conserving_child_nac,
    maybe_spawn_once, run_coupled_spawned_basis_2d
)

from .molecular_backend import (
    AMU_TO_ELECTRON_MASS, MolecularGeometry,
    CartesianElectronicStructurePoint, GeneralizedElectronicStructurePoint,
    LinearGeometryMap, GeneralizedCoordinateProvider,
    geometry_fingerprint
)
from .pyscf_backend_v05 import (
    PySCFSACASSCFConfig, PySCFSACASSCFBackend
)
from .benchmark_provider_nd import LVC2DGeneralizedProvider
from .local_gaussian_nd import (
    overlap_centroid_equal_width, gradient_matrix_element_equal_width,
    kinetic_matrix_element_equal_width, basis_time_matrix_element_equal_width,
    local_d2_matrix, LocalAdiabaticTBF, local_matrices
)
from .direct_dynamics_nd import (
    nac_indicator as backend_nac_indicator,
    energy_conserving_child as backend_energy_conserving_child,
    maybe_spawn_once as backend_maybe_spawn_once,
    run_backend_spawned_gaussians
)
from .backend_cache import DiskCachedGeneralizedProvider

from .state_tracking import (
    StateTrackingResult, maximum_overlap_assignment,
    transform_state_properties, reorder_and_phase_ci_roots,
    energy_degeneracy_clusters, subspace_overlap_singular_values
)
from .pyscf_wavefunction_overlap import (
    CASSCFWavefunctionSnapshot,
    embed_active_ci_with_doubly_occupied_core,
    casscf_state_overlap_matrix,
    correlated_orbital_cross_overlap
)
from .pyscf_tracked_backend_v06 import PySCFTrackedSACASSCFBackend
from .overlap_transport import (
    nearest_unitary, current_to_previous_procrustes,
    directional_nac_from_overlap, overlap_unitarity_defect,
    principal_angles
)
from .tracked_scan import (
    TrackedScanResult, run_tracked_scan, save_tracked_scan,
    TrackedScan1DProvider
)

# v0.7 graph-gauge API
from .gauge_graph import GaugeEdge, ElectronicGaugeGraph
from .graph_electronic import (
    adiabatic_hamiltonian_matrix,
    derivative_hamiltonian_matrices,
    rotate_operator,
    rotate_operator_field,
    rotate_coefficients,
    ElectronicOperatorNode,
    GraphElectronicRegistry,
)
from .pyscf_gauge_graph import (
    build_snapshot_gauge_graph,
    edge_overlap_diagnostics,
    tbf_centroid_edge_pairs,
)
from .graph_gaussian import (
    GraphGaussianTBF,
    pair_overlap_and_hamiltonian,
    build_static_graph_gaussian_matrices,
    generalized_cayley_step,
    generalized_norm,
)

# v0.8 time-dependent graph API
from .temporal_electronic import (
    hermitian_exponential, explicit_nac_step, overlap_strang_step,
    overlap_step_operator, electronic_fidelity
)
from .dynamic_gauge_graph import (
    ElectronicFramePoint, AnalyticCI2DFrameProvider, IncrementalElectronicGraph
)
from .moving_graph_gaussian import (
    nuclear_seed_basis_time_matrix, metric_compatible_basis_connection,
    basis_connection_residual, moving_basis_coefficient_step
)
from .dynamic_graph_aims import (
    DynamicGraphTBF, maybe_spawn_dynamic, run_dynamic_graph_aims
)
from .incremental_snapshot_graph import IncrementalSnapshotGaugeGraph

# v0.9 convergence-controlled graph-AIMS API
from .spa_matrix_elements import (
    SPAResult, real_saddle_point_equal_width, scalar_spa_matrix_element,
    graph_pair_spa_result, build_graph_gaussian_matrices_spa,
    spa1_correction_norm,
)
from .basis_management import (
    BasisConditioningReport, PruningResult, overlap_conditioning,
    project_coefficients_to_subset, prune_redundant_basis,
    canonical_orthogonalizer,
)
from .adaptive_spawning import CouplingExposureTracker, first_order_transfer_bound
from .convergence import (
    RefinementResult, phase_aligned_state_error, observed_order,
    scalar_refinement_study, vector_refinement_study, converged,
    save_convergence_json,
)
from .managed_graph_aims import (
    maybe_spawn_integrated, run_managed_graph_aims,
)
from .exact_benchmark import (
    localized_adiabatic_packet_2d, adiabatic_populations_from_diabatic,
    run_exact_ci_reference,
)
from .benchmark_suite import (
    run_managed_ci_case, compare_managed_to_exact,
    managed_timestep_refinement, spa_order_comparison,
)

from .initial_conditions import (
    GaussianWignerEnsemble, gaussian_wigner_covariances,
    sample_gaussian_wigner
)
from .benchmark_metrics import (
    ManagedRunMetrics, population_l2_error, population_l1_error,
    population_sum_error, norm_error, summarize_managed_run
)
from .error_budget import ErrorBudget, estimate_population_error_budget
from .benchmark_acceptance import (
    BenchmarkThresholds, BenchmarkAcceptance, evaluate_managed_benchmark
)
from .benchmark_campaign import (
    CIPassageConfig, run_managed_passage, run_exact_passage,
    run_managed_parameter_surface, run_exact_grid_timestep_surface,
    select_finest_exact_reference, campaign_settings_dict
)
from .ensemble_benchmark import (
    EnsembleStatistics, ensemble_statistics, run_ci_initial_condition_ensemble
)
from .campaign_io import save_campaign_json
from .electronic_observables import (
    reduced_electronic_density_from_vectors,
    reduced_electronic_density_graph,
    reduced_electronic_density_analytic_ci_diabatic,
    density_matrix_populations, density_matrix_purity,
    density_matrix_linear_entropy, density_matrix_von_neumann_entropy,
    exact_reduced_electronic_density_diabatic,
    rotate_density_to_frame, exact_reference_frame_density
)
from .reference_comparison import (
    compare_managed_exact_common_frame,
    compare_managed_exact_diabatic_density
)

from .release_benchmark import run_compact_v010_release_benchmark


from .gaussian_general import (
    validate_spd, gaussian_overlap_general, gaussian_cross_centroid,
    gaussian_cross_covariance, real_overlap_saddle_point,
    gradient_matrix_element_general, kinetic_matrix_element_general,
    basis_time_matrix_element_general, width_scaled
)
from .spa_matrix_elements_v11 import (
    GeneralSPAResult, graph_pair_spa_result_general,
    build_graph_gaussian_matrices_spa_general,
    spa1_correction_norm_general
)
from .moving_graph_gaussian_v11 import nuclear_seed_basis_time_matrix_general
from .optimized_spawning import (
    SpawnCandidate, classical_energy,
    energy_conserving_momentum_at_position,
    local_spa1_coupling_proxy, generate_spawn_candidates,
    select_spawn_children
)
from .managed_graph_aims_v11 import run_basis_complete_graph_aims
from .basis_completeness import (
    generation_histogram, width_determinants, width_diversity_ratio,
    overlap_spectrum_metrics, canonical_coefficient_weights,
    lineage_depth, basis_completeness_report
)
from .v11_benchmark import (
    V11AcceptanceThresholds, evaluate_v11_acceptance,
    run_v011_case, run_v011_release_benchmark
)


from .lvc_exact_gaussian import (
    ExactLVCPairResult,
    center_adiabatic_spinor,
    center_spinor_time_derivative,
    exact_lvc_potential_matrix_element,
    exact_lvc_pair_result,
    build_exact_lvc_gaussian_matrices,
    exact_lvc_basis_time_matrix,
)
from .moving_basis_v12 import (
    moving_basis_midpoint_cayley_step,
    fixed_basis_cayley_operator,
    endpoint_generalized_norm_error,
    phase_aligned_vector_error,
)
from .coherence_metrics import (
    offdiagonal_coherence,
    coherence_magnitude,
    wrapped_phase_difference,
    coherence_phase_error,
    coherence_magnitude_error,
    density_trace_distance,
    bloch_vector,
    bloch_vector_error,
)
from .coherent_lvc_dynamics_v12 import run_coherent_lvc_gaussians


from .local_diabatic_tbf_v12 import (
    LocalDiabaticTBF,
    from_adiabatic_guided_tbf,
    reset_to_instantaneous_adiabatic_spinor,
    parallel_transport_spinor_full_space,
)
from .spinor_complete_lvc_v12 import (
    build_nuclear_overlap_matrix,
    build_spinor_complete_lvc_matrices,
    build_spinor_complete_time_matrix,
    coefficients_matrix,
    flatten_coefficients,
    spinor_complete_reduced_density,
    spinor_complete_generalized_norm,
)
from .paired_basis_management_v12 import (
    PairedPruningResult,
    spinor_wavefunction_norm,
    project_spinor_coefficients_to_subset,
    prune_nuclear_gaussian_pairs,
)
from .spinor_complete_dynamics_v12 import (
    initialize_spinor_complete_coefficients,
    run_spinor_complete_lvc_gaussians,
)
from .born_huang_grid_v12 import (
    BornHuangGrid2D,
    build_born_huang_grid_2d,
    born_huang_basis_wavefunctions,
    born_huang_basis_fields,
    apply_spectral_kinetic_to_basis_fields,
    build_born_huang_matrices,
    born_huang_basis_time_matrix_grid,
    born_huang_basis_time_matrix,
    reconstruct_born_huang_wavefunction,
    born_huang_reduced_density,
)
from .born_huang_dynamics_v12 import (
    run_born_huang_projected_gaussians,
)


from .initial_projection_v12 import (
    InitialProjectionResult,
    make_shifted_initial_gaussian_bank,
    project_grid_wavefunction_to_spinor_complete_basis,
)


from .v12_benchmark import (
    V12AcceptanceThresholds,
    evaluate_v12_acceptance,
    load_v11_release_context,
    run_v012_release_benchmark,
)


from .residual_basis_v13 import (
    GaussianCandidate,
    ResidualCandidateScore,
    ResidualSelectionStep,
    ResidualBasisBuild,
    PreparedGaussianDictionary,
    cartesian_offsets_2d,
    generate_gaussian_dictionary,
    nuclear_overlap_matrix,
    normalized_grid_density,
    candidate_orthogonal_norm,
    residual_capture_gain,
    rank_residual_candidates,
    build_residual_greedy_basis,
    prepare_gaussian_dictionary,
    build_residual_greedy_basis_prepared,
)
from .tdse_defect_v13 import (
    TDSEDefect,
    DefectCandidateScore,
    DefectEnrichmentResult,
    reconstruct_spinor_complete_wavefunction,
    spinor_complete_coefficient_derivative,
    reconstruct_spinor_complete_time_derivative,
    apply_lvc_grid_hamiltonian,
    compute_tdse_defect,
    defect_candidate_capture,
    rank_defect_candidates,
    enrich_basis_from_tdse_defect,
)
from .v13_benchmark import (
    V13AcceptanceThresholds,
    load_v12_context,
    evaluate_v13_acceptance,
    run_v013_release_benchmark,
)


from .complexity_v14 import (
    AsymptoticComplexity,
    ComplexityLedger,
    asymptotic_complexity,
    dense_dimension_cost_proxy,
    pair_matrix_cost_proxy,
    candidate_ranking_cost_proxy,
)
from .residual_pruning_v14 import (
    LeaveOneOutScore,
    LowLossPruningResult,
    leave_one_out_projection_losses,
    prune_low_loss_gaussian_pair,
)
from .defect_candidates_v14 import (
    DynamicDefectCandidate,
    DynamicDefectScore,
    generate_energy_conserving_defect_candidates,
    rank_dynamic_defect_candidates_prepared,
)
from .adaptive_defect_dynamics_v14 import (
    AdaptiveDefectSettings,
    compute_tdse_defect_with_matrices,
    run_time_adaptive_defect_lvc_gaussians,
)


from .fast_lvc_matrices_v14 import (
    build_spinor_complete_lvc_matrices_symmetric,
    hermitian_pair_evaluation_count,
    ordered_pair_evaluation_count,
    pair_evaluation_reduction,
)


from .v14_benchmark import (
    V14AcceptanceThresholds,
    load_v13_context,
    evaluate_v14_acceptance,
    run_v014_release_benchmark,
)


from .pair_cache_v15 import (
    GaussianPairData,
    PairCacheStats,
    GaussianPairCache,
    build_cached_spinor_lvc_matrices,
    build_cached_spinor_time_matrix,
    expand_cached_spinor_lvc_matrices,
    subset_cached_spinor_lvc_matrices,
    v14_factorization_equivalent_for_sh,
    v14_factorization_equivalent_for_time,
    v15_factorization_count_for_snapshot,
)
from .complexity_v15 import (
    AsymptoticComplexityV15,
    ComplexityLedgerV15,
    asymptotic_complexity_v15,
    dense_cubic_units,
    canonical_pair_count,
    incremental_pair_count_for_add,
)
from .cost_aware_adaptation_v15 import (
    IncrementalCostEstimate,
    CostAwareCandidateScore,
    estimate_one_tbf_incremental_cost,
    rank_candidates_by_cost_aware_utility,
)
from .defect_candidates_v15 import (
    DynamicDefectCandidateV15,
    CachedDynamicDefectScore,
    generate_energy_conserving_defect_candidates_v15,
    rank_dynamic_defect_candidates_cached,
)
from .adaptive_defect_dynamics_v15 import (
    AdaptiveDefectSettingsV15,
    reduced_density_from_snuc,
    compute_tdse_defect_cached_v15,
    run_time_adaptive_cost_aware_lvc_gaussians,
)


from .v15_benchmark import (
    V15AcceptanceThresholds,
    load_v14_context,
    evaluate_v15_acceptance,
    run_v015_release_benchmark,
)


from .locality_graph_v16 import (
    conservative_position_overlap_bound,
    LocalityGraphSettings,
    LocalityGraphUpdate,
    PersistentGaussianLocalityGraph,
)
from .sparse_pair_matrices_v16 import (
    SparseSpinorMatrices,
    build_sparse_spinor_lvc_matrices,
    build_sparse_spinor_time_matrix,
    sparse_metric_compatible_connection,
    sparse_moving_basis_midpoint_cayley_step,
    sparse_generalized_norm,
    sparse_reduced_density,
    sparse_matrix_relative_difference,
)
from .electronic_cost_v16 import (
    ElectronicCostEstimate,
    UniformElectronicCostModel,
    GeometryCacheElectronicCostModel,
)
from .local_cost_aware_v16 import (
    LocalSparseCostEstimate,
    LocalSparseUtilityScore,
    predicted_local_degree,
    estimate_local_sparse_incremental_cost,
    rank_local_sparse_candidates,
)
from .sparse_complexity_v16 import (
    SparseComplexityModelV16,
    SparseComplexityLedgerV16,
    sparse_complexity_model_v16,
)
from .sparse_adaptive_dynamics_v16 import (
    SparseAdaptiveSettingsV16,
    compute_sparse_tdse_defect_v16,
    run_sparse_cost_aware_lvc_gaussians,
)


from .v16_benchmark import (
    V16AcceptanceThresholds,
    load_v15_context,
    run_sparse_scaling_benchmark,
    fit_sparse_scaling_exponents,
    electronic_cost_demo,
    evaluate_v16_acceptance,
    run_v016_release_benchmark,
)

from .sparse_pair_matrices_v16 import audit_sparse_lvc_matrices_against_dense


from .edge_importance_v17 import (
    EdgeImportanceSettingsV17,
    EdgeImportance,
    EdgeControlledGraphUpdateV17,
    safe_global_overlap_radius,
    pair_specific_overlap_upper_bound,
    exact_edge_importance,
    ErrorControlledGaussianLocalityGraphV17,
)
from .sparse_error_complexity_v17 import (
    SparseErrorComplexityLedgerV17,
)
from .error_controlled_sparse_dynamics_v17 import (
    ErrorControlledSparseSettingsV17,
    compute_sparse_tdse_defect_v16 as compute_error_controlled_tdse_defect_v17,
    run_error_controlled_sparse_lvc_gaussians,
)


from .sparse_error_budget_v17 import (
    monotone_nonincreasing,
    score_threshold_snapshot_sweep,
    local_score_budget_snapshot_sweep,
    summarize_snapshot_convergence,
)
from .v17_benchmark import (
    V17AcceptanceThresholds,
    load_v16_context,
    run_edge_controlled_scaling_benchmark,
    fit_v17_scaling_exponents,
    evaluate_v17_acceptance,
    run_v017_release_benchmark,
)


from .wavefunction_metrics_v18 import (
    grid_inner_product,
    grid_wavefunction_norm,
    normalize_grid_wavefunction,
    phase_align_wavefunction,
    phase_aligned_fidelity,
    phase_aligned_l2_error,
    nuclear_density,
    nuclear_density_l2_error,
    nuclear_density_total_variation,
    spatial_moments,
    moment_errors,
    gaussian_wavefunction_on_grid,
    compare_wavefunctions,
)
from .defect_candidates_v18 import (
    BatchedRankingDiagnosticsV18,
    rank_dynamic_defect_candidates_batched_v18,
)
from .sampled_sparse_audit_v18 import (
    SampledSparseAuditV18,
    sampled_omitted_edge_audit_v18,
)
from .convergence_complexity_v18 import (
    ConvergenceComplexityLedgerV18,
)
from .convergence_complete_dynamics_v18 import (
    ConvergenceCompleteSettingsV18,
    compute_sparse_tdse_defect_v16 as compute_convergence_tdse_defect_v18,
    run_convergence_complete_lvc_gaussians,
)


from .convergence_campaign_v18 import (
    ConvergenceCoordinatesV18,
    evaluate_convergence_run_v18,
    compare_snapshot_trajectory_v18,
    observed_order_from_dt,
    axis_sensitivity_summary,
    refinement_ladder_summary,
    successive_self_convergence_order,
)
from .v18_benchmark import (
    V18AcceptanceThresholds,
    load_v17_context,
    release_settings_v18,
    sampled_audit_scaling_v18,
    evaluate_v18_acceptance,
    run_v018_release_benchmark,
    assemble_v018_campaign_from_partials,
)

from .convergence_worker_v18 import run_coordinate_worker_v18


from .molecular_snapshot_v19 import (
    MolecularElectronicSnapshotV19,
    TrackedGeneralizedSnapshotV19,
)
from .analytic_molecular_backend_v19 import (
    AnalyticMolecularLVCConfigV19,
    AnalyticMolecularLVCBackendV19,
    default_diatomic_two_mode_map_v19,
)
from .state_tracking_v19 import (
    scalable_maximum_overlap_assignment_v19,
)
from .molecular_direct_provider_v19 import (
    BackendEvaluationPolicyV19,
    MolecularTrackingSettingsV19,
    MolecularProviderDiagnosticsV19,
    TrackedMolecularDirectProviderV19,
)
from .molecular_gauge_graph_v19 import (
    MolecularCentroidGraphV19,
    build_molecular_centroid_graph_v19,
)
from .molecular_direct_dynamics_v19 import (
    run_molecular_direct_dynamics_v19,
)
from .pyscf_molecular_bridge_v19 import (
    PySCFRawSnapshotBackendV19,
    pyscf_snapshot_overlap_engine_v19,
)
from .v19_benchmark import (
    V19AcceptanceThresholds,
    run_v019_release_benchmark,
)


from .indexed_molecular_provider_v20 import (
    IndexedCacheDiagnosticsV20,
    BufferedKDTreeIndexV20,
    IndexedTrackedMolecularDirectProviderV20,
)
from .sparse_molecular_matrices_v20 import (
    SparseMolecularTBFV20,
    MolecularSparseSettingsV20,
    MolecularPairDataV20,
    MolecularSparseUpdateV20,
    SparseMolecularMatricesV20,
    SparseMolecularEdgeGraphV20,
    molecular_pair_data_v20,
    build_sparse_molecular_matrices_v20,
    build_dense_molecular_reference_v20,
    dense_audit_sparse_molecular_v20,
)
from .sampled_molecular_audit_v20 import (
    SampledMolecularAuditV20,
    sampled_molecular_edge_audit_v20,
)
from .sparse_molecular_dynamics_v20 import (
    SparseMolecularDynamicsSettingsV20,
    run_sparse_molecular_dynamics_v20,
    run_dense_molecular_reference_dynamics_v20,
)
from .v20_benchmark import (
    V20AcceptanceThresholds,
    run_v020_release_benchmark,
)

from .electronic_operator_v21 import (
    ElectronicOperatorPointV21,
    ElectronicOperatorSnapshotV21,
    adiabatic_point_to_operator_v21,
    ElectronicOperatorProviderAdapterV21,
)
from .complex_gauge_v21 import (
    random_unitary_v21,
    PhaseMixingGaugeV21,
    transform_operator_point_v21,
    GaugeTransformedOperatorProviderV21,
)
from .subspace_tracking_v21 import (
    SubspaceTrackingResultV21,
    procrustes_subspace_alignment_v21,
    transform_subspace_operator_v21,
)
from .wilson_loop_v21 import (
    wilson_product_v21,
    gauge_transform_cycle_links_v21,
    sorted_wilson_eigenphases_v21,
)
from .synthetic_operator_provider_v21 import (
    SyntheticLinearOperatorConfigV21,
    SyntheticLinearOperatorProviderV21,
)
from .block_sparse_molecular_v21 import (
    BlockMolecularTBFV21,
    BlockSparseSettingsV21,
    BlockPairDataV21,
    BlockSparseUpdateV21,
    BlockSparseMatricesV21,
    BlockSparseMolecularGraphV21,
    block_pair_data_v21,
    build_block_sparse_matrices_v21,
    build_dense_block_reference_v21,
    block_diagonal_gauge_v21,
)
from .block_dynamics_v21 import (
    PrescribedBlockDynamicsSettingsV21,
    prescribed_linear_basis_v21,
    run_prescribed_block_dynamics_v21,
    gauge_block_matrices_v21,
    gauge_covariance_errors_v21,
    gauge_mapped_coefficient_error_v21,
)
from .v21_benchmark import (
    V21AcceptanceThresholds,
    run_v021_release_benchmark,
)

# v0.21.2 pre-SOC integration hardening
from .block_basis_lifecycle_v212 import (
    BlockPruneResultV212,
    insert_zero_block_v212,
    prune_block_projected_v212,
)
from .electronic_observables_v212 import (
    ElectronicObservableV212,
    build_electronic_observable_matrix_v212,
    observable_expectation_v212,
)
from .subspace_provider_v212 import (
    SubspaceTrackingSettingsV212,
    SubspaceProviderDiagnosticsV212,
    SubspaceAwareOperatorProviderV212,
)
from .self_consistent_block_v212 import (
    MeanFieldGuidanceSettingsV212,
    BlockMeanFieldGuidanceV212,
    SelfConsistentBlockSettingsV212,
    run_self_consistent_block_dynamics_v212,
)
from .complex_dtype_audit_v212 import (
    ComplexDtypeAuditResultV212,
    audit_pre_soc_complex_core_v212,
)
from .v212_benchmark import (
    V212AcceptanceThresholds,
    run_v0212_release_benchmark,
)

# v0.21.3 SOC-contract freeze
from .matrix_invariants_v213 import (
    MatrixInvariantTolerancesV213,
    scaled_matrix_residual_v213,
    hermiticity_residual_v213,
    antihermiticity_residual_v213,
    isometry_residual_v213,
    symmetry_residual_v213,
    require_residual_v213,
)
from .electronic_contract_v213 import (
    HARTREE_PER_WAVENUMBER_V213,
    wavenumber_to_hartree_v213,
    hartree_to_wavenumber_v213,
    ElectronicStateDescriptorV213,
    ElectronicModelSpaceV213,
    ElectronicOperatorProvenanceV213,
    compose_electronic_operator_v213,
    validate_electronic_contract_v213,
    ContractedElectronicOperatorProviderV213,
)
from .density_guidance_v213 import (
    DensityMatrixGuidanceSettingsV213,
    normalized_density_from_vector_v213,
    validate_guide_density_v213,
    density_force_v213,
    BlockDensityMatrixGuidanceV213,
)
from .self_consistent_block_v213 import (
    SelfConsistentBlockSettingsV213,
    run_self_consistent_block_dynamics_v213,
)
from .initial_projection_v213 import (
    GridProjectionResultV213,
    block_metric_fixed_frame_v213,
    transform_electronic_vector_to_local_frame_v213,
    initialize_separable_block_state_v213,
    project_grid_wavefunction_fixed_frame_v213,
)
from .complex_operator_cache_v213 import FixedFrameComplexOperatorCacheV213
from .v213_benchmark import (
    V213AcceptanceThresholds,
    run_v0213_release_benchmark,
)

# v0.21.4 differential-provider and deterministic-restart certification
from .provider_differential_audit_v214 import (
    ProviderDifferentialAuditSettingsV214,
    CoordinateDifferentialAuditV214,
    ProviderDifferentialAuditV214,
    audit_provider_differentials_v214,
    require_provider_differential_contract_v214,
)
from .checkpoint_restart_v214 import (
    SelfConsistentBlockSettingsV214,
    SelfConsistentBlockCheckpointV214,
    settings_fingerprint_v214,
    save_self_consistent_checkpoint_v214,
    load_self_consistent_checkpoint_v214,
    run_self_consistent_block_dynamics_v214,
)
from .zero_soc_rehearsal_v214 import (
    ZeroSOCRehearsalProviderV214,
    ZeroSOCEquivalenceReportV214,
    audit_zero_soc_equivalence_v214,
)
from .v214_benchmark import (
    V214AcceptanceThresholds,
    run_v0214_release_benchmark,
)

# v0.22.0 first physical analytic-SOC release
from .analytic_soc_models_v220 import (
    SOCOperatorComponentsV220,
    SingletTripletSOCConfigV220,
    AnalyticSingletTripletSOCProviderV220,
    DoubletSOCConfigV220,
    AnalyticDoubletSOCProviderV220,
    singlet_triplet_time_reversal_matrix_v220,
    kramers_time_reversal_matrix_v220,
    singlet_triplet_projectors_v220,
    doublet_root_projectors_v220,
)
from .physical_soc_validation_v220 import (
    PhysicalSOCAuditSettingsV220,
    ComponentDerivativeRowV221,
    PhysicalSOCAuditV220,
    KramersAuditV220,
    time_reversal_residual_v220,
    time_reversal_square_residual_v220,
    transform_time_reversal_matrix_v220,
    transform_projector_v220,
    projector_population_v220,
    audit_physical_soc_provider_v220,
    require_physical_soc_contract_v220,
    audit_kramers_degeneracy_v220,
)
from .spinor_exact_grid_v220 import (
    SpinorGridSettingsV220,
    normalize_spinor_grid_v220,
    initial_gaussian_spinor_v220,
    spinor_split_operator_step_v220,
    spinor_grid_energy_v220,
    spinor_grid_projector_population_v220,
    phase_aligned_spinor_grid_error_v220,
    run_spinor_exact_grid_v220,
)
from .v220_benchmark import (
    V220AcceptanceThresholds,
    run_v0220_release_benchmark,
)

# v0.22.1 corrective SOC admission and convergence hardening
from .soc_admission_v221 import (
    SOCSymmetryContractV221,
    SOCSymmetryAuditV221,
    audit_soc_symmetry_contract_v221,
    soc_symmetry_contract_from_provider_v221,
    require_soc_symmetry_contract_v221,
)
from .v221_benchmark import (
    V221AcceptanceThresholds,
    run_v0221_release_benchmark,
)

# v0.23.0 molecular-SOC backend admission and deterministic replay
from .molecular_soc_contract_v230 import (
    MolecularSOCCapabilitiesV230,
    MolecularSOCBackendIdentityV230,
    MolecularSOCValidationEvidenceV230,
    MolecularSOCAdmissionContractV230,
    provenance_with_molecular_soc_contract_v230,
    molecular_soc_contract_from_provider_v230,
    require_trajectory_ready_molecular_soc_v230,
)
from .molecular_soc_replay_v230 import (
    REPLAY_MANIFEST_NAME_V230,
    REPLAY_ARRAYS_NAME_V230,
    ReplayWavefunctionTokenV230,
    MolecularSOCReplayDatasetV230,
    capture_molecular_soc_replay_v230,
    load_molecular_soc_replay_v230,
    FileBackedMolecularSOCProviderV230,
)
from .molecular_soc_admission_v230 import (
    MolecularSOCAdmissionSettingsV230,
    MolecularSOCAdmissionAuditV230,
    audit_molecular_soc_provider_v230,
    require_molecular_soc_protocol_v230,
    require_real_molecular_soc_backend_v230,
)
from .pyscf_soc_bridge_v230 import (
    PySCFSOCRuntimeProbeV230,
    probe_pyscf_soc_runtime_v230,
    require_pyscf_soc_runtime_v230,
    PySCFMolecularSOCBridgeV230,
)
from .v230_benchmark import (
    V230AcceptanceThresholds,
    v230_reference_coordinates,
    build_v230_reference_replay,
    build_v230_doublet_reference_replay,
    run_v0230_release_benchmark,
)

# v0.23.1 raw-evidence dossiers and executable backend attestation
from .molecular_soc_evidence_v231 import (
    IndependentReferenceObservationV231,
    ConvergenceLadderObservationV231,
    FrameInvarianceObservationV231,
    TrackingSpecificationV231,
    DerivedEvidenceBundleV231,
)
from .molecular_soc_dossier_v231 import (
    DOSSIER_NAME_V231,
    RawArtifactRecordV231,
    CalculationReceiptV231,
    BackendRuntimeAttestationV231,
    MolecularSOCAdmissionDossierV231,
    write_raw_json_artifact_v231,
    write_molecular_soc_dossier_v231,
    load_molecular_soc_dossier_v231,
)
from .molecular_soc_admission_v231 import (
    MolecularSOCAdmissionAuditV231,
    audit_molecular_soc_provider_v231,
    require_molecular_soc_protocol_v231,
    require_external_molecular_soc_snapshot_v231,
    require_live_molecular_soc_backend_v231,
)
from .pyscf_soc_adapter_v231 import (
    PYSCF_NAC_CONVENTION_V231,
    PySCFMethodSpecificCapabilitiesV231,
    PySCFSOCAdapterProbeV231,
    probe_pyscf_soc_adapter_v231,
    require_pyscf_soc_adapter_v231,
    PySCFMethodSpecificSOCAdapterV231,
)
from .v231_benchmark import (
    V231AcceptanceThresholds,
    build_v231_admission_bundle,
    run_v0231_release_benchmark,
)

# v0.23.2 real PySCF runtime, physical overlap contract, and trusted admission
from .molecular_soc_replay_v230 import MolecularSOCReplayOverlapDiagnosticsV232
from .pyscf_nac_convention_v232 import (
    PYSCF_REQUIRED_VERSION_V232,
    PYSCF_NAC_UPSTREAM_DOCUMENTATION_V232,
    PYSCF_NAC_EMPIRICAL_MAPPING_V232,
    require_exact_pyscf_version_v232,
    pyscf_state_tuple_for_internal_dij_v232,
)
from .pyscf_runtime_v232 import (
    RUNTIME_SCHEMA_V232,
    PySCFRuntimeProbeV232,
    PySCFRuntimeFingerprintV232,
    PySCFRuntimeContextV232,
    PySCFRuntimeEvidenceV232,
    probe_pyscf_runtime_v232,
    require_pyscf_runtime_v232,
    guarded_pyscf_runtime_v232,
    build_pyscf_runtime_fingerprint_v232,
    run_pyscf_runtime_evidence_v232,
)
from .molecular_soc_runtime_v232 import (
    CONVERGENCE_VOCABULARY_V232,
    CONVERGENCE_METADATA_KEY_V232,
    CONVERGENCE_STAGES_V232,
    RUNTIME_PROBE_FORMAT_V232,
    ConvergenceMetadataV232,
    BackendMethodIdentityV232,
    RuntimeProbeRecordV232,
    ReceiptExecutionEvidenceV232,
    BackendArtifactValidationProofV232,
    BackendAdmissionPolicyV232,
    convergence_from_snapshot_v232,
    load_runtime_probe_v232,
)
from .molecular_soc_admission_v232 import (
    MolecularSOCRuntimeAdmissionAuditV232,
    audit_molecular_soc_provider_v232,
    require_external_molecular_soc_snapshot_v232,
    require_live_molecular_soc_backend_v232,
)
from .pyscf_soc_adapter_v232 import (
    PYSCF_NAC_CONVENTION_V232,
    PySCFMethodSpecificCapabilitiesV232,
    validate_pyscf_engine_contract_v232,
    PySCFMethodSpecificSOCAdapterV232,
)
from .v232_benchmark import (
    V232AcceptanceThresholds,
    run_v0232_release_benchmark,
)

# v0.23.3 finite-manifold transport and compatibility hardening
from .finite_manifold_transport_v233 import (
    OVERLAP_CONTRACT_ID_V233,
    TRANSPORT_CONTRACT_ID_V233,
    CONSUMER_OVERLAP_POLICY_V233,
    FiniteManifoldOverlapPolicyV233,
    FiniteManifoldTransportV233,
    ReciprocalTransportPairV233,
    analyze_finite_manifold_overlap_v233,
    certified_transport_from_overlap_v233,
    certify_reciprocal_transport_pair_v233,
)
from .nac_compatibility_v233 import (
    NAC_CONVENTION_SCHEMA_V233,
    INTERNAL_NAC_DEFINITION_V233,
    PYSCF_NAC_MAPPING_ID_V233,
    PYSCF_ETF_NAC_MAPPING_ID_V233,
    ANALYTIC_NAC_MAPPING_ID_V233,
    LEGACY_NAC_DISPOSITIONS_V233,
    DerivativeCouplingConventionV233,
    LegacyReplayMigrationAttestationV233,
    corrected_pyscf_nac_convention_v233,
    analytic_nac_convention_v233,
    derivative_coupling_convention_from_dict_v233,
    require_snapshot_nac_identity_v233,
)
from .molecular_soc_replay_v233 import (
    REPLAY_MANIFEST_NAME_V233,
    REPLAY_ARRAYS_NAME_V233,
    REPLAY_FORMAT_V233,
    REPLAY_FORMAT_VERSION_V233,
    ReplayWavefunctionTokenV233,
    MolecularSOCReplayDatasetV233,
    capture_molecular_soc_replay_v233,
    migrate_molecular_soc_replay_v230_to_v233,
    load_molecular_soc_replay_v233,
    FileBackedMolecularSOCProviderV233,
)
from .runtime_compatibility_v233 import (
    RuntimeCompatibilityProfileV233,
    RuntimeCompatibilityReportV233,
    release_locked_runtime_profile_v233,
    scientifically_compatible_runtime_profile_v233,
    assess_runtime_compatibility_v233,
)
from .manifold_transport_v233 import (
    ManifoldTransportPolicyV233,
    ManifoldBlockTransportV233,
    CompleteManifoldTransportAuditV233,
    audit_complete_manifold_transport_v233,
    require_complete_manifold_transport_v233,
)
from .provider_numerical_identity_v233 import (
    PROVIDER_IDENTITY_SCHEMA_V233,
    ProviderNumericalIdentityV233,
    build_provider_numerical_identity_v233,
    require_provider_numerical_identity_v233,
    run_convention_bound_dynamics_v233,
)
from .molecular_soc_convention_v233 import (
    MolecularSOCMatrixConventionV233,
    MolecularSOCConventionAuditV233,
    molecular_soc_convention_from_dict_v233,
    analytic_soc_convention_v233,
    audit_molecular_soc_convention_v233,
    require_molecular_soc_convention_v233,
    require_exact_molecular_soc_convention_v233,
)
from .v233_benchmark import (
    V233AcceptanceThresholds,
    run_v0233_release_benchmark,
)

# v0.24.0 fail-closed OpenMolcas RASSI-SO external snapshot intake
from .openmolcas_rassi_protocol_v240 import (
    OPENMOLCAS_PROTOCOL_SCHEMA_V240,
    OPENMOLCAS_EXPORT_SCHEMA_V240,
    OPENMOLCAS_MANIFEST_SCHEMA_V240,
    OpenMolcasRASSIProtocolV240,
    water_rassi_so_protocol_v240,
    openmolcas_protocol_from_dict_v240,
)
from .openmolcas_rassi_snapshot_v240 import (
    OPENMOLCAS_MANIFEST_NAME_V240,
    NATIVE_OPENMOLCAS_NUMERIC_CROSSCHECK_V240,
    OpenMolcasArtifactRecordV240,
    OpenMolcasBundleManifestV240,
    ParsedOpenMolcasRecordV240,
    ParsedOpenMolcasBundleV240,
    OpenMolcasRASSISnapshotParserV240,
    sha256_file_v240,
)
from .soc_derivative_evidence_v240 import (
    ExternalSOCDerivativePolicyV240,
    ExternalSOCDerivativeEvidenceV240,
    audit_external_soc_derivatives_v240,
    require_external_soc_derivatives_v240,
)
from .external_soc_validation_v240 import (
    ExternalSOCValidationAuditV240,
    audit_external_soc_validation_v240,
)
from .external_soc_admission_v240 import (
    ExternalSOCAdmissionPolicyV240,
    ExternalSOCAdmissionAuditV240,
    audit_external_soc_snapshot_v240,
    require_external_soc_snapshot_v240,
)
from .external_soc_dynamics_v240 import (
    FrozenSnapshotCheckpointV240,
    FrozenSnapshotDynamicsV240,
    preview_frozen_snapshot_dynamics_v240,
    run_admitted_external_soc_dynamics_v240,
)
from .v240_benchmark import (
    V240AcceptanceThresholds,
    build_v240_protocol_fixture,
    run_v0240_release_benchmark,
)

# v0.24.1 direct static PySCF BP-SOMF state-interaction SOC
from .pyscf_state_interaction_soc_v241 import (
    PYSCF_REQUIRED_VERSION_V241,
    PYSCF_BP_SOMF_PROVIDER_NAME_V241,
    PYSCF_BP_SOMF_PROVIDER_VERSION_V241,
    BP_SOMF_OPERATOR_FAMILY_V241,
    BP_SOMF_ONE_ELECTRON_INTEGRAL_V241,
    BP_SOMF_TWO_ELECTRON_INTEGRAL_V241,
    BP_SOMF_STATIC_LIMITATION_V241,
    PySCFStaticSOCProbeV241,
    probe_pyscf_static_soc_runtime_v241,
    require_pyscf_static_soc_runtime_v241,
    clebsch_gordan_twice_v241,
    SpinFreeRootV241,
    SpinMicrostateV241,
    complete_spin_microstates_v241,
    time_reversal_matrix_v241,
    root_projectors_v241,
    BPSOMFIntegralsV241,
    build_pyscf_bp_somf_integrals_v241,
    state_average_density_mo_from_pyscf_ci_v241,
    spin_ladder_ci_v241,
    wigner_reduced_transition_density_from_pyscf_ci_v241,
    StateInteractionSOCMatricesV241,
    assemble_state_interaction_soc_v241,
    PySCFStateInteractionSOCResultV241,
    PySCFStateInteractionSOCProviderV241,
    tampered_soc_result_v241,
)
from .pyscf_soc_runtime_v241 import (
    PYSCF_SOC_RUNTIME_SCHEMA_V241,
    OH_BOND_LENGTH_BOHR_V241,
    OH_ISOTOPE_MASSES_AMU_V241,
    PySCFStaticSOCAuditV241,
    PySCFStaticSOCRuntimeEvidenceV241,
    crosscheck_pyscf_somf_jk_v241,
    audit_pyscf_static_soc_v241,
    run_pyscf_oh_static_soc_evidence_v241,
    save_pyscf_oh_static_soc_evidence_v241,
)
from .v241_benchmark import (
    V241AcceptanceThresholds,
    run_v0241_release_benchmark,
)

# v0.24.2 connected-geometry direct-JK SOC differential preview
from .pyscf_differential_soc_v242 import (
    PYSCF_DIFFERENTIAL_SOC_SCHEMA_V242,
    PYSCF_DIRECT_JK_SOMF_STRATEGY_V242,
    PYSCF_DIFFERENTIAL_SOC_CAPABILITY_V242,
    OH_BOND_LENGTH_BOHR_V242,
    OH_ISOTOPE_MASSES_AMU_V242,
    OH_BOND_STEPS_BOHR_V242,
    build_pyscf_bp_somf_integrals_direct_jk_v242,
    PySCFSOCGeometrySnapshotV242,
    build_pyscf_soc_geometry_snapshot_v242,
    complete_multiplet_overlap_v242,
    phase_align_complete_multiplet_overlap_v242,
    TransportedSOCDerivativeV242,
    transported_soc_central_difference_v242,
    PySCFSOCDifferentialScanV242,
    run_pyscf_oh_bond_differential_soc_v242,
    PySCFSOCDifferentialAuditV242,
    audit_pyscf_oh_bond_differential_soc_v242,
    PySCFSOCDifferentialEvidenceV242,
    run_pyscf_oh_bond_differential_evidence_v242,
    save_pyscf_oh_bond_differential_evidence_v242,
)
from .v242_benchmark import (
    V242AcceptanceThresholds,
    run_v0242_release_benchmark,
)

# v0.25.0 restricted time-dependent-variational SOC dynamics
from .variational_soc_dynamics_v250 import (
    VARIATIONAL_SOC_SCHEMA_V250,
    RESTRICTED_TDVP_ANSATZ_V250,
    RESTRICTED_NUCLEAR_INTEGRATOR_V250,
    GENERAL_TDVP_INTEGRATOR_V250,
    ELECTRONIC_INTEGRATOR_V250,
    POLAR_ALGORITHM_V250,
    VariationalSOCIntegratorSettingsV250,
    CanonicalVariationalSOCStateV250,
    SymmetricVariationalSOCStepV250,
    V250_TRAJECTORY_CLAIMS,
    SymmetricVariationalSOCTrajectoryV250,
    symmetric_variational_soc_step_v250,
    run_symmetric_variational_soc_dynamics_v250,
    reverse_variational_soc_trajectory_v250,
)
from .variational_soc_validation_v250 import (
    VARIATIONAL_SOC_VALIDATION_SCHEMA_V250,
    V250_CONVERGENCE_DT_AU,
    V250_CONVERGENCE_FINAL_TIME_AU,
    VariationalSOCValidationAuditV250,
    VariationalSOCValidationEvidenceV250,
    run_variational_soc_validation_evidence_v250,
    save_variational_soc_validation_evidence_v250,
)
from .v250_benchmark import (
    V250AcceptanceThresholds,
    run_v0250_release_benchmark,
)

# v0.25.1 frozen-width multi-Gaussian TDVP metric layer
from .multigaussian_tdvp_v251 import (
    MULTIGAUSSIAN_TDVP_SCHEMA_V251,
    MULTIGAUSSIAN_TDVP_ANSATZ_V251,
    VARIATIONAL_PRINCIPLE_V251,
    VARIATIONAL_INTEGRATOR_V251,
    VARIATIONAL_METRIC_SOLVER_V251,
    POTENTIAL_CONTRACT_V251,
    VariationalMetricSettingsV251,
    QuadraticSpinHamiltonianV251,
    FrozenGaussianSpinorStateV251,
    MetricSolveReceiptV251,
    VariationalMetricSystemV251,
    ImplicitMidpointTDVPStepV251,
    FrozenWidthMultiGaussianTrajectoryV251,
    V251_TDVP_CLAIMS,
    quadratic_spin_hamiltonian_from_provider_v251,
    pack_variational_parameters_v251,
    state_from_variational_parameters_v251,
    build_frozen_gaussian_spinor_matrices_v251,
    variational_energy_v251,
    solve_variational_metric_v251,
    build_variational_metric_system_v251,
    implicit_midpoint_tdvp_step_v251,
    run_frozen_width_multigaussian_tdvp_v251,
    reverse_frozen_width_multigaussian_tdvp_v251,
)
from .multigaussian_tdvp_validation_v251 import (
    MULTIGAUSSIAN_TDVP_VALIDATION_SCHEMA_V251,
    V251_CONVERGENCE_DT_AU,
    V251_CONVERGENCE_FINAL_TIME_AU,
    MultiGaussianTDVPValidationAuditV251,
    MultiGaussianTDVPValidationEvidenceV251,
    run_multigaussian_tdvp_validation_evidence_v251,
    save_multigaussian_tdvp_validation_evidence_v251,
)
from .v251_benchmark import (
    V251AcceptanceThresholds,
    run_v0251_release_benchmark,
)

# v0.25.2 adaptive log-width/quadratic-chirp multi-Gaussian TDVP
from .adaptive_multigaussian_tdvp_v252 import (
    ADAPTIVE_MULTIGAUSSIAN_TDVP_SCHEMA_V252,
    ADAPTIVE_MULTIGAUSSIAN_TDVP_ANSATZ_V252,
    VARIATIONAL_PRINCIPLE_V252,
    VARIATIONAL_INTEGRATOR_V252,
    VARIATIONAL_METRIC_SOLVER_V252,
    WIDTH_COORDINATES_V252,
    POTENTIAL_CONTRACT_V252,
    QuadraticSpinHamiltonianV252,
    AdaptiveVariationalSettingsV252,
    ThawedGaussianSpinorStateV252,
    AdaptiveVariationalMetricSystemV252,
    AdaptiveImplicitMidpointTDVPStepV252,
    AdaptiveWidthMultiGaussianTrajectoryV252,
    V252_TDVP_CLAIMS,
    quadratic_spin_hamiltonian_from_provider_v252,
    pack_adaptive_variational_parameters_v252,
    state_from_adaptive_variational_parameters_v252,
    build_adaptive_gaussian_spinor_matrices_v252,
    adaptive_variational_energy_v252,
    build_adaptive_variational_metric_system_v252,
    adaptive_implicit_midpoint_tdvp_step_v252,
    run_adaptive_width_multigaussian_tdvp_v252,
    reverse_adaptive_width_multigaussian_tdvp_v252,
)
from .adaptive_multigaussian_tdvp_validation_v252 import (
    ADAPTIVE_MULTIGAUSSIAN_VALIDATION_SCHEMA_V252,
    V252_CONVERGENCE_DT_AU,
    V252_CONVERGENCE_FINAL_TIME_AU,
    AdaptiveMultiGaussianValidationAuditV252,
    AdaptiveMultiGaussianValidationEvidenceV252,
    run_adaptive_multigaussian_validation_evidence_v252,
    save_adaptive_multigaussian_validation_evidence_v252,
)
from .v252_benchmark import (
    V252AcceptanceThresholds,
    run_v0252_release_benchmark,
)

# v0.25.3 controlled residual-driven adaptive-basis lifecycle
from .controlled_basis_adaptation_v253 import (
    CONTROLLED_BASIS_SCHEMA_V253,
    CONTROLLED_BASIS_EVENT_SCHEMA_V253,
    SPAWN_SCORE_V253,
    PROJECTION_POLICY_V253,
    EVENT_ORDER_V253,
    POTENTIAL_CONTRACT_V253,
    ControlledBasisSettingsV253,
    SpawnCandidateV253,
    BasisProjectionReceiptV253,
    SpawnCandidateEvaluationV253,
    BasisLifecycleEventV253,
    ControlledMetricSystemV253,
    CoefficientActivationStepV253,
    ControlledBasisStepV253,
    ControlledBasisTrajectoryV253,
    V253_CONTROLLED_BASIS_CLAIMS,
    generate_spawn_candidates_v253,
    project_adaptive_state_v253,
    evaluate_spawn_candidate_v253,
    adapt_basis_once_v253,
    build_controlled_metric_system_v253,
    coefficient_activation_implicit_step_v253,
    controlled_tdvp_step_v253,
    run_controlled_basis_dynamics_v253,
)
from .controlled_basis_validation_v253 import (
    CONTROLLED_BASIS_VALIDATION_SCHEMA_V253,
    ControlledBasisValidationAuditV253,
    ControlledBasisValidationEvidenceV253,
    run_controlled_basis_validation_evidence_v253,
    save_controlled_basis_validation_evidence_v253,
)
from .v253_benchmark import (
    V253AcceptanceThresholds,
    run_v0253_release_benchmark,
    save_v0253_release_benchmark,
)

# v0.26.0 reference-first multidimensional CI+SOC dynamics
from .multidimensional_soc_v260 import (
    MULTIDIMENSIONAL_SOC_MODEL_SCHEMA_V260,
    EXACT_GRID_SCHEMA_V260,
    QUADRATIC_CONVENTION_V260,
    KINETIC_CONVENTION_V260,
    ELECTRONIC_FRAME_CONVENTION_V260,
    GRID_INTEGRATOR_V260,
    QuadraticSpinHamiltonianNDV260,
    UniformGrid2DV260,
    ExactGridSettingsV260,
    ExactGridTrajectoryV260,
    V260_EXACT_GRID_CLAIMS,
    two_state_ci_soc_model_v260,
    kramers_doublet_ci_soc_model_v260,
    singlet_triplet_ci_soc_model_v260,
    normalize_spinor_grid_v260,
    initial_gaussian_spinor_2d_v260,
    exact_grid_split_step_v260,
    exact_grid_norm_v260,
    exact_grid_overlap_v260,
    exact_grid_boundary_probability_v260,
    exact_grid_reduced_density_v260,
    exact_grid_energy_v260,
    phase_aligned_grid_error_v260,
    run_exact_grid_ci_soc_v260,
)
from .multidimensional_gaussian_tdvp_v260 import (
    MULTIDIMENSIONAL_TDVP_SCHEMA_V260,
    MULTIDIMENSIONAL_TDVP_ANSATZ_V260,
    VARIATIONAL_PRINCIPLE_V260,
    VARIATIONAL_INTEGRATOR_V260,
    VARIATIONAL_METRIC_SOLVER_V260,
    WIDTH_CONVENTION_V260,
    MultidimensionalVariationalSettingsV260,
    DiagonalGaussianSpinorStateV260,
    MultidimensionalMetricSystemV260,
    MultidimensionalImplicitMidpointStepV260,
    MultidimensionalTDVPTrajectoryV260,
    V260_MULTIDIMENSIONAL_TDVP_CLAIMS,
    pack_multidimensional_parameters_v260,
    state_from_multidimensional_parameters_v260,
    active_parameter_indices_v260,
    build_multidimensional_gaussian_matrices_v260,
    multidimensional_variational_energy_v260,
    multidimensional_reduced_density_v260,
    build_multidimensional_metric_system_v260,
    multidimensional_implicit_midpoint_step_v260,
    run_multidimensional_tdvp_v260,
    evaluate_multidimensional_state_v260,
    multidimensional_state_on_grid_v260,
    residual_coupling_at_geometry_v260,
)
from .multidimensional_basis_adaptation_v260 import (
    MULTIDIMENSIONAL_BASIS_SCHEMA_V260,
    MULTIDIMENSIONAL_EVENT_SCHEMA_V260,
    SPAWN_SCORE_V260,
    PROJECTION_POLICY_V260,
    EVENT_ORDER_V260,
    NEWBORN_ACTIVATION_V260,
    ControlledMultidimensionalBasisSettingsV260,
    MultidimensionalSpawnCandidateV260,
    BasisProjectionReceiptV260,
    SpawnCandidateEvaluationV260,
    MultidimensionalBasisEventV260,
    ControlledMultidimensionalStepV260,
    ControlledMultidimensionalTrajectoryV260,
    V260_MULTIDIMENSIONAL_BASIS_CLAIMS,
    generate_multidimensional_spawn_candidates_v260,
    project_multidimensional_state_v260,
    evaluate_multidimensional_spawn_candidate_v260,
    metric_compatible_activation_mask_v260,
    adapt_multidimensional_basis_once_v260,
    run_controlled_multidimensional_dynamics_v260,
)
from .multidimensional_validation_v260 import (
    MULTIDIMENSIONAL_VALIDATION_SCHEMA_V260,
    MultidimensionalValidationEvidenceV260,
    run_multidimensional_validation_evidence_v260,
    save_multidimensional_validation_evidence_v260,
)
from .v260_benchmark import (
    V260AcceptanceThresholds,
    run_v0260_release_benchmark,
    save_v0260_release_benchmark,
)

# v0.27.0 full complex-symmetric correlated Gaussian widths and chirps
from .correlated_gaussian_tdvp_v270 import (
    CORRELATED_TDVP_SCHEMA_V270,
    CORRELATED_TDVP_ANSATZ_V270,
    VARIATIONAL_PRINCIPLE_V270,
    VARIATIONAL_INTEGRATOR_V270,
    VARIATIONAL_METRIC_SOLVER_V270,
    WIDTH_CONVENTION_V270,
    CorrelatedVariationalSettingsV270,
    CorrelatedGaussianSpinorStateV270,
    CorrelatedMetricSystemV270,
    CorrelatedImplicitMidpointStepV270,
    CorrelatedTDVPTrajectoryV270,
    V270_CORRELATED_TDVP_CLAIMS,
    symmetric_pairs_v270,
    symmetric_size_v270,
    symmetric_basis_v270,
    svec_v270,
    smat_v270,
    log_spd_v270,
    exp_symmetric_v270,
    exp_frechet_symmetric_v270,
    rotate_symmetric_v270,
    cross_correlated_gaussian_data_v270,
    correlated_moment_table_v270,
    integrate_correlated_polynomial_v270,
    pack_correlated_parameters_v270,
    state_from_correlated_parameters_v270,
    active_correlated_parameter_indices_v270,
    build_correlated_gaussian_matrices_v270,
    correlated_variational_energy_v270,
    correlated_reduced_density_v270,
    build_correlated_metric_system_v270,
    correlated_implicit_midpoint_step_v270,
    run_correlated_tdvp_v270,
    evaluate_correlated_state_v270,
    correlated_state_on_grid_v270,
    residual_coupling_correlated_v270,
    rotate_correlated_velocity_v270,
    gauge_correlated_velocity_v270,
    permute_correlated_velocity_v270,
)
from .correlated_basis_adaptation_v270 import (
    CORRELATED_BASIS_SCHEMA_V270,
    CORRELATED_BASIS_EVENT_SCHEMA_V270,
    SPAWN_SCORE_V270,
    SPAWN_DIRECTIONS_V270,
    PROJECTION_POLICY_V270,
    EVENT_ORDER_V270,
    NEWBORN_ACTIVATION_V270,
    ControlledCorrelatedBasisSettingsV270,
    CorrelatedSpawnCandidateV270,
    CorrelatedBasisProjectionReceiptV270,
    CorrelatedSpawnCandidateEvaluationV270,
    CorrelatedBasisEventV270,
    ControlledCorrelatedStepV270,
    ControlledCorrelatedTrajectoryV270,
    V270_MULTIDIMENSIONAL_BASIS_CLAIMS,
    generate_correlated_spawn_candidates_v270,
    project_correlated_state_v270,
    evaluate_correlated_spawn_candidate_v270,
    metric_compatible_activation_mask_v270,
    adapt_correlated_basis_once_v270,
    run_controlled_correlated_dynamics_v270,
)
from .correlated_validation_v270 import (
    CORRELATED_VALIDATION_SCHEMA_V270,
    CorrelatedValidationEvidenceV270,
    run_correlated_validation_evidence_v270,
    save_correlated_validation_evidence_v270,
)
from .v270_benchmark import (
    V270AcceptanceThresholds,
    run_v0270_release_benchmark,
    save_v0270_release_benchmark,
)


# v0.28.0 development: flat coordinate-dependent electronic frames
from .moving_frame_v280 import (
    CLAIM_BOUNDARY_V280, CONNECTION_CONVENTION_V280, TRANSPORT_CONVENTION_V280,
    FlatMovingFrameV280, MovingFrameCorrelatedStateV280, MovingFrameImplicitMidpointStepV280,
    MovingFrameBasisEventV280, adapt_moving_frame_basis_once_v280,
    evaluate_moving_physical_v280, evaluate_moving_section_v280, fixed_to_moving_state_v280,
    moving_frame_hamiltonian_v280, moving_frame_implicit_midpoint_step_v280,
    moving_frame_velocity_v280, moving_to_fixed_state_v280, reference_wavefunction_error_v280,
    require_flat_moving_frame_v280,
)
from .moving_frame_validation_v280 import (
    LatticeGaugeOracleV280, build_lattice_gauge_oracle_v280,
    finite_difference_connection_residual_v280, finite_difference_curvature_residual_v280,
    lattice_action_covariance_v280, lattice_propagation_covariance_v280,
)
from .moving_frame_evidence_v280 import (
    MOVING_FRAME_EVIDENCE_SCHEMA_V280, MovingFrameEvidenceV280, run_moving_frame_evidence_v280,
)
