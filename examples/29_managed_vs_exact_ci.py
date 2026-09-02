from gaussian_dynamics.benchmark_suite import compare_managed_to_exact

out=compare_managed_to_exact(
    managed_dt=2e-4,
    exact_dt=5e-4,
    final_time=0.006,
    grid_n=32,
    spa_order=0,
)

print('Exact final adiabatic populations:  ',out['exact_populations'])
print('Managed graph-AIMS populations:    ',out['managed_populations'])
print('Population L2 error:                ',out['population_l2_error'])
print('Exact norm:                         ',out['exact_norm'])
print('Managed generalized norm:           ',out['managed_norm'])
print('Spawn events:                       ',out['spawn_events'])
